"""Traveling Telescope Problem (TTP) solver

The model implements the MILP of Handley et al. 2024 (arXiv:2310.18497).
"""

# Standard library imports
import logging
import time
from itertools import product

# Third-party imports
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.table import QTable
import gurobipy as gp
from gurobipy import GRB

logs = logging.getLogger(__name__)

#: Columns required on the ``requests`` QTable passed to :class:`TTPModel`.
REQUIRED_COLUMNS = (
    "unique_id",       # str, primary key
    "coord",           # SkyCoord column (ICRS)
    "first_available", # Time column
    "last_available",  # Time column
    "t_visit",         # Quantity column with time units
    "n_intra_max",     # int column
    "tau_intra",       # Quantity column with time units
    "priority",        # float column
)


def _validate_requests(requests):
    """Type-check the input ``requests`` QTable."""
    if not isinstance(requests, QTable):
        raise TypeError(
            "`requests` must be an astropy.table.QTable; got "
            f"{type(requests).__name__}"
        )
    missing = [col for col in REQUIRED_COLUMNS if col not in requests.colnames]
    if missing:
        raise ValueError(f"`requests` missing required columns: {missing}")
    if len(requests) == 0:
        return
    if not isinstance(requests["coord"], SkyCoord):
        raise TypeError("`coord` column must be an astropy SkyCoord")
    if not isinstance(requests["first_available"], Time):
        raise TypeError("`first_available` column must be an astropy Time")
    if not isinstance(requests["last_available"], Time):
        raise TypeError("`last_available` column must be an astropy Time")
    for col in ("t_visit", "tau_intra"):
        unit = requests[col].unit
        if unit is None or not unit.is_equivalent(u.s):
            raise TypeError(
                f"`{col}` column must be a Quantity with time units; "
                f"got unit={unit}"
            )


class TTPModel:
    """MILP solver for the Traveling Telescope Problem (Handley+ 2024).

    Args:
        requests (astropy.table.QTable): one row per request, with native
            astropy-typed columns:

            unique_id        str                            primary key
            coord            SkyCoord column                ICRS
            first_available  Time column                    earliest start of accessibility window
            last_available   Time column                    latest end of accessibility window
            t_visit          Quantity column (time)         per-visit duration
            n_intra_max      int column                     max visits per night
            tau_intra        Quantity column (time)         min spacing between visits within a night
            priority         float column                   objective weight; higher = more important

        night_start (astropy.time.Time): start of the observing interval.
        night_end (astropy.time.Time): end of the observing interval.

    Keyword Args:
        slew_fn (callable): ``slew_fn(coord_a, coord_b, time) -> np.ndarray``
            of slew minutes on the cartesian grid ``(pairs, time)``. ``coord_a``
            and ``coord_b`` are pair-aligned 1-D ``SkyCoord`` arrays of length
            ``P`` (pair ``k`` is ``(coord_a[k], coord_b[k])``); ``time`` is a
            1-D ``Time`` of length ``T``. The return has shape ``(P, T)`` and
            its ``[k, t]`` entry is the slew time from ``coord_a[k]`` to
            ``coord_b[k]`` evaluated at ``time[t]``. ``(P, T)`` matches
            astroplan's ``Observer.altaz(..., grid_times_targets=True)``
            convention.
        n_slots (int): number of TTP slew slots ``M`` (Handley+ 2024 §2.2).
            ``n_slots=1`` is the recommended default.
        slew_sample_cadence_min (int): max spacing in minutes at which to
            sample arcs within a slot.

    The input QTable ``self.requests`` is the single source of truth for
    the per-request inputs (native astropy / numpy types). Downstream
    relational state lives in three DataFrames built by the solver:

    * ``self.nodes`` -- one row per MILP node (:meth:`build_nodes`).
      Scalar-only dtypes (no object columns) so it is round-trip-safe
      through ``to_hdf``/``to_csv``.
    * ``self.arcs`` -- precomputed arc catalog ``(i, j, m)`` (:meth:`build_arcs`).
    * ``self.schedule`` -- post-solve output parallel to ``nodes``
      (:meth:`build_schedule`); ``None`` if Gurobi finds no incumbent
      within the time limit.

    Notes:
        Internal naming is aligned with Handley+ 2024 (``N``, ``M``, ``Yi``,
        ``Xijm``, ``arcs``) and AstroQ vocabulary (``t_visit``, ``tau_intra``).

    Usage:
        tm = TTPModel(...)
        tm.build_nodes() # builds nodes from requests
        tm.build_arcs()  # builds arcs from nodes
        tm.build_model() # builds the MILP

        # set the Gurobi parameters
        tm.model.params.TimeLimit = 300  # set the time limit
        tm.model.params.MIPGap = 0.05  # set the MIP gap
        tm.model.params.OutputFlag = 0 # set the output flag
        tm.model.update() # update the model

        tm.run_model()      # run the model
        tm.build_schedule() # builds the schedule
        tm.to_string()      # summarizes the solve for logging.
    """

    #: Slew tie-breaker weight in the objective. Same constant is used to
    #: recover total slew time inside :meth:`build_schedule`.

    # Constant that balances slew time vs. number of targets.
    # Interpretation: if dropping the highest priority target saves this many
    # minutes of slew time, drop it
    _SLEW_MINUTES_FOR_TOP_TARGET = 30
    _SLEW_IDLE_PENALTY_RATIO = 0.5  # idle-between weight, as fraction of slew penalty

    def __init__(
        self,
        requests,
        night_start,
        night_end,
        *,
        slew_fn,
        n_slots=1,
        slew_sample_cadence_min=30,
    ):
        _validate_requests(requests)
        self.requests = requests
        self.night_start = night_start
        self.night_end = night_end

        if len(self.requests) > 0:
            dfa = self.requests["first_available"] - night_start
            if dfa.min().to_value(u.s) > 0:
                logs.warning(
                    "min(first_available) is {:.1f} after night_start".format(
                        dfa.min().to(u.min)
                    )
                )

        self.slew_fn = slew_fn
        self.n_slots = n_slots
        self.slew_sample_cadence_min = slew_sample_cadence_min
        self.M = self.n_slots
        self.schedule = None
        self.stats = {}

    # ---------------------------------------------------------- helpers

    def _minutes_from_start(self, t):
        """``astropy.time.Time`` -> minutes since ``self.night_start``."""
        return (t.jd - self.night_start.jd) * 24 * 60

    def _time_from_minutes(self, mins):
        """Minutes since ``self.night_start`` -> ``Time`` (broadcast-friendly)."""
        return self.night_start + TimeDelta(np.asarray(mins) * 60, format="sec")

    # ------------------------------------------------------------- node setup

    def build_nodes(self):
        """Expand the request frame into the ``self.nodes`` DataFrame.

        Sets attributes:
            N, M, dur, nodes, multi_visit_groups.

        ``self.nodes`` carries only scalar dtypes. All times are in minutes
        from ``night_start``; sky coordinates appear as ``ra`` / ``dec``
        degrees. The :class:`~astropy.coordinates.SkyCoord` array used for
        slew evaluation is read from ``self.requests["coord"]`` and indexed
        by ``request_idx``.
        """
        self.dur = float(np.round(self._minutes_from_start(self.night_end), 0))

        # Local scratch DataFrame: scalar-only projection of self.requests
        # for the relational joins below. Built fresh here so the QTable
        # remains the single source of truth on the model.
        r = self.requests
        t_early = (r["first_available"] - self.night_start).to_value(u.min).astype(float)
        t_late = (r["last_available"] - self.night_start).to_value(u.min).astype(float)
        reqs = pd.DataFrame(
            {
                "unique_id": np.asarray(r["unique_id"], dtype=object),
                "ra": np.asarray(r["coord"].ra.deg, dtype=float),
                "dec": np.asarray(r["coord"].dec.deg, dtype=float),
                "t_early": t_early,
                "t_late": t_late,
                "t_visit": np.asarray(r["t_visit"].to_value(u.min), dtype=float),
                "n_intra_max": np.asarray(r["n_intra_max"], dtype=int),
                "tau_intra": np.asarray(r["tau_intra"].to_value(u.min), dtype=float),
                "priority": np.asarray(r["priority"], dtype=float),
                "request_idx": np.arange(len(r), dtype=np.int64),
                "is_anchor": False,
            }
        )

        # Attach visit_seq via a simple cross-join + filter.
        max_intra = int(reqs.n_intra_max.max()) if len(reqs) else 0
        visit_seq_table = pd.DataFrame(
            {
                "visit_seq": np.arange(max_intra, dtype=np.int64),
            }
        )
        visits = (
            reqs.merge(visit_seq_table, how="cross")
            .query("visit_seq < n_intra_max")
            .sort_values(["request_idx", "visit_seq"], kind="stable")
            .reset_index(drop=True)
        )

        anchor_template = {
            "unique_id": "",
            "request_idx": -1,
            "visit_seq": 0,
            "is_anchor": True,
            "t_early": 0.0,
            "t_late": self.dur,
            "t_visit": 0.0,
            "tau_intra": 0.0,
            "priority": 0.0,
            "n_intra_max": 0,
            "ra": np.nan,
            "dec": np.nan,
        }
        anchor_df = pd.DataFrame([anchor_template, anchor_template])

        self.nodes = pd.concat(
            [anchor_df.iloc[[0]], visits, anchor_df.iloc[[1]]],
            ignore_index=True,
        )
        self.N = len(self.nodes)

        # Multi-visit groups: keyed by unique_id, value is the list of node
        # ids (in visit_seq order) for requests with n_intra_max > 1.
        gb = self.nodes[~self.nodes["is_anchor"]].groupby("unique_id", sort=False)
        sizes = gb.size()
        self.multi_visit_groups = {
            uid: gb.get_group(uid).sort_values("visit_seq").index.tolist()
            for uid in sizes.index[sizes > 1]
        }

    def build_arcs(self, *, samples_per_slot=3):
        """Build precomputed arc catalog ``self.arcs``.

        For each ordered pair of real (non-anchor) nodes ``(i, j)`` and each
        slew slot ``m``, evaluate ``self.slew_fn`` at ``n_samples`` time points
        spread across the slot and store the worst-case slew. ``n_samples`` is
        ``max(dur / (M * slew_sample_cadence_min), samples_per_slot)``.

        Args:
            samples_per_slot (int): floor on the number of temporal samples per
                slew slot; the used count is ``int(max(...))`` of this and the
                value implied by ``dur``, ``self.M``, and
                ``slew_sample_cadence_min``.
        """
        if not hasattr(self, "nodes"):
            raise RuntimeError("call build_nodes() before build_arcs()")

        n_samples = int(
            max(
                self.dur / (self.M * self.slew_sample_cadence_min),
                samples_per_slot,
            )
        )
        total_samples = self.M * n_samples

        fractions = np.linspace(0, 1, total_samples)
        times = self.night_start + (self.night_end - self.night_start) * fractions

        # Index the per-request SkyCoord array by the real-node request_idx.
        non_anchor = self.nodes[~self.nodes.is_anchor]
        node_coords = self.requests["coord"][non_anchor["request_idx"].to_numpy()]
        n_real = len(node_coords)

        # Pair list: ordered (i, j) with i != j over real-node positions
        # (0-based). `id` in self.nodes for real nodes is `position + 1`.
        ii, jj = np.meshgrid(np.arange(n_real), np.arange(n_real), indexing="ij")
        mask = ii != jj
        i_pos = ii[mask]  # starting node of arc
        j_pos = jj[mask]  # ending node of arc
        n_pair = len(i_pos)

        # slew_fn returns shape (n_pair, T) = (n_pair, M*n_samples).
        tau = self.slew_fn(node_coords[i_pos], node_coords[j_pos], times)

        # Worst-case slew per slot: max over the n_samples axis -> shape (n_pair, M).
        tau_per_slot = tau.reshape(-1, self.M, n_samples).max(axis=2)

        p_idx, m_lev = np.mgrid[0:n_pair, 0 : self.M]
        i_id = i_pos[p_idx.ravel()] + 1  # account for anchor nodes
        j_id = j_pos[p_idx.ravel()] + 1  # account for anchor nodes

        uid = self.nodes["unique_id"].to_numpy()
        self.arcs = pd.DataFrame(
            {
                "i_id": uid[i_id],
                "j_id": uid[j_id],
                "t_slew": tau_per_slot.ravel(),
            },
            index=pd.MultiIndex.from_arrays(
                [i_id, j_id, m_lev.ravel()], names=("i", "j", "m")
            ),
        )

        # Slot bounds (minutes from night_start) used by build_model.
        slot_bound_fractions = np.linspace(0.0, 1.0, self.M + 1)
        self.w = self.dur * slot_bound_fractions

    # ------------------------------------------------------------- MILP build
    def build_model(self):
        """Construct the TTP MILP (Handley+ 2024, eqs. 2-9, B3, 10).

        Constraints are added in the order they appear in Handley+ 2024 §2.4.
        """
        if not hasattr(self, "arcs"):
            raise RuntimeError("call build_arcs() before build_model()")
        self.model = gp.Model("TTP")

        N, M = self.N, self.M
        nodes = self.nodes
        # O(1) slew lookup for hot loops (anchor arcs absent => .get(..., 0.0)).
        arcs_lookup = self.arcs["t_slew"].to_dict()

        self.Yi = self.model.addVars(
            range(N),
            vtype=GRB.BINARY,
            name="Yi",
        )
        self.Xijm = self.model.addVars(
            range(N),
            range(N),
            range(M),
            vtype=GRB.BINARY,
            name="Xijm",
        )
        self.tijm = self.model.addVars(
            range(N),
            range(N),
            range(M),
            vtype=GRB.CONTINUOUS,
            name="tijm",
        )
        self.ti = self.model.addVars(
            range(N),
            vtype=GRB.CONTINUOUS,
            lb=0,
            name="ti",
        )

        # Anchor visitation is always true
        self.model.addConstr(self.Yi[0] == 1, "start_anchor_visit")
        self.model.addConstr(self.Yi[N - 1] == 1, "end_anchor_visit")
        self.model.addConstr(self.ti[0] == 0.0, "anchor_start_time")

        # Not in Handley et al. We wish to force the first exposure to occur at the
        # start of the night or at the earliest feasible time, whichever is later.
        # This prevents idle time from being at the front of underfilled schedules.
        t_earliest = nodes[~nodes.is_anchor].t_early.min()
        t_start = max(0.0, t_earliest)
        for j in range(1, N - 1):
            t_visit_j = float(nodes.at[j, "t_visit"])
            for m in range(M):
                self.model.addGenConstrIndicator(
                    self.Xijm[0, j, m],
                    1,
                    self.ti[j],
                    GRB.EQUAL,
                    t_start + t_visit_j,
                    name=f"first_exposure_at_start_{j}_{m}",
                )

        # eq. 2 - exactly one arc out of the start anchor.
        self.model.addConstr(
            gp.quicksum(self.Xijm[0, j, m] for j in range(1, N) for m in range(M)) == 1,
            "start_anchor",
        )

        # eq. 3 - exactly one arc into the end anchor.
        self.model.addConstr(
            gp.quicksum(self.Xijm[i, N - 1, m] for i in range(N - 1) for m in range(M))
            == 1,
            "end_anchor",
        )

        # eq. 4 - visit indicator (one constraint per non-start node).
        for j in range(1, N):
            self.model.addConstr(
                gp.quicksum(self.Xijm[i, j, m] for i in range(N - 1) for m in range(M))
                == self.Yi[j],
                f"visit_once_{j}",
            )

        # eq. 5 - flow conservation (one constraint per real node).
        for k in range(1, N - 1):
            self.model.addConstr(
                gp.quicksum(self.Xijm[i, k, m] for i in range(N - 1) for m in range(M))
                - gp.quicksum(self.Xijm[k, j, m] for j in range(1, N) for m in range(M))
                == 0,
                f"flow_constr_{k}",
            )

        # eq. 6 - link ti to tijm (one constraint per non-end node).
        for i in range(N - 1):
            self.model.addConstr(
                self.ti[i]
                == gp.quicksum(
                    self.tijm[i, j, m] for j in range(1, N) for m in range(M)
                ),
                f"tijm_def_{i}",
            )

        # eq. 7 - exposure/slew time linking (one constraint per non-start
        # node). Anchor arcs (i=0) are MILP bookkeeping with t_slew = 0; the
        # .get(..., 0.0) fallback keeps arcs restricted to real arcs.
        for j in range(1, N):
            t_visit_j = nodes.at[j, "t_visit"]
            self.model.addConstr(
                self.ti[j]
                >= gp.quicksum(
                    self.tijm[i, j, m]
                    + (arcs_lookup.get((i, j, m), 0.0) + t_visit_j) * self.Xijm[i, j, m]
                    for i in range(N - 1)
                    for m in range(M)
                ),
                f"exp_constr_{j}",
            )

        # eq. 8 - slot bounds on tijm (one pair per (i, j, m)).
        for i, j, m in product(range(N), range(N), range(M)):
            self.model.addConstr(
                self.tijm[i, j, m] >= self.w[m] * self.Xijm[i, j, m],
                f"t_min_{i}_{j}_{m}",
            )
            self.model.addConstr(
                self.tijm[i, j, m] <= self.w[m + 1] * self.Xijm[i, j, m],
                f"t_max_{i}_{j}_{m}",
            )

        # eq. 9 - node accessibility (one pair per node).
        for i in range(N):
            row = nodes.loc[i]
            self.model.addConstr(
                self.ti[i] >= (row.t_early + row.t_visit) * self.Yi[i],
                f"rise_constr_{i}",
            )
            self.model.addConstr(
                self.ti[i] <= row.t_late * self.Yi[i],
                f"set_constr_{i}",
            )

        # eq. B3 - intra-night separation (multi-visit only).
        for indices in self.multi_visit_groups.values():
            for k in range(1, len(indices)):
                self.model.addConstr(
                    gp.quicksum(
                        self.tijm[indices[k], j, m]
                        for j in range(1, N)
                        for m in range(M)
                    )
                    >= gp.quicksum(
                        self.tijm[indices[k - 1], j, m]
                        for j in range(1, N)
                        for m in range(M)
                    )
                    + self.Yi[indices[k]] * nodes.at[indices[k], "tau_intra"],
                    f"intra_sep_constr_{indices[k - 1]}_{indices[k]}",
                )

        self.t_slew = self.model.addVar(lb=0.0, name="t_slew")
        self.model.addConstr(
            self.t_slew
            == gp.quicksum(
                arcs_lookup.get((i, j, m), 0.0) * self.Xijm[i, j, m]
                for i in range(1, N - 1)
                for j in range(1, N - 1)
                for m in range(M)
            ),
            "t_slew_def",
        )

        self.t_visit = self.model.addVar(lb=0.0, name="t_visit")
        self.model.addConstr(
            self.t_visit
            == gp.quicksum(
                nodes.at[j, "t_visit"] * self.Yi[j] for j in range(1, N - 1)
            ),
            "t_visit_def",
        )

        self.t_idle_between = self.model.addVar(lb=0.0, name="t_idle_between")
        self.model.addConstr(
            self.t_idle_between == self.ti[N - 1] - self.t_visit - self.t_slew,
            name="t_idle_between_def",
        )

        # eq. 10 - objective: priority-weighted visit count minus slew tie-breaker.
        P_max = float(self.nodes.loc[1 : N - 1, "priority"].max())
        slew_penalty = P_max / self._SLEW_MINUTES_FOR_TOP_TARGET
        self.model.setObjective(
            gp.quicksum(nodes.at[j, "priority"] * self.Yi[j] for j in range(1, N - 1))
            - slew_penalty * self.t_slew
            - slew_penalty * self._SLEW_IDLE_PENALTY_RATIO * self.t_idle_between,
            GRB.MAXIMIZE,
        )

        self.model.update()

    def run_model(self):
        """Solve the MILP and build ``schedule`` / ``stats``.

        On success, ``schedule`` is a DataFrame and ``stats`` is populated.
        If Gurobi has no incumbent (``SolCount == 0``), logs a warning and
        leaves ``schedule`` as ``None``.
        """
        if not hasattr(self, "model"):
            raise RuntimeError("call build_model() before run_model()")
        logs.info(f"Solving TTP for {self.N - 2} visits")
        t0 = time.time()
        self.model.optimize()
        if self.model.Status == GRB.INFEASIBLE:
            logs.critical("TTP infeasible; computing IIS.")
            self.model.computeIIS()
            for c in self.model.getConstrs():
                if c.IISConstr:
                    logs.critical(c.ConstrName)

        logs.info(f"TTP solve finished in {time.time() - t0:.3f}s")
        if self.model.SolCount == 0:
            logs.warning(
                "No incumbent TTP solution within time limit. "
                "Try raising ``self.model.params.TimeLimit`` before solving."
            )

    # ---------------------------------------------------------- post-process

    def build_schedule(self):
        """Walk the chosen path and populate ``schedule``.

        ``schedule`` is parallel to ``nodes`` (same index = node id ``i``),
        with ``scheduled``, ``t_slew``, and solve-time columns added.
        """
        if not hasattr(self, "model"):
            raise RuntimeError("call build_model() before build_schedule()")

        arcs_selected = []
        for (i, j, m), var in self.Xijm.items():
            if var.X > 0.5 and j != 0 and i != self.N - 1:
                arcs_selected.append({"i": i, "j": j, "m": m, "ti": self.ti[i].X})
        arcs_selected = pd.DataFrame(arcs_selected)

        # merge selected arcs with nodes, unvisited nodes will have NaN for ti
        schedule = pd.merge(
            self.nodes.query("~is_anchor"),
            arcs_selected,
            left_index=True,
            right_on=["i"],
            how="left",
        )
        schedule["t_start"] = schedule["ti"] - schedule["t_visit"]
        schedule["t_end"] = schedule["ti"]
        schedule["scheduled"] = ~schedule["ti"].isna()

        schedule = pd.merge(
            schedule,
            self.arcs["t_slew"],
            left_on=["i", "j", "m"],
            right_index=True,
            how="left",
        ).sort_values(by="t_start", na_position="last")
        schedule["order"] = range(len(schedule))

        # self.nodes has only scalar columns (ra/dec as floats, no SkyCoord);
        # schedule inherits that and is round-trip-safe through to_csv/to_hdf.
        self.schedule = schedule
        scheduled = self.schedule[self.schedule["scheduled"]]
        stats = {
            "dur": self.dur,
            "n_requested": self.N - 2,
            "n_scheduled": len(scheduled),
            "t_first_start": scheduled["t_start"].min(),
            "t_last_end": scheduled["t_end"].max(),
            "t_visit_sum": scheduled["t_visit"].sum(),
            "t_slew_sum": scheduled["t_slew"].sum(),
            "t_idle_sum": self.dur
            - scheduled["t_visit"].sum()
            - scheduled["t_slew"].sum(),
        }
        stats["t_idle_after_last"] = self.dur - stats["t_last_end"]
        stats["t_idle_before_last"] = stats["t_idle_sum"] - stats["t_idle_after_last"]
        self.stats = stats

    def to_string(self, *, header="Stats for TTP Solution"):
        """Return a human-readable summary of the solve from ``self.stats``."""
        if not self.stats:
            raise RuntimeError("call build_schedule() before to_string()")
        s = self.stats
        rows = [
            ("Observations Requested:", s["n_requested"], "d"),
            ("Observations Scheduled:", s["n_scheduled"], "d"),
            ("Observing Duration (min):", s["dur"], ".1f"),
            ("First Exposure Start (min):", s["t_first_start"], ".1f"),
            ("Last Exposure End (min):", s["t_last_end"], ".1f"),
            ("Visit Time (min):", s["t_visit_sum"], ".1f"),
            ("Slew Time (min):", s["t_slew_sum"], ".1f"),
            ("Idle Time (min):", s["t_idle_sum"], ".1f"),
            ("Idle After Last (min):", s["t_idle_after_last"], ".1f"),
            ("Idle Before Last (min):", s["t_idle_before_last"], ".1f"),
        ]
        label_w = max(len(label) for label, _, _ in rows)
        value_w = 7
        divider = "-" * (2 + label_w + 1 + value_w)
        fmt = f"  {{:<{label_w}}} {{:>{value_w}{{spec}}}}"

        lines = [header, divider]
        lines.extend(fmt.format(label, value, spec=spec) for label, value, spec in rows)
        lines.append(divider)
        return "\n".join(lines) + "\n"
