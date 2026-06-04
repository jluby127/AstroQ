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
        slew_fn (callable): ``slew_fn(coord_a, coord_b, window_start,
            window_end) -> np.ndarray`` of worst-case slew minutes per
            (pair, window). ``coord_a`` and ``coord_b`` are pair-aligned
            1-D ``SkyCoord`` arrays of length ``P``; ``window_start`` and
            ``window_end`` are 1-D ``Time`` arrays of length ``M`` giving
            the bounds of each slew slot. The return has shape ``(P, M)``
            and its ``[k, m]`` entry is the slew time from ``coord_a[k]``
            to ``coord_b[k]`` over ``[window_start[m], window_end[m]]``.
            The implementation owns the per-window sampling cadence and
            reduction (max, mean, percentile, ...).
        n_slots (int): number of TTP slew slots ``M`` (Handley+ 2024 §2.2).
            ``n_slots=1`` is the recommended default.

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

    #: Per-column spec for the ``requests`` QTable. Value is either an astropy
    #: type the column must be an instance of, or a time-equivalent unit the
    #: column's ``.unit`` must be convertible to. ``None`` skips type checking.
    _COLUMN_SPECS = {
        "unique_id":       None,      # str, primary key
        "coord":           SkyCoord,  # ICRS
        "first_available": Time,
        "last_available":  Time,
        "t_visit":         u.s,       # Quantity with time units
        "n_intra_max":     None,      # int
        "tau_intra":       u.s,       # Quantity with time units
        "priority":        None,      # float
    }

    def __init__(
        self,
        requests,
        night_start,
        night_end,
        *,
        slew_fn,
        n_slots=1,
    ):
        if not isinstance(requests, QTable):
            raise TypeError(
                "`requests` must be an astropy.table.QTable; got "
                f"{type(requests).__name__}"
            )
        missing = [c for c in self._COLUMN_SPECS if c not in requests.colnames]
        if missing:
            raise ValueError(f"`requests` missing required columns: {missing}")
        for col, spec in self._COLUMN_SPECS.items():
            if spec is None or len(requests) == 0:
                continue
            column = requests[col]
            if isinstance(spec, u.UnitBase):
                unit = column.unit
                if unit is None or not unit.is_equivalent(spec):
                    raise TypeError(
                        f"`{col}` column must be a Quantity with units "
                        f"equivalent to {spec}; got unit={unit}"
                    )
            elif not isinstance(column, spec):
                raise TypeError(
                    f"`{col}` column must be an astropy {spec.__name__}; "
                    f"got {type(column).__name__}"
                )

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
        self.dur_min = float(
            np.round((night_end - night_start).to_value(u.min), 0)
        )
        self.M = self.n_slots
        self.schedule = None
        self.stats = {}

    def build_nodes(self):
        """Expand the request frame into the ``self.nodes`` DataFrame.

        Sets attributes:
            N, M, nodes, multi_visit_groups.

        ``self.nodes`` carries only scalar dtypes. All times are in minutes
        from ``night_start``; sky coordinates appear as ``ra`` / ``dec``
        degrees. The :class:`~astropy.coordinates.SkyCoord` array used for
        slew evaluation is read from ``self.requests["coord"]`` and indexed
        by ``request_idx``.
        """
        # Local scratch DataFrame: scalar-only projection of self.requests
        # for the relational joins below. Built fresh here so the QTable
        # remains the single source of truth on the model.
        r = self.requests
        rdf = self.requests["unique_id","n_intra_max","priority"].to_pandas()
        rdf["t_early"] = (r["first_available"] - self.night_start).to_value(u.min).astype(float)
        rdf["t_late"] =  (r["last_available"] - self.night_start).to_value(u.min).astype(float)
        rdf["tau_intra"] = r["tau_intra"].to_value(u.min).astype(float)
        rdf["t_visit"] = r["t_visit"].to_value(u.min).astype(float)
        rdf["request_idx"] = np.arange(len(r), dtype=np.int64)
        rdf["is_anchor"] = False

        # Attach visit_seq via a simple cross-join + filter.
        max_intra = int(rdf.n_intra_max.max()) if len(rdf) else 0
        visit_seq_table = pd.DataFrame(
            {
                "visit_seq": np.arange(max_intra, dtype=np.int64),
            }
        )
        visits = (
            rdf.merge(visit_seq_table, how="cross")
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
            "t_late": self.dur_min,
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

    def build_arcs(self):
        """Build precomputed arc catalog ``self.arcs``.

        Splits ``[night_start, night_end]`` into ``M`` equal-length windows
        and asks ``self.slew_fn`` for the worst-case slew (minutes) on each
        ordered pair of real nodes within each window. The per-window
        sampling cadence and reduction policy live inside ``slew_fn``; see
        the class docstring for the callable contract.
        """
        if not hasattr(self, "nodes"):
            raise RuntimeError("call build_nodes() before build_arcs()")

        # Equal-length window bounds across the night.
        fractions = np.linspace(0.0, 1.0, self.M + 1)
        slot_bounds = self.night_start + (self.night_end - self.night_start) * fractions
        window_start = slot_bounds[:-1]
        window_end = slot_bounds[1:]

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

        # slew_fn returns shape (n_pair, M): worst-case slew per (pair, window).
        tau_per_slot = self.slew_fn(
            node_coords[i_pos], node_coords[j_pos], window_start, window_end
        )

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
        self.w = self.dur_min * fractions

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
            "dur_min": self.dur_min,
            "n_requested": self.N - 2,
            "n_scheduled": len(scheduled),
            "t_first_start": scheduled["t_start"].min(),
            "t_last_end": scheduled["t_end"].max(),
            "t_visit_sum": scheduled["t_visit"].sum(),
            "t_slew_sum": scheduled["t_slew"].sum(),
            "t_idle_sum": self.dur_min
            - scheduled["t_visit"].sum()
            - scheduled["t_slew"].sum(),
        }
        stats["t_idle_after_last"] = self.dur_min - stats["t_last_end"]
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
            ("Observing Duration (min):", s["dur_min"], ".1f"),
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
