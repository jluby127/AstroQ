"""Generate params.tex from the recorded run artifacts under data/.

Every measured number quoted in the paper goes through \\param{key} so that no
figure is ever hand-typed into the prose. This script derives those keys from
the artifacts in data/01/ and writes params.tex in the expl3 \\str_case:nnF
format inherited from the KPF-CC paper, including its red XX fallback so an
unresolved key is visible in the PDF instead of silently disappearing.

Run from the paper directory:
    python tools/make_params.py

Key namespaces follow the prior paper's ns12-/ns6-/ns3- convention:
    base-   the 29-block 2026B baseline allocation
    drop2-  the 27-block allocation with two donated U258 nights removed
"""

import configparser
import csv
import datetime
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
DATA = os.path.join(PAPER, "data", "01")
DATA02 = os.path.join(PAPER, "data", "02")
DATA03 = os.path.join(PAPER, "data", "03")

# Programs that are bookkeeping rather than science: the filler backstop and a
# program with no requests in the CPS queue. Excluded from subscription totals.
NON_SCIENCE = {"2026B_E475", "2026B_C362"}

# U258 observes only from December onward; blocks before this date are the ones
# it cannot use and therefore donates to the rest of the queue.
U258_SEASON_START = datetime.date(2026, 12, 4)


def read(name):
    with open(os.path.join(DATA, name)) as fh:
        return fh.read()


def allocation_stats(name):
    blocks = list(csv.DictReader(open(os.path.join(DATA, name))))
    hours = sum(
        (datetime.datetime.fromisoformat(b["stop"])
         - datetime.datetime.fromisoformat(b["start"])).total_seconds() / 3600.0
        for b in blocks
    )
    return len(blocks), hours


def programs(name):
    return {r["program"]: r for r in csv.DictReader(open(os.path.join(DATA, name)))}


def report_table(name):
    """Parse the per-program statistics block out of a run report."""
    rows = {}
    for line in read(name).splitlines():
        parts = line.split()
        if not parts or not parts[0].startswith("2026B_"):
            continue
        rows[parts[0]] = {
            "aw": float(parts[1]), "req": float(parts[2]),
            "past": float(parts[3]), "proj": float(parts[4]),
            "proj_pct": parts[8], "maff_pct": parts[9], "fsf_pct": parts[10],
        }
    return rows


def report_scalar(name, label):
    m = re.search(rf"^{re.escape(label)}\s+(\d+)$", read(name), re.M)
    return int(m.group(1)) if m else None


def final_solve(name, stage):
    """Last solve line for a stage: (status, objective, bound, gap, seconds)."""
    pattern = rf"{stage} solve: (\w[\w ]*), objective=([\d.]+), bound=([\d.]+), gap=([\d.]+)%, (\d+)s"
    matches = re.findall(pattern, read(name))
    return matches[-1] if matches else None


def worst_fsf(name):
    m = re.search(
        r"balance: worst fill shortfall ([\d.]+) \((\S+)\) -> ([\d.]+) \((\S+)\)",
        read(name),
    )
    return m.groups() if m else None


def usage_totals(name):
    m = re.search(
        r"total allocated ([\d.]+) h\s+scheduled ([\d.]+) h\s+idle ([\d.]+) h", read(name)
    )
    return tuple(float(x) for x in m.groups())


def u258_donated(name):
    """Hours in U258-owned blocks that fall before it can observe anything."""
    total = 0.0
    count = 0
    for b in csv.DictReader(open(os.path.join(DATA, name))):
        if b["ProjCode"].strip(",") != "U258":
            continue
        start = datetime.datetime.fromisoformat(b["start"])
        if start.date() < U258_SEASON_START:
            total += (datetime.datetime.fromisoformat(b["stop"]) - start).total_seconds() / 3600.0
            count += 1
    return count, total


def collect():
    p = {}

    for tag, prefix in (("baseline", "base"), ("test", "drop2")):
        nblocks, capacity = allocation_stats(f"{tag}-allocation.csv")
        p[f"{prefix}-nblocks"] = nblocks
        p[f"{prefix}-capacity"] = f"{capacity:.2f}"

        progs = programs(f"{tag}-programs.csv")
        p[f"{prefix}-u258-award"] = f"{float(progs['2026B_U258']['hours']):.2f}"
        p[f"{prefix}-u258-maff"] = f"{float(progs['2026B_U258']['max_feasible_fill']):.2f}"

        shortfall = report_table(f"{tag}-report-shortfall.txt")
        science = {k: v for k, v in shortfall.items() if k not in NON_SCIENCE}
        req = sum(v["req"] for v in science.values())
        aw = sum(v["aw"] for v in science.values())
        p[f"{prefix}-requested"] = f"{req:.1f}"
        p[f"{prefix}-awarded"] = f"{aw:.1f}"
        p[f"{prefix}-subscription"] = f"{100 * req / aw:.0f}"
        p[f"{prefix}-nprograms"] = len(science)
        p[f"{prefix}-u258-proj"] = f"{shortfall['2026B_U258']['proj']:.2f}"

        p[f"{prefix}-fill-shortfall-stage"] = report_scalar(
            f"{tag}-report-shortfall.txt", "Future fill factor")
        p[f"{prefix}-fill-balance-stage"] = report_scalar(
            f"{tag}-report-balance.txt", "Future fill factor")
        p[f"{prefix}-alloc-slots"] = report_scalar(
            f"{tag}-report-shortfall.txt", "Total allocated slots")

        alloc, sched, idle = usage_totals(f"u258-night-usage-{tag}.txt")
        p[f"{prefix}-scheduled"] = f"{sched:.2f}"
        p[f"{prefix}-idle"] = f"{idle:.2f}"

        solves = f"solve-lines-{tag}.txt"
        status, obj, bound, gap, secs = final_solve(solves, "shortfall")
        p[f"{prefix}-shortfall-gap"] = gap
        p[f"{prefix}-shortfall-seconds"] = secs
        p[f"{prefix}-shortfall-objective"] = obj
        status, obj, bound, gap, secs = final_solve(solves, "balance")
        p[f"{prefix}-balance-gap"] = gap
        p[f"{prefix}-balance-seconds"] = secs

        # How much the four tie-breaking stages actually degraded the stage-1
        # objective, versus how much the slack would have permitted.
        m = re.search(
            r"shortfall objective=([\d.]+) cap=([\d.]+) "
            r"post-step shortfall objective=([\d.]+)",
            read(solves),
        )
        start, allowed, final = (float(x) for x in m.groups())
        p[f"{prefix}-shortfall-cap"] = f"{allowed:.1f}"
        p[f"{prefix}-shortfall-final"] = f"{final:.0f}"
        p[f"{prefix}-shortfall-drift"] = f"{100 * (final - start) / start:+.1f}"

        before, before_prog, after, after_prog = worst_fsf(solves)
        p[f"{prefix}-worst-fsf-before"] = before
        p[f"{prefix}-worst-fsf-before-prog"] = before_prog.replace("2026B_", "")
        p[f"{prefix}-worst-fsf-after"] = after
        p[f"{prefix}-worst-fsf-after-prog"] = after_prog.replace("2026B_", "")

    # Pipeline constants describing the method rather than a result. Read from
    # the recorded config so the paper cannot drift from what was actually run.
    cfg = configparser.ConfigParser()
    cfg.read(os.path.join(DATA, "test-config.ini"))
    p["pipeline-stages"] = len(cfg.get("semester", "mode").split(","))
    p["pipeline-slack"] = cfg.get("semester.balance", "global_shortfall_slack")
    p["pipeline-slack-pct"] = f"{100 * (float(p['pipeline-slack']) - 1):.0f}"
    p["pipeline-hold-alpha"] = cfg.get("semester.prioritize", "hold_fill_alpha")
    p["pipeline-mipgap"] = cfg.get("semester.default.gurobi", "MIPGap")
    p["pipeline-slot-minutes"] = cfg.get("semester", "slot_size")

    caps = {float(r["max_fill"]) for r in programs("test-programs.csv").values()}
    p["pipeline-maxfill-late"] = f"{max(caps):.2f}"

    nights, hours = u258_donated("baseline-keck-blocks.csv")
    p["base-u258-donated-nights"] = nights
    p["base-u258-donated-hours"] = f"{hours:.2f}"

    p["drop2-dropped-hours"] = f"{float(p['base-capacity']) - float(p['drop2-capacity']):.2f}"
    p["drop2-dropped-nights"] = p["base-nblocks"] - p["drop2-nblocks"]

    p.update(collect_weather())
    p.update(collect_paired())
    return p


def collect_paired():
    """Metrics from the paired balance/no-balance weather run (data/03/)."""
    p = {}
    # The filler backstop's award is nominal rather than allocated time, so it
    # is not a program whose shortfall means anything.
    exclude = {"2026B_E475"}
    arms = {}
    for arm in ("balance", "nobalance"):
        path = os.path.join(DATA03, f"metrics-seed1-{arm}.csv")
        if not os.path.exists(path):
            return p
        d = pd.read_csv(path)
        arms[arm] = d[(d["awarded_hr"] > 0) & (~d["program"].isin(exclude))]

    p["pair-seed"] = "1"
    p["pair-weather-p"] = "25"
    ref = arms["balance"]
    nights = ref.groupby("night_index")["clear"].first()
    p["pair-nights"] = len(nights)
    p["pair-lost"] = int((~nights).sum())
    p["pair-lost-frac"] = f"{100 * (~nights).mean():.0f}"

    last = ref["night_index"].max()
    for arm, d in arms.items():
        tag = "bal" if arm == "balance" else "nobal"
        worst = d.groupby("night_index")["fsf"].max()
        frozen = d.groupby("night_index")["fsf_frozen"].max()
        p[f"pair-fsf-mean-{tag}"] = f"{worst.mean():.3f}"
        p[f"pair-fsf-peak-{tag}"] = f"{worst.max():.3f}"
        p[f"pair-fsf-final-{tag}"] = f"{worst.iloc[-1]:.3f}"
        p[f"pair-frozen-peak-{tag}"] = f"{frozen.max():.3f}"
        p[f"pair-frozen-final-{tag}"] = f"{frozen.iloc[-1]:.3f}"

        final = d[d["night_index"] == last]
        p[f"pair-hours-{tag}"] = f"{final['proj_hr'].sum():.1f}"

        # Per-program floor of projected completion: the worst point of the
        # weather, which is where the arms differ most.
        by_program = d.groupby("program")["fill"]
        for code in ("C275", "U258"):
            p[f"pair-{code.lower()}-dip-{tag}"] = (
                f"{by_program.min()[f'2026B_{code}']:.2f}"
            )

        gap = d.groupby("night_index")["shortfall_gap"].first()
        optimal = d.groupby("night_index")["shortfall_optimal"].first()
        p[f"pair-shortfall-converged-{tag}"] = int(optimal.sum())
        p[f"pair-shortfall-gap-max-{tag}"] = f"{100 * gap.max():.1f}"

        if arm == "balance":
            bgap = d.groupby("night_index")["balance_gap"].first()
            p["pair-balance-gap-max"] = f"{100 * bgap.max():.0f}"
            p["pair-balance-loose"] = int((bgap > 0.01).sum())
            p["pair-balance-timelimit"] = "300"

    p["pair-awarded"] = f"{ref[ref['night_index'] == last]['awarded_hr'].sum():.1f}"

    # The night on which the arms first separate durably, and the split of the
    # queue into programs balance squeezes and programs it protects.
    worst = {arm: d.groupby("night_index")["fsf"].max() for arm, d in arms.items()}
    dates = ref.groupby("night_index")["date"].first()
    apart = (worst["nobalance"] - worst["balance"]).abs() > 0.005
    # An isolated night crosses the threshold well before the arms really
    # separate, so require the start of a run of consecutive nights.
    idx = list(apart.index)
    runs = [
        i
        for n, i in enumerate(idx[:-1])
        if apart[i] and apart[idx[n + 1]] and not (n and apart[idx[n - 1]])
    ]
    p["pair-split-date"] = dates[runs[0]]
    p["pair-blip-date"] = dates[apart[apart].index.min()]

    dip = {arm: d.groupby("program")["fill"].min() for arm, d in arms.items()}
    delta = dip["balance"] - dip["nobalance"]
    squeezed = delta < -0.005
    p["pair-squeezed"] = int(squeezed.sum())
    p["pair-protected"] = int((delta > 0.005).sum())
    p["pair-squeezed-lo"] = f"{dip['balance'][squeezed].min():.2f}"
    p["pair-squeezed-hi"] = f"{dip['balance'][squeezed].max():.2f}"

    # How concentrated the unbalanced shortfall is: programs that carry
    # essentially none of it without balance but pick some up with it.
    means = {arm: d.groupby("program")["fsf"].mean() for arm, d in arms.items()}
    spared = means["nobalance"] < 0.005
    p["pair-spared-nobal"] = int(spared.sum())
    p["pair-nprograms"] = len(means["nobalance"])
    p["pair-spared-mean-bal"] = f"{means['balance'][spared].mean():.3f}"
    return p


def collect_weather():
    """Metrics from the rolling weather simulation (data/02/)."""
    p = {}
    p["wx-weather-p"] = "0.25"

    toy_path = os.path.join(DATA02, "toy-combined-metrics.csv")
    if os.path.exists(toy_path):
        toy = pd.read_csv(toy_path)
        p["wx-toy-seeds"] = len(toy)
        p["wx-toy-nights"] = int(toy["n_nights"].iloc[0])
        p["wx-toy-replans-mean"] = f"{toy['n_replans'].mean():.0f}"
        p["wx-toy-fill-spread-mean"] = f"{toy['fill_realized_spread'].mean():.3f}"
        p["wx-toy-fill-spread-max"] = f"{toy['fill_realized_spread'].max():.3f}"
        p["wx-toy-fill-realized-mean"] = f"{toy['fill_realized_mean'].mean():.3f}"
        p["wx-toy-weather-frac-mean"] = f"{100 * toy['weather_fraction'].mean():.0f}"
        p["wx-toy-balance-gap-p90-max"] = f"{toy['balance_gap_p90'].max():.1f}"
        p["wx-toy-fsf-start-max"] = f"{toy['fsf_start_max'].max():.3f}"

    hires_path = os.path.join(DATA02, "hires-combined-metrics.csv")
    if os.path.exists(hires_path):
        hires = pd.read_csv(hires_path)
        row = hires.iloc[0]
        p["wx-hires-seed"] = int(row["seed"])
        p["wx-hires-nights"] = int(row["n_nights"])
        p["wx-hires-replans"] = int(row["n_replans"])
        p["wx-hires-weather-frac"] = f"{100 * row['weather_fraction']:.0f}"
        p["wx-hires-fill-spread"] = f"{row['fill_realized_spread']:.3f}"
        p["wx-hires-fill-realized-mean"] = f"{row['fill_realized_mean']:.3f}"
        p["wx-hires-fsf-start-max"] = f"{row['fsf_start_max']:.3f}"
        p["wx-hires-balance-gap-p90"] = f"{row['balance_gap_p90']:.1f}"
        p["wx-hires-shortfall-gap-p90"] = f"{row['shortfall_gap_p90']:.1f}"

    p["wx-shortfall-timelimit"] = "900"
    p["wx-shortfall-norel"] = "300"
    p["wx-balance-timelimit"] = "600"

    return p


def write(params):
    lines = [
        "% GENERATED by tools/make_params.py -- do not edit by hand.",
        "% Regenerate with: make params",
        "\\usepackage{xparse}",
        "\\usepackage{xcolor}",
        "\\usepackage{etoolbox}",
        "\\ExplSyntaxOn",
        "\\NewDocumentCommand{\\param}{m}{",
        "\\str_case:nnF {#1} {",
    ]
    for key in sorted(params):
        lines.append(f"{{{key}}}{{{params[key]}}}%")
    lines += [
        "}{{\\color{red}XX}}  % fallback",
        "}",
        "\\ExplSyntaxOff",
        "",
    ]
    out = os.path.join(PAPER, "params.tex")
    with open(out, "w") as fh:
        fh.write("\n".join(lines))
    return out, len(params)


def check(params):
    """Verify every \\param{} key used in the .tex sources is defined.

    Checking the sources rather than the rendered PDF means this needs no
    external tooling and reports which key is missing, instead of just noting
    that a red XX appeared somewhere.
    """
    used = {}
    for name in os.listdir(PAPER):
        if not name.endswith(".tex") or name == "params.tex":
            continue
        text = open(os.path.join(PAPER, name)).read()
        # Drop comments first, or documentation that mentions \param{key} in
        # prose would be reported as an unresolved key.
        text = re.sub(r"(?<!\\)%.*", "", text)
        for key in re.findall(r"\\param\{([^}]*)\}", text):
            used.setdefault(key, set()).add(name)

    missing = {k: v for k, v in used.items() if k not in params}
    for key in sorted(missing):
        print(f"MISSING {key}  (used in {', '.join(sorted(missing[key]))})")

    unused = sorted(set(params) - set(used))
    print(f"{len(used)} keys used, {len(params)} defined, {len(unused)} defined but unused")
    if missing:
        print(f"FAIL: {len(missing)} unresolved \\param key(s)")
        return 1
    print("OK: every \\param key resolves")
    return 0


if __name__ == "__main__":
    import sys

    params = collect()
    if "--check" in sys.argv:
        raise SystemExit(check(params))

    path, n = write(params)
    print(f"wrote {path} with {n} keys")
    for key in sorted(params):
        print(f"  {key:32} {params[key]}")
