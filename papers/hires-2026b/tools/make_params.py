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

import csv
import datetime
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
DATA = os.path.join(PAPER, "data", "01")

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

        before, before_prog, after, after_prog = worst_fsf(solves)
        p[f"{prefix}-worst-fsf-before"] = before
        p[f"{prefix}-worst-fsf-before-prog"] = before_prog.replace("2026B_", "")
        p[f"{prefix}-worst-fsf-after"] = after
        p[f"{prefix}-worst-fsf-after-prog"] = after_prog.replace("2026B_", "")

    nights, hours = u258_donated("baseline-keck-blocks.csv")
    p["base-u258-donated-nights"] = nights
    p["base-u258-donated-hours"] = f"{hours:.2f}"

    p["drop2-dropped-hours"] = f"{float(p['base-capacity']) - float(p['drop2-capacity']):.2f}"
    p["drop2-dropped-nights"] = p["base-nblocks"] - p["drop2-nblocks"]

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
