"""A/B test: does Yi-continuous + Symmetry=2 close the two-state gap faster?

Baseline = adversarial-2 v2 two-state at 600s with the *old* model (Yi binary,
default symmetry). This script solves the SAME v2 target set with:
  - the model change already in astroq/ttp/model.py (Yi continuous / implied-int)
  - Gurobi ``Symmetry=2`` (set here, not in production defaults)
at the same 600s budget, and prints presolved integer count, final gap, bound,
node count, and physical slew so the two can be compared directly.

Usage (astroq-testing env, repo root):
    python examples/ttp_twostate/adversarial-2_2026-06-12_symmetry-test/symmetry_experiment.py
"""

import importlib.util
import os

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
from astroq.scripts.demo_twostate import physical_slew_minutes

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.join(os.path.dirname(HERE), "adversarial-2_2026-06-12_runtime-600")
REQUEST_CSV = os.path.join(V2, "request_selected.csv")
NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")
RUNTIME_S = int(os.environ.get("TTP_RUNTIME_S", "600"))
SYMMETRY = int(os.environ.get("TTP_SYMMETRY", "2"))
LOGFILE = os.path.join(HERE, "gurobi_two_state_symmetry.log")

# Reuse the v2 availability builder so the model is identical to the baseline.
spec = importlib.util.spec_from_file_location(
    "advc", os.path.join(V2, "adversarial_compare.py")
)
advc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(advc)


def main():
    queue = HIRESCPS()
    df = pd.read_csv(REQUEST_CSV)
    requests = advc.build_requests_access(df, queue, NIGHT_START, NIGHT_END)

    tm = TTPModel(
        requests=requests,
        night_start=NIGHT_START,
        night_end=NIGHT_END,
        slew_fn=queue.slew_fn_state,
        n_slots=queue.nSlots,
        n_states=2,
    )
    tm.build_nodes()
    tm.build_arcs()
    tm.build_model()
    if os.path.exists(LOGFILE):
        os.remove(LOGFILE)
    tm.model.params.LogFile = LOGFILE
    tm.model.params.LogToConsole = 0
    tm.model.params.TimeLimit = RUNTIME_S
    tm.model.params.MIPGap = 0.005
    tm.model.params.Symmetry = SYMMETRY
    tm.model.update()
    n_int = sum(1 for v in tm.model.getVars() if v.VType != "C")
    n_cont = sum(1 for v in tm.model.getVars() if v.VType == "C")
    tm.run_model()
    tm.build_schedule()

    sc = tm.schedule[tm.schedule["scheduled"]]
    phys = physical_slew_minutes(tm.schedule, queue, NIGHT_START)
    lines = [
        "Two-state symmetry experiment (Yi continuous + Symmetry=%d)" % SYMMETRY,
        f"  time limit       : {RUNTIME_S}s",
        f"  vars (pre-solve) : {n_int} integer / {n_cont} continuous "
        f"(old model had Yi binary; presolved integer count is in the log)",
        f"  scheduled        : {len(sc)}/{tm.stats['n_requested']}",
        f"  modeled slew min : {float(tm.stats['t_slew_sum']):.2f}",
        f"  physical slew min: {float(phys):.2f}",
        f"  gurobi status    : {int(tm.model.Status)}",
        f"  mip gap          : {float(tm.model.MIPGap) * 100:.2f}%   "
        f"(baseline 600s: 7.58%)",
        f"  objective        : {float(tm.model.ObjVal):.3f}",
        f"  objective bound  : {float(tm.model.ObjBound):.3f}   (baseline: 513.98)",
        f"  nodes explored   : {int(tm.model.NodeCount)}   (baseline: 32204)",
        f"  solve runtime s  : {float(tm.model.Runtime):.1f}",
        f"  log              : {os.path.basename(LOGFILE)}",
    ]
    txt = "\n".join(lines)
    with open(os.path.join(HERE, "result.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
