# Finding 02: weather simulation pilot, paired arms

## What was run

Sequential night-by-night replan of the 2026B HIRES semester under simulated
weather, one weather realization (seed 1), two arms.

    papers/hires-2026b/tools/weather_sim.py --seed 1 --threads 6 \
        --source 2026B/2026-08-01/band1_baseline \
        --arm {balance,nobalance} --out <rundir>

Code: `5f99587` (adds the four-stage no-balance mode and the harness).
Run: 2026-07-31 23:26 through 2026-08-01 02:51 local.

Solver budgets: shortfall `TimeLimit=900` / `NoRelHeurTime=300`, balance
`TimeLimit=300`, `Threads=6`. Prioritize and fill-empty ran at the ops default
of 60 s.

## Arms

| arm | `[semester] mode` |
|---|---|
| balance | `shortfall,balance,prioritize,fill-empty,fill-current-day` |
| nobalance | `shortfall,prioritize,fill-empty,fill-current-day` |

The two differ only by the balance stage. The weighted-shortfall cap applies in
both, so the control is held to the same shortfall tolerance.

Weather is drawn from the seed and the night list alone, never from the arm, so
both arms lost an identical set of nights. Verified in the analysis:
`weather identical: True`.

## Files

- `metrics-seed1-{balance,nobalance}.csv` — one row per (night, program):
  awarded/past/projected hours, fill, maximum feasible fill, and two shortfall
  readings (`fsf` against nightly-recomputed maff, `fsf_frozen` against the
  night-1 value), plus per-stage solver gaps and runtimes.
- `console-seed1-{balance,nobalance}.txt` — per-night progress lines.

## Caveats

- **One realization.** This draw lost 5 of 29 allocated nights (17%), below the
  25% target, so it was a mild weather year. Conclusions about the endpoint
  need more seeds.
- **Balance did not always converge.** Its stage hit the 300 s limit on at
  least one night with a 62% gap. This biases against balance, so the measured
  advantage is a lower bound.
- The shortfall stage converged on 29/29 nights in both arms (max gap 0.3%), so
  the confound recorded in finding 01 is absent here.
- Only allocated nights are stepped; nights with no allocation cannot change
  the state.
- What counts as "observed" is the semester plan's visit assignment for that
  night, not a TTP night plan, so execution losses within a clear night are not
  modeled.
