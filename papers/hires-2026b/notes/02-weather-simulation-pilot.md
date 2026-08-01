# 02 — Balance bounds the worst excursion, not the endpoint

Evidence: `data/03/`. Code `5f99587`. One weather realization (seed 1), paired
arms differing only by the balance stage.

## Claim

Under simulated weather, the balance stage roughly halves the worst fill
shortfall averaged over the semester and caps its worst excursion at about a
third of the unbalanced value, but in this realization it does not change where
programs finish.

## Numbers

Worst fill shortfall over all programs, per night, balance vs no balance:

| | balance | no balance |
|---|---|---|
| mean over 29 nights | 0.017 | 0.030 |
| worst excursion | 0.050 | 0.140 |
| final night | 0.030 | 0.030 |

The trajectory has three regimes:

- **Nights 0–16 (Aug–early Oct).** Identical, both at 0.000. Weather has not
  yet removed enough time for programs to compete.
- **Nights 17–23 (mid Oct–mid Dec).** Balance is clearly better, peaking at
  night 19 with 0.050 against 0.140, a factor of 2.8.
- **Nights 24–28 (late Dec–Jan).** The arms converge and finish identical.

Final per-program fill differs by at most 0.03. Balance ends with U258 at 0.62
and the strong programs at 0.97–1.00; no balance ends with U258 also at 0.62
and the strong programs at 1.00. Total delivered time is the same to 0.1 hr
(144.5 vs 144.6 of 164.5 awarded).

## Reading

The mid-semester advantage is real and it works the way the design intends:
balance holds up the weakest program by shaving points off programs that would
otherwise run to exactly 1.00. At the peak of the effect, U258 sits at 0.66
under balance against 0.58 without it.

That advantage does not survive to the end of the semester here. The reason is
the same seasonal mismatch recorded in finding 01: U258's ceiling is set by
when its targets rise, not by contention with other programs. Balance can
change the path but not the destination for a program whose limit is
astronomical rather than competitive.

Whether that generalizes is the open question. This draw lost only 17% of
nights against a 25% target, so it was a mild year and the queue had slack to
recover in January. A harsher realization, or one that concentrates losses late
in the semester when there is no room left to recover, is where the endpoint
would be expected to separate.

## Confounds checked

- Shortfall converged on 29/29 nights in both arms, max gap 0.3%. The
  non-convergence confound from finding 01 is absent.
- Weather sequences verified identical across arms.
- Balance hit its own 300 s time limit on at least one night at a 62% gap,
  which biases against balance. The measured advantage is a lower bound.

## Next

More seeds, and seeds chosen or filtered for heavier and later-weighted losses.
A single realization cannot distinguish "balance does not change the endpoint"
from "this particular year was too easy to tell".
