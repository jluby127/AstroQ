# 01. Fill shortfall has a floor, and U258's seasonal mismatch was the only way past it

Recorded 2026-07-31. Artifacts in [`../data/01/`](../data/01/), provenance in
[`../data/01/PROVENANCE.md`](../data/01/PROVENANCE.md).

## Summary

The 2026B HIRES queue was subscribed at 120% in requested hours and still
produced a worst-case fill shortfall of 0.5%. Getting `fsf` above zero at all
required removing nights from the allocation. The reason is structural rather
than accidental: awards are derived from the allocation, so they sum to
capacity by construction, and fill is capped at the award while shortfall is
being optimized. Excess demand therefore never becomes contention. The one
program with a genuine shortfall, U258, was short not because it lost a
competition but because a third of its award fell in a part of the year when
its targets were not observable.

This matters for the paper because it undercuts the natural framing of program
balancing as a fairness problem between competing PIs. In this queue, balance
had almost nothing to arbitrate.

## Why programs cannot compete

Prep computes each program's award by summing the Keck blocks assigned to that
program (`format_keck_allocation_info` to `hours_by_program` to
`write_programs_csv`, in `astroq/driver.py`). The scheduler then pools all
blocks: any target may be scheduled in any allocated slot. Two consequences
follow.

First, the sum of awards equals the pooled capacity identically. There is no
slack and no overcommitment at the level of awarded time, regardless of how
much PIs requested.

Second, `F[p]` is capped at 1.0 during the shortfall and balance stages. The
125% ceiling only appears later, at prioritize and fill-empty. So while
shortfall is being optimized, no program may take more than its award.

Together these mean a solution where every program sits at exactly 100% is
always feasible in principle, and no program can improve by taking time from
another. `fsf > 0` therefore requires a program that *cannot physically use its
own allocation*: a structural mismatch, not a contested one.

The demand is real, and that is what makes the result interesting rather than
trivial. Every science program asked for more than it was given:

| program | awarded (h) | requested (h) | requested / awarded |
|---------|------------:|--------------:|--------------------:|
| C275    | 33.2 | 42.1 | 127% |
| D466    |  9.5 | 11.1 | 117% |
| H402    | 10.7 | 14.4 | 135% |
| N040    |  5.4 |  7.0 | 130% |
| N062    | 14.1 | 18.9 | 134% |
| N144    | 25.5 | 30.7 | 120% |
| U252    | 19.7 | 22.6 | 115% |
| U258    | 36.1 | 36.5 | 101% |
| Y063    | 10.3 | 14.8 | 144% |
| **total** | **164.5** | **198.1** | **120%** |

Yet the worst fill shortfall in the baseline run was 0.005, and the balance
stage could not improve on it.

## The U258 diagnosis

U258 held six nights totalling 36.11 h, but scheduled only 25.70 h, giving a
maximum feasible fill of 0.71. Per-night accounting
([`u258-night-usage-baseline.txt`](../data/01/u258-night-usage-baseline.txt),
regenerable with [`night_usage.py`](../data/01/night_usage.py)) shows every one
of its observations falling between Dec 7 and Jan 24, and four of its six
nights carrying zero U258 time:

| UT night | hours | U258 used | consumed instead by |
|----------|------:|----------:|---------------------|
| Aug 7  | 7.03 | 0.0 | C275 2.3, N062 1.1, N144 1.0, others |
| Aug 20 | 7.22 | 0.0 | C275 2.0, N144 1.3, U252 1.0, others |
| Sep 19 | 5.10 | 0.0 | N062 1.3, C275 1.0, N144 0.8, others |
| Oct 24 | 5.42 | 0.0 | C275 2.0, U252 1.7, N144 0.9, others |
| Dec 14 | 5.67 | 3.4 | shared |
| Dec 29 | 5.67 | 2.5 | shared |

So 24.77 h, roughly two thirds of the award, sat in the wrong part of the year
and was donated wholesale to the rest of the queue. The other programs were
not competing with U258; they were living off it. This is also why `fsf` was
zero everywhere: the surplus was large enough that everyone else reached their
ceiling, and U258's own shortfall was invisible to the balance objective
because `maff` had already been lowered to 0.71 to match what it could actually
do.

## The intervention

Dropping the two nights the PI identified as misplaced, Keck dates 2026-08-06
and 2026-09-18 (12.13 h), removed the capacity and the corresponding award
together, since both derive from the same Keck blocks.

| quantity | baseline | after dropping 2 nights |
|----------|---------:|------------------------:|
| allocation blocks | 29 | 27 |
| pooled capacity | 164.55 h | 152.42 h |
| U258 award | 36.10 h | 23.97 h |
| U258 `maff` | 0.71 | 1.00 |
| subscription | 120% | 130% |
| idle time in final plan | 8.78 h | 0.62 h |
| worst `fsf` after shortfall | 0.005 (C275) | 0.019 (U252) |
| worst `fsf` after balance | 0.005 (C275) | 0.005 (H402) |

Two things stand out. U258's `maff` reached exactly 1.00, meaning the program
could now fill its whole (smaller) award, and it stopped being a donor. And
idle time collapsed from 8.78 h to 0.62 h: the surplus had not been absorbed by
the other programs at all, it had simply gone unused, because they too were
capped at their awards.

## The sharper claim: 0.5% looks like a floor

The tempting reading is that balance rescued a skewed queue, cutting the worst
shortfall from 1.9% to 0.5%. The baseline argues against that. There the worst
`fsf` was *already* 0.005 coming out of the shortfall stage, and balance left
it at 0.005, changing nothing. Both allocations converge on the same 0.005.

That value therefore looks like a floor set by packing friction, slot
granularity, and cadence spacing, rather than a balance-quality result. The
honest statement is that balance moved the contended instance down to the same
floor the uncontended instance was already sitting on.

## What else changed between the two runs

This was an operational comparison, not a controlled experiment. Three things
differ at once, and the second one matters:

1. **Allocation.** 29 blocks / 164.55 h versus 27 blocks / 152.42 h. This is
   the intended variable.
2. **Shortfall solver budget.** Baseline had `TimeLimit 900`,
   `NoRelHeurTime 300`, and proved optimality in 326 s at a 0.0% gap. The test
   had `TimeLimit 300`, `NoRelHeurTime 120`, and stopped at the time limit with
   a 2.3% gap.
3. **Balance solver budget.** 60 s ending at a 62% gap versus 600 s ending at
   1.2%.

Because the test's shortfall stage never proved optimality, some part of the
1.9% to 0.5% improvement attributed to balance may simply be the balance stage
finishing optimization that shortfall left undone. The full config delta is in
[`config-gurobi-delta.txt`](../data/01/config-gurobi-delta.txt).

**This is the main threat to the finding and the first thing to fix.**

There is now direct evidence for it. The `fill-current-day` stage logs the
weighted shortfall objective before and after the tie-breaking stages:

| | stage-1 objective | cap (1.1x) | after all five stages |
|---|---:|---:|---:|
| baseline | 1001 | 1101.1 | 1001 |
| test | 1074 | 1181.4 | **1063** |

In the baseline the objective did not move: the later stages worked entirely
within the set of stage-1 optima. In the test the objective *fell* by 1.0%,
which is only possible because the stage-1 incumbent was not optimal, its 2.3%
gap leaving exactly this much room. So the later stages demonstrably performed
optimization that stage 1 should have done, which is the confound made
visible rather than inferred.

A useful side result: the `global_shortfall_slack` cap did not bind in either
run, so the 10% degradation the pipeline is licensed to accept was never spent.

## A prerequisite correction

None of this was interpretable until a prep bug was fixed. The U258 Nov 16 for
Dec 28 swap was injecting a full night (11.35 h) in place of a half night,
inflating total allocation by about 5.7 h and making U258 look better endowed
than it was. The half-for-half fix restored the total to 164.55 h and dropped
the U258 award from 41.78 h to 36.10 h. Worth a paragraph in the paper as an
example of how much ops data hygiene gates any quantitative claim about
completion.

## Open questions

1. **Resolve the confound.** Rerun the 27-block allocation with the shortfall
   stage at `TimeLimit 900` / `NoRelHeurTime 300` so it proves optimality, then
   see what balance actually contributes. Until this is done the balance result
   is not quotable.
2. **Is 0.5% a real floor?** Test whether it is set by slot granularity and
   cadence, or is just an artifact of the 0.005 `MIPGap`. The coincidence
   between the floor value and the MIPGap is suspicious and needs ruling out.
3. **What if programs could actually compete?** Raising `max_fill` above 1.0
   during shortfall and balance would let programs bid for more than their
   award and turn the 120% subscription into genuine contention. That is
   probably the more interesting scheduler to study, and it is a design
   question for the queue, not just a parameter.
4. **How common is seasonal mismatch?** U258 was the only 2026B program with a
   structural mismatch. Checking 2026A and earlier semesters would establish
   whether this is a recurring feature of how Keck assigns nights, which would
   make it a finding about the TAC process rather than about one PI.
