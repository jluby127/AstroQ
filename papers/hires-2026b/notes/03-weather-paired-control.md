# 03 — Balance redistributes shortfall across the queue; it does not create time

Evidence: [`../data/03/`](../data/03/), provenance in
[`../data/03/PROVENANCE.md`](../data/03/PROVENANCE.md). Code `5f99587`.
Written up in the paper as \S "Balance Under Weather".

Companion to [note 02](02-weather-rolling-simulation.md), which runs the same
kind of simulation on a symmetric toy instance and replans only after weathered
nights. This note is the real 2026B instance, replanning every allocated night,
with a paired control arm.

## Claim

Under simulated weather, the balance stage takes fill shortfall off the two
programs that would otherwise absorb all of it and spreads it across the whole
queue. It does not change how much time the queue delivers, and in this
realization it does not change where the worst-off program finishes.

## Design

Step through the 29 allocated nights in order. On each night:

1. rebuild the semester plan from scratch using only what is known that day,
   including a fresh `compute-max-fill`
2. flip a Bernoulli coin, `p = 0.25`, on whether the night is wholly lost
3. if clear, append every visit the plan assigned to that night to `past.csv`
   and carry it into the next night's replan; if lost, bank nothing

Weather is therefore modeled as the *absence of data* rather than as a penalty
term, which is what the scheduler actually sees.

Two design choices carry the interpretation:

**Paired arms.** Run twice on the same weather realization: the five-stage
pipeline, and a four-stage pipeline identical except that balance is omitted.
The weighted-shortfall cap applies in both, so the control is held to the same
tolerance on `Z`, and balance is the only difference. The weather draw depends
only on the seed and the night list, never on the arm, so both arms lose an
identical set of nights (verified: `weather identical: True`).

This control did not exist before `5f99587` — the only available modes were
stage 1 alone and all five stages. Without it the experiment is unfalsifiable,
because every program's fill falls when nights disappear regardless of balance.
That is the lesson of [note 01](01-fsf-floor-and-u258-mismatch.md).

**maff recomputed nightly.** As operations does. This means maff falls as
nights are lost, so `fsf` measures the gap to what a program can *still*
achieve, not to what it was promised. Both readings are recorded: `fsf` against
nightly maff, and `fsf_frozen` against the night-1 value.

## Result

The realization lost 5 of 29 nights (17%). Shortfall proved optimality on 29/29
nights in both arms, max gap 0.3%, so note 01's confound is absent.

Worst fill shortfall over all programs:

| | balance | no balance |
|---|---|---|
| mean over 29 nights | 0.017 | 0.030 |
| peak | 0.050 | 0.140 |
| final night | 0.030 | 0.030 |
| final, vs frozen night-1 maff | 0.130 | 0.140 |

Three regimes: both at zero through mid-September apart from one isolated night
(2026-09-18, nobalance 0.010), durable separation from 2026-10-18 through
2026-12-11, then convergence from 2026-12-13 on.

### The part that matters: where the shortfall goes

The envelope above hides the real finding. Per program, mean `fsf` over the
semester:

| | C275 | U258 | U252 | H402 | N040 | N062 | N144 | Y063 | D466 |
|---|---|---|---|---|---|---|---|---|---|
| balance | 0.013 | 0.016 | 0.016 | 0.016 | 0.014 | 0.016 | 0.016 | 0.013 | 0.000 |
| no balance | 0.020 | 0.019 | 0.012 | 0.000 | 0.000 | 0.001 | 0.000 | 0.000 | 0.000 |

**Without balance, six of nine programs carry no shortfall at all** and the
entire cost of the lost nights lands on C275 and U258. With balance it is
spread across the queue at a fairly uniform 0.04–0.05 per night.

In completion terms, the floor each program reaches at the worst of the weather:

| | C275 | U258 | D466 | U252 | H402 | N144 | Y063 | N040 | N062 |
|---|---|---|---|---|---|---|---|---|---|
| balance | 0.87 | 0.62 | 1.00 | 0.93 | 0.95 | 0.95 | 0.96 | 0.96 | 0.95 |
| no balance | 0.81 | 0.58 | 0.98 | 0.92 | 1.00 | 0.99 | 1.00 | 0.99 | 0.98 |
| delta | +0.06 | +0.04 | +0.02 | +0.01 | −0.05 | −0.04 | −0.04 | −0.03 | −0.03 |

The queue splits cleanly by the sign: **4 programs protected, 5 squeezed**, the
squeezed ones all pulled to 0.95–0.96 from at or near their full award. That is
the signature of a max-min objective and nothing else in the pipeline produces
it.

Total delivered time is the same either way: 144.5 vs 144.6 hr of 164.5
awarded. Balance decides *who* gives up time; it cannot create time. Realized
completion trajectories are nearly identical between arms, since accumulation
is set by which nights the weather took.

## Why the endpoints converge

U258 finishes at 0.62 in both arms because its ceiling is astronomical rather
than competitive — the seasonal mismatch of note 01. Balance can change the
path but not the destination for a program whose limit is when its targets
rise.

## Confounds and limits

- **One realization**, and a mild one: 17% of nights lost against a 25%
  target, leaving room to recover in January. A single draw cannot distinguish
  "balance does not change the endpoint" from "this year was too easy to tell".
  The redistribution result does not depend on the draw; the convergence result
  might.
- **Balance usually did not converge.** Its gap exceeded 1% on 8 of 29 nights,
  reaching 62% on one, against a 300 s limit. This biases against balance, so
  the reported advantage is a lower bound on what the objective would deliver
  if solved properly. Worth rerunning with a larger budget before trusting the
  magnitudes.
- Shortfall converged everywhere; verified, not assumed.
- Only allocated nights are stepped, since nights with no allocation cannot
  change the state.
- What counts as observed on a clear night is the semester plan's assignment,
  not an executed TTP night plan, so losses within an otherwise clear night are
  not modeled.

## Next

More seeds, weighted toward heavier and later-concentrated losses, which is
where the endpoint would be expected to separate. Worth consolidating with the
note 02 harness first rather than maintaining two.
