# HIRES 2026B paper notes

Running record of results from operating the California Planet Search HIRES
queue under AstroQ during semester 2026B. These notes exist so the paper can
be drafted later from evidence rather than memory.

## Working thesis

Moving an automated cadence scheduler from simulation onto a real queue
exposes a class of problem the methods paper could not: the queue's *structure*
determines whether programs compete at all. On paper the 2026B queue is
oversubscribed. In practice programs never contend for the same time, because
awards are derived from the allocation and fill is capped at the award. What
looks like a balancing problem turns out to be an allocation problem.

## Findings

| # | Title | Status |
|---|-------|--------|
| 01 | [Fill shortfall has a floor, and U258's seasonal mismatch was the only way past it](01-fsf-floor-and-u258-mismatch.md) | recorded |

## Conventions

**Every measured number is an artifact, not a memory.** Anything quoted in a
finding must be traceable to a file under `../data/<nn>/`. The run directories
that produced them are gitignored and get overwritten by the next `make`, so a
number that is not copied here is effectively lost.

**Numbers reach the paper through `\param{}`.** `../tools/make_params.py` reads
`../data/<nn>/` and generates `../params.tex`. The prose in `../main.tex` cites
`\param{base-u258-maff}`, never `0.71`. An unresolved key renders as a red XX
in the PDF, so a stale or missing number is visible rather than silent. Run
`make params` after adding artifacts, and `make check` before circulating a
draft.

**Log excerpts must be saved as `.txt`.** The repository `.gitignore` has a
bare `*.log` pattern, so an artifact named `astroq-excerpt.log` would be
silently untracked.

**Record the confounds in the finding itself.** Several of these results come
from operational runs rather than controlled experiments, and the differences
between runs are rarely limited to the one variable of interest. Each finding
carries an explicit section listing what else changed.

**Provenance lives with the data.** Each `../data/<nn>/PROVENANCE.md` records
the git SHA, branch, run directory, timestamp, and command for every run the
finding depends on. Where a finding depends on a throwaway branch, that commit
is tagged so the SHA stays reachable.

## Layout

```
notes/          findings, one numbered markdown file each
data/<nn>/      artifacts backing finding <nn>, plus PROVENANCE.md
tools/          make_params.py, generates params.tex from data/
main.tex        the paper; prose cites \param{} keys only
params.tex      GENERATED, do not hand-edit
```
