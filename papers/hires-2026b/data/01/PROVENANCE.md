# Provenance for finding 01

Backs [`../../notes/01-fsf-floor-and-u258-mismatch.md`](../../notes/01-fsf-floor-and-u258-mismatch.md).

## Runs

Both runs were driven by the same command from the repository root of
`AstroQ-testing`:

```bash
make -C ops/HIRES/2026B DATE=2026-08-01 BAND=band1 plan-semester
```

| | baseline | test |
|---|---|---|
| allocation | 29 blocks, 164.55 h | 27 blocks, 152.42 h |
| commit | `45776c77f8f280726a092deda10c52671b8e2ae5` | `e84bb3bb6628b14ae1e2813260217d33a5738401` |
| branch | `feature/balance-fill-shortfall-opus-5` | `test/u258-drop-two-nights` |
| tag | — | `exp/u258-drop-two-nights` |
| run directory | `AstroQ-testing/2026B/2026-08-01/band1_baseline/` | `AstroQ-hires26b/2026B/2026-08-01/band1/` |
| wall clock | 2026-07-31 20:47 to 20:58 | 2026-07-31 21:11 to 21:30 |

The two runs landed in different directories because `WORKDIR` comes from the
`CC_OUTPUT_PATH` environment variable, which pointed at `AstroQ-hires26b` for
the second run. The baseline directory is a manual copy made before the second
run, since `prep-run` truncates `astroq.log` and every stage overwrites
`outputs/`.

The test branch is a throwaway that is not intended to merge. It is tagged so
the commit stays reachable if the branch is deleted.

## Difference between the runs

The test commit changes exactly two files relative to the baseline commit:

- `astroq/queue/hirescps/prep.py` — drops Keck dates `2026-08-06` and
  `2026-09-18` from the U258 crossmatch, which removes both the pooled capacity
  and the matching award, because awards are derived from these same blocks.
- `ops/HIRES/common/config_template.ini` — sets `TimeLimit = 600` and
  `OutputFlag = 1` under `[semester.balance.gurobi]`.

The shortfall-stage budget also differs (900/300 versus 300/120) because the
template was edited between the two runs. This is unintended and is recorded as
a confound in the finding. See `config-gurobi-delta.txt`.

## Files

| file | description |
|------|-------------|
| `baseline-allocation.csv`, `test-allocation.csv` | `allocation.csv` as consumed by the scheduler; pooled capacity, no program tags |
| `baseline-keck-blocks.csv`, `test-keck-blocks.csv` | crossmatched Keck schedule with `ProjCode`, the source of both the allocation and the awards |
| `baseline-programs.csv`, `test-programs.csv` | awards and fill bounds, including `max_feasible_fill` from `astroq compute-max-fill` |
| `baseline-report-shortfall.txt`, `baseline-report-balance.txt` | per-program statistics tables after each stage |
| `test-report-shortfall.txt`, `test-report-balance.txt` | same, for the test run |
| `solve-lines-baseline.txt`, `solve-lines-test.txt` | one line per Gurobi solve: status, objective, bound, gap, runtime, plus the balance worst-shortfall transition |
| `u258-night-usage-baseline.txt`, `u258-night-usage-test.txt` | per-night allocated versus scheduled hours broken out by program |
| `config-gurobi-delta.txt` | diff of the `[semester*]` config sections between the runs |
| `night_usage.py` | regenerates the night-usage tables from a run directory |

Saved as `.txt` rather than `.log` deliberately: the repository `.gitignore`
has a bare `*.log` pattern that would silently exclude them.

## Regenerating

The night-usage tables can be rebuilt from either run directory, as long as it
still exists:

```bash
python night_usage.py /path/to/run_dir > u258-night-usage-<tag>.txt
```

The paper's `params.tex` is derived from the files in this directory:

```bash
cd ../.. && python tools/make_params.py
```

Nothing here depends on the run directories surviving. They are gitignored and
will be overwritten by the next `make`.
