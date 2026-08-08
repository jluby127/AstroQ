# MILP/ACS dimensionality sweep

80 runs: 8 target sets × 10 solver configurations.

## Target sets

| Key | N | Source |
|-----|---|--------|
| full-band1 | 58 | `full_band1_2026-06-12/` |
| sphere10 … sphere125 | 10–125 | `sphere_datasets/sphere{N}_2026-02-01/` |

Generate sphere datasets::

    python examples/ttp_twostate/sphere_datasets/gen_sphere.py --n 50

## Methods

Each row uses a benchmark `method_key` and a production-style `solve_method` (`milp`, `acs8+milp`, or `norel+milp`).

| method_key | states | solve_method | Gurobi budget |
|------------|--------|--------------|---------------|
| single_milp_120/600 | 1 | milp | TimeLimit 120 / 600 |
| single_acs8_120/600 | 1 | acs8+milp | 8×30s ACS + TimeLimit 120 / 600 |
| single_norel_120 | 1 | norel+milp | NoRel 30s + TimeLimit 90 (120s nominal) |
| two_milp_120/600 | 2 | milp | TimeLimit 120 / 600 |
| two_acs8_120/600 | 2 | acs8+milp | 8×30s ACS + TimeLimit 120 / 600 |
| two_norel_120 | 2 | norel+milp | NoRel 30s + TimeLimit 90 (120s nominal) |

`norel+milp` matches production [`nplan.py`](../../../astroq/nplan.py): `warmstart_time=30`, `max_solve_time=120` → 30s NoRel heuristic, 90s MILP TimeLimit.

## Run

```bash
python examples/ttp_twostate/milp_acs_matrix/run_matrix.py
python examples/ttp_twostate/milp_acs_matrix/plot_scaling.py
python examples/ttp_twostate/milp_acs_matrix/render_gifs.py   # GIFs for existing runs
```

Run only the norel methods (skip per-run HTML plots)::

    MATRIX_METHODS=single_norel_120,two_norel_120 MATRIX_SKIP_PLOTS=1 python ...

Env vars: `MATRIX_DATASETS`, `MATRIX_METHODS`, `MATRIX_SKIP_PLOTS=1`, `MATRIX_FORCE=1`.

## Outputs

- `runs/{dataset}/{method}/` — schedule.csv, slew_path.html, **slew_animation.html + slew_animation.gif**, ladder.html, run.json
- `all_results.csv` / `all_results.json` — flat table
- `summary/*.png` — scaling plots

## run.json fields

Identifiers: `dataset`, `n_targets`, `method`, `solve_method`, `n_states`, `milp_s`, `acs_starts`.

NoRel (null for other methods): `norel_heur_time`, `milp_time_limit`.

ACS checkpoint (null for plain MILP / norel): `acs_wall_s`, `acs_objective`, `acs_scheduled`, `acs_physical_slew`, `acs_modeled_slew`.

Post-MILP: `scheduled`, `physical_slew`, `modeled_slew`, `objective`, `bound`, `gap_pct`, `milp_solve_s`, `nodes`.
