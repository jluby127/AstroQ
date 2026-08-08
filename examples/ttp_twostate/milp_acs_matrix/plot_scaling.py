"""Summary scaling plots from ``all_results.csv``.

Usage::

    python examples/ttp_twostate/milp_acs_matrix/plot_scaling.py
    PLOT_SPHERE_ONLY=1 python examples/ttp_twostate/milp_acs_matrix/plot_scaling.py
    PLOT_ALL=1 python examples/ttp_twostate/milp_acs_matrix/plot_scaling.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_CSV = os.path.join(HERE, "all_results.csv")

N_ORDER_ALL = [10, 20, 40, 50, 58, 75, 100, 125]
N_ORDER_SPHERE = [10, 20, 40, 50, 75, 100, 125]

# Paul Tol–style distinct colors; solid=120s, dashed=600s within each family.
METHOD_STYLE = {
    "single_milp_120": ("#0072B2", "-"),
    "single_milp_600": ("#0072B2", "--"),
    "single_acs8_120": ("#E69F00", "-"),
    "single_acs8_600": ("#E69F00", "--"),
    "single_norel_120": ("#56B4E9", "-"),
    "two_milp_120": ("#009E73", "-"),
    "two_milp_600": ("#009E73", "--"),
    "two_acs8_120": ("#CC79A7", "-"),
    "two_acs8_600": ("#CC79A7", "--"),
    "two_norel_120": ("#D55E00", "-"),
}

METHOD_ORDER = [
    "single_milp_120",
    "single_milp_600",
    "single_acs8_120",
    "single_acs8_600",
    "single_norel_120",
    "two_milp_120",
    "two_milp_600",
    "two_acs8_120",
    "two_acs8_600",
    "two_norel_120",
]

METHODS_120 = [m for m in METHOD_ORDER if m.endswith("_120")]
METHODS_600 = [m for m in METHOD_ORDER if m.endswith("_600")]
ACS_METHODS = [m for m in METHOD_ORDER if "acs8" in m]

BASELINE = "single_milp_120"


def _with_baseline_delta(df, ycol, *, baseline=BASELINE, n_order=None):
    base = df.loc[df["method"] == baseline, ["n_targets", ycol]].copy()
    if n_order is not None:
        base = base[base["n_targets"].isin(n_order)]
    base = base.rename(columns={ycol: "_baseline"}).drop_duplicates("n_targets")
    out = df.merge(base, on="n_targets", how="inner")
    out["delta"] = out[ycol] - out["_baseline"]
    return out


def _plot_delta_metric(
    df, ycol, title, ylabel, out_name, *, summary_dir, n_order, methods=None,
    baseline=BASELINE,
):
    os.makedirs(summary_dir, exist_ok=True)
    delta_df = _with_baseline_delta(df, ycol, n_order=n_order, baseline=baseline)
    fig, ax = plt.subplots(figsize=(10, 6))
    if methods is None:
        methods = [
            m for m in METHOD_ORDER
            if m in delta_df["method"].unique() and m != baseline
        ]
    else:
        methods = [m for m in METHOD_ORDER if m in methods and m != baseline]

    for method in methods:
        sub = delta_df[delta_df["method"] == method].copy()
        sub = sub[sub["n_targets"].isin(n_order)]
        if sub.empty:
            continue
        sub["n_sort"] = sub["n_targets"].map({n: i for i, n in enumerate(n_order)})
        sub = sub.sort_values("n_sort")
        color, ls = METHOD_STYLE.get(method, ("#333333", "-"))
        ax.plot(
            sub["n_targets"],
            sub["delta"],
            marker="o",
            label=method,
            color=color,
            linestyle=ls,
            linewidth=2.0 if ls == "-" else 1.8,
            markersize=6,
        )

    ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle=":", zorder=0)
    ax.set_xlabel("Number of targets (dimensionality)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\n(relative to {baseline})")
    ax.set_xticks(n_order)
    ax.set_xlim(min(n_order) - 5, max(n_order) + 5)
    ax.legend(fontsize=8, ncol=2, framealpha=0.95)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(summary_dir, out_name), dpi=150)
    plt.close(fig)


def _plot_metric(
    df, ycol, title, ylabel, out_name, *, summary_dir, n_order, methods=None
):
    os.makedirs(summary_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    if methods is None:
        methods = [m for m in METHOD_ORDER if m in df["method"].unique()]
    else:
        methods = [m for m in METHOD_ORDER if m in methods]

    for method in methods:
        sub = df[df["method"] == method].copy()
        sub = sub[sub["n_targets"].isin(n_order)]
        if sub.empty:
            continue
        sub["n_sort"] = sub["n_targets"].map({n: i for i, n in enumerate(n_order)})
        sub = sub.sort_values("n_sort")
        color, ls = METHOD_STYLE.get(method, ("#333333", "-"))
        ax.plot(
            sub["n_targets"],
            sub[ycol],
            marker="o",
            label=method,
            color=color,
            linestyle=ls,
            linewidth=2.0 if ls == "-" else 1.8,
            markersize=6,
        )

    ax.set_xlabel("Number of targets (dimensionality)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(n_order)
    ax.set_xlim(min(n_order) - 5, max(n_order) + 5)
    ax.legend(fontsize=8, ncol=2, framealpha=0.95)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(summary_dir, out_name), dpi=150)
    plt.close(fig)


def _plot_milp_improvement_over_acs(
    df,
    metric_final,
    metric_acs,
    title,
    ylabel,
    out_name,
    *,
    summary_dir,
    n_order,
    methods=None,
    higher_is_better=True,
):
    """Plot post-MILP gain relative to the ACS warm-start checkpoint."""
    os.makedirs(summary_dir, exist_ok=True)
    if methods is None:
        methods = [m for m in ACS_METHODS if m in df["method"].unique()]
    else:
        methods = [m for m in ACS_METHODS if m in methods]

    fig, ax = plt.subplots(figsize=(10, 6))
    for method in methods:
        sub = df[df["method"] == method].copy()
        sub = sub[sub["n_targets"].isin(n_order)]
        if sub.empty:
            continue
        if higher_is_better:
            sub["delta"] = sub[metric_final] - sub[metric_acs]
        else:
            sub["delta"] = sub[metric_acs] - sub[metric_final]
        sub["n_sort"] = sub["n_targets"].map({n: i for i, n in enumerate(n_order)})
        sub = sub.sort_values("n_sort")
        color, ls = METHOD_STYLE.get(method, ("#333333", "-"))
        ax.plot(
            sub["n_targets"],
            sub["delta"],
            marker="o",
            label=method,
            color=color,
            linestyle=ls,
            linewidth=2.0 if ls == "-" else 1.8,
            markersize=6,
        )

    ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle=":", zorder=0)
    ax.set_xlabel("Number of targets (dimensionality)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(n_order)
    ax.set_xlim(min(n_order) - 5, max(n_order) + 5)
    ax.legend(fontsize=8, ncol=2, framealpha=0.95)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(summary_dir, out_name), dpi=150)
    plt.close(fig)


def _plot_acs_vs_final(
    df, metric_final, metric_acs, title, out_name, *, summary_dir, n_order
):
    os.makedirs(summary_dir, exist_ok=True)
    seeded = [m for m in METHOD_ORDER if "acs8" in m and m in df["method"].unique()]
    if not seeded:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    for method in seeded:
        sub = df[df["method"] == method].copy()
        sub = sub[sub["n_targets"].isin(n_order)]
        sub["n_sort"] = sub["n_targets"].map({n: i for i, n in enumerate(n_order)})
        sub = sub.sort_values("n_sort")
        color, _ = METHOD_STYLE.get(method, ("#333333", "-"))
        ax.plot(
            sub["n_targets"],
            sub[metric_final],
            marker="o",
            label=f"{method} (MILP)",
            color=color,
            linestyle="-",
            linewidth=2.0,
            markersize=6,
        )
        ax.plot(
            sub["n_targets"],
            sub[metric_acs],
            marker="x",
            label=f"{method} (ACS)",
            color=color,
            linestyle="--",
            linewidth=1.8,
            markersize=6,
        )

    ax.set_xlabel("Number of targets")
    ax.set_ylabel(title)
    ax.set_title(f"{title}: ACS checkpoint vs post-MILP")
    ax.set_xticks(n_order)
    ax.set_xlim(min(n_order) - 5, max(n_order) + 5)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(summary_dir, out_name), dpi=150)
    plt.close(fig)


def _make_plots(df, *, summary_dir, n_order, title_suffix=""):
    suffix = title_suffix
    acs_df = df[df["acs_starts"] > 0].copy()
    acs_methods = [m for m in METHOD_ORDER if m in acs_df["method"].unique()]
    kw = dict(summary_dir=summary_dir, n_order=n_order)

    _plot_metric(
        df, "scheduled", f"Scheduled targets vs dimensionality{suffix}", "Scheduled",
        "scheduled_vs_n.png", **kw,
    )
    _plot_metric(
        df, "physical_slew", f"Physical slew vs dimensionality{suffix}",
        "Physical slew (min)", "physical_slew_vs_n.png", **kw,
    )
    _plot_metric(
        df, "modeled_slew", f"Modeled slew vs dimensionality{suffix}",
        "Modeled slew (min)", "modeled_slew_vs_n.png", **kw,
    )
    _plot_metric(
        df, "objective", f"Objective vs dimensionality{suffix}", "Objective",
        "objective_vs_n.png", **kw,
    )

    if not acs_df.empty:
        _plot_metric(
            acs_df, "acs_scheduled",
            f"ACS checkpoint: scheduled vs dimensionality{suffix}",
            "Scheduled (ACS)", "acs_scheduled_vs_n.png",
            methods=acs_methods, **kw,
        )
        _plot_metric(
            acs_df, "acs_objective",
            f"ACS checkpoint: objective vs dimensionality{suffix}",
            "Objective (ACS)", "acs_objective_vs_n.png",
            methods=acs_methods, **kw,
        )
        _plot_metric(
            acs_df, "acs_physical_slew",
            f"ACS checkpoint: physical slew vs dimensionality{suffix}",
            "Physical slew (min, ACS)", "acs_physical_slew_vs_n.png",
            methods=acs_methods, **kw,
        )
        _plot_acs_vs_final(
            df, "scheduled", "acs_scheduled", "Scheduled count",
            "scheduled_acs_vs_milp.png", **kw,
        )
        _plot_acs_vs_final(
            df, "objective", "acs_objective", "Objective",
            "objective_acs_vs_milp.png", **kw,
        )
        _plot_acs_vs_final(
            df, "physical_slew", "acs_physical_slew", "Physical slew (min)",
            "physical_slew_acs_vs_milp.png", **kw,
        )

        acs_improve_kw = dict(summary_dir=summary_dir, n_order=n_order)
        _plot_milp_improvement_over_acs(
            acs_df, "scheduled", "acs_scheduled",
            f"MILP improvement over ACS warm start: scheduled{suffix}",
            "Δ scheduled (MILP − ACS)",
            "milp_improvement_scheduled_vs_n.png",
            methods=acs_methods,
            **acs_improve_kw,
        )
        _plot_milp_improvement_over_acs(
            acs_df, "objective", "acs_objective",
            f"MILP improvement over ACS warm start: objective{suffix}",
            "Δ objective (MILP − ACS)",
            "milp_improvement_objective_vs_n.png",
            methods=acs_methods,
            **acs_improve_kw,
        )
        _plot_milp_improvement_over_acs(
            acs_df, "physical_slew", "acs_physical_slew",
            f"MILP improvement over ACS warm start: physical slew{suffix}",
            "Δ physical slew (ACS − MILP, min)",
            "milp_improvement_physical_slew_vs_n.png",
            higher_is_better=False,
            methods=acs_methods,
            **acs_improve_kw,
        )
        _plot_milp_improvement_over_acs(
            acs_df, "modeled_slew", "acs_modeled_slew",
            f"MILP improvement over ACS warm start: modeled slew{suffix}",
            "Δ modeled slew (ACS − MILP, min)",
            "milp_improvement_modeled_slew_vs_n.png",
            higher_is_better=False,
            methods=acs_methods,
            **acs_improve_kw,
        )
        acs_120 = [m for m in ACS_METHODS if m in METHODS_120]
        acs_600 = [m for m in ACS_METHODS if m in METHODS_600]
        for runtime, acs_rt, tag in (
            ("120s", acs_120, "120s"),
            ("600s", acs_600, "600s"),
        ):
            rt_kw = {**acs_improve_kw, "methods": acs_rt}
            _plot_milp_improvement_over_acs(
                acs_df, "scheduled", "acs_scheduled",
                f"MILP improvement over ACS warm start: scheduled ({runtime}){suffix}",
                "Δ scheduled (MILP − ACS)",
                f"milp_improvement_scheduled_vs_n_{tag}.png",
                **rt_kw,
            )
            _plot_milp_improvement_over_acs(
                acs_df, "objective", "acs_objective",
                f"MILP improvement over ACS warm start: objective ({runtime}){suffix}",
                "Δ objective (MILP − ACS)",
                f"milp_improvement_objective_vs_n_{tag}.png",
                **rt_kw,
            )
            _plot_milp_improvement_over_acs(
                acs_df, "physical_slew", "acs_physical_slew",
                f"MILP improvement over ACS warm start: physical slew ({runtime}){suffix}",
                "Δ physical slew (ACS − MILP, min)",
                f"milp_improvement_physical_slew_vs_n_{tag}.png",
                higher_is_better=False,
                **rt_kw,
            )

    _plot_delta_metric(
        df, "scheduled", "Δ scheduled", "Δ scheduled (targets)",
        "delta_scheduled_vs_n.png", **kw,
    )
    _plot_delta_metric(
        df, "scheduled", "Δ scheduled (120s runs)", "Δ scheduled (targets)",
        "delta_scheduled_vs_n_120s.png", methods=METHODS_120, **kw,
    )
    _plot_delta_metric(
        df, "scheduled", "Δ scheduled (600s runs)", "Δ scheduled (targets)",
        "delta_scheduled_vs_n_600s.png", methods=METHODS_600,
        baseline="single_milp_600", **kw,
    )
    _plot_delta_metric(
        df, "physical_slew", "Δ physical slew", "Δ physical slew (min)",
        "delta_physical_slew_vs_n.png", **kw,
    )
    _plot_delta_metric(
        df, "modeled_slew", "Δ modeled slew", "Δ modeled slew (min)",
        "delta_modeled_slew_vs_n.png", **kw,
    )
    _plot_delta_metric(
        df, "objective", "Δ objective", "Δ objective",
        "delta_objective_vs_n.png", **kw,
    )
    _plot_delta_metric(
        df, "objective", "Δ objective (120s runs)", "Δ objective",
        "delta_objective_vs_n_120s.png", methods=METHODS_120, **kw,
    )
    _plot_delta_metric(
        df, "objective", "Δ objective (600s runs)", "Δ objective",
        "delta_objective_vs_n_600s.png", methods=METHODS_600,
        baseline="single_milp_600", **kw,
    )
    _plot_delta_metric(
        df, "gap_pct", "Δ MIP gap", "Δ gap (percentage points)",
        "delta_gap_vs_n.png", **kw,
    )

    if not acs_df.empty:
        _plot_delta_metric(
            acs_df, "acs_scheduled",
            "Δ scheduled (ACS checkpoint)", "Δ scheduled (targets)",
            "delta_acs_scheduled_vs_n.png", methods=acs_methods, **kw,
        )
        _plot_delta_metric(
            acs_df, "acs_objective",
            "Δ objective (ACS checkpoint)", "Δ objective",
            "delta_acs_objective_vs_n.png", methods=acs_methods, **kw,
        )
        _plot_delta_metric(
            acs_df, "acs_physical_slew",
            "Δ physical slew (ACS checkpoint)", "Δ physical slew (min)",
            "delta_acs_physical_slew_vs_n.png", methods=acs_methods, **kw,
        )


def main():
    if not os.path.isfile(RESULTS_CSV):
        raise SystemExit(f"Missing {RESULTS_CSV}; run run_matrix.py first")

    df = pd.read_csv(RESULTS_CSV)
    plot_all = os.environ.get("PLOT_ALL", "").strip() in ("1", "true", "yes")
    sphere_only = not plot_all

    if sphere_only:
        df = df[df["dataset"].str.startswith("sphere")].copy()
        summary_dir = os.path.join(HERE, "summary")
        n_order = N_ORDER_SPHERE
        title_suffix = " (sphere benchmarks)"
        print(f"Sphere-only filter: N={sorted(df['n_targets'].unique())}")
    else:
        summary_dir = os.path.join(HERE, "summary")
        n_order = N_ORDER_ALL
        title_suffix = ""

    df = df[df["n_targets"].isin(n_order)].copy()

    _make_plots(df, summary_dir=summary_dir, n_order=n_order, title_suffix=title_suffix)
    print(f"Wrote summary plots to {summary_dir}/")


if __name__ == "__main__":
    main()
