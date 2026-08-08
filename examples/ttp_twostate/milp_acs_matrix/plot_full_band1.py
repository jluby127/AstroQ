"""Method comparison plots for the full-band1 matrix runs (58 targets).

Usage::

    python examples/ttp_twostate/milp_acs_matrix/plot_full_band1.py
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from plot_scaling import (  # noqa: E402
    ACS_METHODS,
    BASELINE,
    METHODS_120,
    METHODS_600,
    METHOD_ORDER,
    METHOD_STYLE,
)

RESULTS_CSV = os.path.join(HERE, "all_results.csv")
SUMMARY_DIR = os.path.join(HERE, "summary", "full_band1")
DATASET = "full-band1"

SHORT_LABEL = {
    "single_milp_120": "single\nMILP 120s",
    "single_milp_600": "single\nMILP 600s",
    "single_acs8_120": "single\nACS+MILP 120s",
    "single_acs8_600": "single\nACS+MILP 600s",
    "single_norel_120": "single\nNoRel+MILP 120s",
    "two_milp_120": "two-state\nMILP 120s",
    "two_milp_600": "two-state\nMILP 600s",
    "two_acs8_120": "two-state\nACS+MILP 120s",
    "two_acs8_600": "two-state\nACS+MILP 600s",
    "two_norel_120": "two-state\nNoRel+MILP 120s",
}


def _load_full_band1():
    if not os.path.isfile(RESULTS_CSV):
        raise SystemExit(f"Missing {RESULTS_CSV}; run run_matrix.py first")
    df = pd.read_csv(RESULTS_CSV)
    df = df[df["dataset"] == DATASET].copy()
    if df.empty:
        raise SystemExit(f"No rows for dataset={DATASET!r} in {RESULTS_CSV}")
    order = {m: i for i, m in enumerate(METHOD_ORDER)}
    df["_sort"] = df["method"].map(order)
    return df.sort_values("_sort")


def _colors(methods):
    return [METHOD_STYLE.get(m, ("#333333", "-"))[0] for m in methods]


def _baseline_value(df, metric, baseline=BASELINE):
    row = df.loc[df["method"] == baseline, metric]
    if row.empty:
        raise ValueError(f"Missing baseline method {baseline!r}")
    return float(row.iloc[0])


def _delta_values(df, methods, metric, *, baseline=BASELINE, invert=False):
    """Return per-method delta vs baseline (invert: positive = lower metric is better)."""
    base = _baseline_value(df, metric, baseline=baseline)
    sub = df[df["method"].isin(methods)].set_index("method").reindex(methods)
    vals = sub[metric].astype(float).values
    delta = vals - base
    if invert:
        delta = -delta
    return delta


def _plot_metric_bars(
    df,
    methods,
    metric,
    *,
    title,
    ylabel,
    out_name,
    lower_is_better=False,
    baseline=None,
    invert_delta=False,
):
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    sub = df[df["method"].isin(methods)].set_index("method").reindex(methods)
    values = sub[metric].astype(float).values
    if baseline is not None:
        if baseline in sub.index:
            base = float(sub.loc[baseline, metric])
        else:
            base = _baseline_value(df, metric, baseline=baseline)
        values = base - values if invert_delta else values - base

    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(methods)), 5))
    x = np.arange(len(methods))
    bars = ax.bar(x, values, color=_colors(methods), edgecolor="#333333", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT_LABEL.get(m, m) for m in methods], fontsize=8)
    ax.set_ylabel(ylabel)
    n_targets = int(sub["n_targets"].iloc[0])
    if baseline is not None:
        ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle=":", zorder=0)
        ax.set_title(f"{title}\n(relative to {baseline}; full-band1, N={n_targets})")
    else:
        ax.set_title(f"{title}\n(full-band1, N={n_targets})")
    if lower_is_better and baseline is None:
        best = np.nanmin(values)
        for bar, val in zip(bars, values):
            if np.isclose(val, best):
                bar.set_edgecolor("#000000")
                bar.set_linewidth(2.0)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(SUMMARY_DIR, out_name), dpi=150)
    plt.close(fig)


def _plot_delta_grid(df, methods, out_name, *, baseline=BASELINE):
    """Four-metric panel: deltas vs baseline plus raw MIP gap."""
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    panels = [
        ("objective", "Δ objective", "delta", False),
        ("physical_slew", "Δ physical slew (min)", "delta", False),
        ("modeled_slew", "Δ modeled slew (min)", "delta", False),
        ("gap_pct", "MIP gap (%)", "absolute", False),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    colors = _colors(methods)
    x = np.arange(len(methods))
    n_targets = int(df["n_targets"].iloc[0])
    sub = df[df["method"].isin(methods)].set_index("method").reindex(methods)
    for ax, (col, ylabel, mode, invert) in zip(axes.ravel(), panels):
        if mode == "absolute":
            vals = sub[col].astype(float).values
        else:
            vals = _delta_values(df, methods, col, baseline=baseline, invert=invert)
        ax.bar(x, vals, color=colors, edgecolor="#333333", linewidth=0.6)
        if mode == "delta":
            ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle=":", zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels([SHORT_LABEL.get(m, m) for m in methods], fontsize=7)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        f"full-band1 deltas relative to {baseline} (N={n_targets})",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(SUMMARY_DIR, out_name), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_grouped_runtime_delta(df, out_name, *, baseline=BASELINE):
    """Grouped 120s vs 600s deltas vs baseline, by method family."""
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    families = [
        ("single MILP", "single_milp_120", "single_milp_600"),
        ("single ACS+MILP", "single_acs8_120", "single_acs8_600"),
        ("two-state MILP", "two_milp_120", "two_milp_600"),
        ("two-state ACS+MILP", "two_acs8_120", "two_acs8_600"),
    ]
    panels = [
        ("objective", "Δ objective", False),
        ("physical_slew", "Δ physical slew (min)", False),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    idx = np.arange(len(families))
    width = 0.35
    n_targets = int(df["n_targets"].iloc[0])
    for ax, (col, ylabel, invert) in zip(axes, panels):
        base = _baseline_value(df, col, baseline=baseline)
        vals_120, vals_600 = [], []
        for _, m120, m600 in families:
            v120 = float(df.loc[df["method"] == m120, col].iloc[0]) - base
            v600 = float(df.loc[df["method"] == m600, col].iloc[0]) - base
            if invert:
                v120, v600 = -v120, -v600
            vals_120.append(v120)
            vals_600.append(v600)
        ax.bar(idx - width / 2, vals_120, width, label="120s", color="#56B4E9")
        ax.bar(idx + width / 2, vals_600, width, label="600s", color="#0072B2")
        ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle=":", zorder=0)
        ax.set_xticks(idx)
        ax.set_xticklabels([f[0] for f in families], fontsize=8)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.suptitle(
        f"full-band1: 120s vs 600s (relative to {baseline}, N={n_targets})",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(SUMMARY_DIR, out_name), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_runtime_grid(df, methods, runtime_label, out_name):
    """Four-metric panel for one time limit."""
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    sub = df[df["method"].isin(methods)].set_index("method").reindex(methods)
    panels = [
        ("objective", "Objective", False),
        ("physical_slew", "Physical slew (min)", True),
        ("modeled_slew", "Modeled slew (min)", True),
        ("gap_pct", "MIP gap (%)", True),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    colors = _colors(methods)
    x = np.arange(len(methods))
    for ax, (col, ylabel, lower_better) in zip(axes.ravel(), panels):
        vals = sub[col].astype(float).values
        bars = ax.bar(x, vals, color=colors, edgecolor="#333333", linewidth=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels([SHORT_LABEL.get(m, m) for m in methods], fontsize=7)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.3)
        if lower_better:
            best = np.nanmin(vals)
            for bar, val in zip(bars, vals):
                if np.isclose(val, best):
                    bar.set_edgecolor("#000000")
                    bar.set_linewidth(2.0)
    fig.suptitle(
        f"full-band1 method comparison ({runtime_label})",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(SUMMARY_DIR, out_name), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_acs_improvement(df, methods, runtime_label, out_name):
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    sub = df[df["method"].isin(methods)].set_index("method").reindex(methods)
    panels = [
        ("scheduled", "acs_scheduled", "Δ scheduled\n(MILP − ACS)"),
        ("objective", "acs_objective", "Δ objective\n(MILP − ACS)"),
        ("physical_slew", "acs_physical_slew", "Δ physical slew\n(ACS − MILP, min)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
    colors = _colors(methods)
    x = np.arange(len(methods))
    for ax, (final_col, acs_col, ylabel) in zip(axes, panels):
        if final_col == "physical_slew":
            delta = sub[acs_col] - sub[final_col]
        else:
            delta = sub[final_col] - sub[acs_col]
        vals = delta.astype(float).values
        ax.bar(x, vals, color=colors, edgecolor="#333333", linewidth=0.6)
        ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle=":", zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels([SHORT_LABEL.get(m, m) for m in methods], fontsize=7)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        f"MILP improvement over ACS warm start ({runtime_label}, full-band1)",
        fontsize=12,
        y=1.05,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(SUMMARY_DIR, out_name), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_grouped_runtime(df, out_name):
    """Grouped bars: 120s vs 600s for each method family."""
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    families = [
        ("single MILP", "single_milp_120", "single_milp_600"),
        ("single ACS+MILP", "single_acs8_120", "single_acs8_600"),
        ("two-state MILP", "two_milp_120", "two_milp_600"),
        ("two-state ACS+MILP", "two_acs8_120", "two_acs8_600"),
    ]
    panels = [
        ("objective", "Objective", False),
        ("physical_slew", "Physical slew (min)", True),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    idx = np.arange(len(families))
    width = 0.35
    for ax, (col, ylabel, lower_better) in zip(axes, panels):
        vals_120, vals_600 = [], []
        for _, m120, m600 in families:
            row120 = df.loc[df["method"] == m120, col]
            row600 = df.loc[df["method"] == m600, col]
            vals_120.append(float(row120.iloc[0]))
            vals_600.append(float(row600.iloc[0]))
        ax.bar(idx - width / 2, vals_120, width, label="120s", color="#56B4E9")
        ax.bar(idx + width / 2, vals_600, width, label="600s", color="#0072B2")
        ax.set_xticks(idx)
        ax.set_xticklabels([f[0] for f in families], fontsize=8)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.3)
        if lower_better:
            ax.set_title(f"{ylabel} (lower is better)")
        else:
            ax.set_title(f"{ylabel} (higher is better)")
    axes[0].legend(fontsize=8)
    fig.suptitle("full-band1: 120s vs 600s by method family", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(SUMMARY_DIR, out_name), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    df = _load_full_band1()
    n = int(df["n_targets"].iloc[0])
    print(f"full-band1: N={n}, methods={list(df['method'])}")

    _plot_metric_bars(
        df, METHOD_ORDER, "objective",
        title="Objective by method",
        ylabel="Objective",
        out_name="objective_by_method.png",
    )
    _plot_metric_bars(
        df, METHOD_ORDER, "physical_slew",
        title="Physical slew by method",
        ylabel="Physical slew (min)",
        out_name="physical_slew_by_method.png",
        lower_is_better=True,
    )
    _plot_metric_bars(
        df, METHOD_ORDER, "modeled_slew",
        title="Modeled slew by method",
        ylabel="Modeled slew (min)",
        out_name="modeled_slew_by_method.png",
        lower_is_better=True,
    )

    _plot_runtime_grid(df, METHODS_120, "120s", "method_panel_120s.png")
    _plot_runtime_grid(df, METHODS_600, "600s", "method_panel_600s.png")

    _plot_delta_grid(df, METHOD_ORDER, "delta_panel_vs_single_milp_120.png")

    for metric, title, ylabel, invert in (
        ("objective", "Δ objective", "Δ objective", False),
        ("physical_slew", "Δ physical slew", "Δ physical slew (min)", False),
        ("modeled_slew", "Δ modeled slew", "Δ modeled slew (min)", False),
    ):
        suffix = " (positive = better)" if invert else ""
        _plot_metric_bars(
            df,
            METHOD_ORDER,
            metric,
            title=title,
            ylabel=f"{ylabel}{suffix}",
            out_name=f"delta_{metric}_vs_single_milp_120.png",
            baseline=BASELINE,
            invert_delta=invert,
        )

    _plot_metric_bars(
        df,
        METHOD_ORDER,
        "gap_pct",
        title="MIP gap",
        ylabel="MIP gap (%)",
        out_name="gap_pct_by_method.png",
        lower_is_better=True,
    )

    acs_120 = [m for m in ACS_METHODS if m in METHODS_120]
    acs_600 = [m for m in ACS_METHODS if m in METHODS_600]
    _plot_acs_improvement(df, acs_120, "120s", "milp_improvement_acs_120s.png")
    _plot_acs_improvement(df, acs_600, "600s", "milp_improvement_acs_600s.png")

    _plot_grouped_runtime(df, "grouped_120s_vs_600s.png")
    _plot_grouped_runtime_delta(df, "grouped_120s_vs_600s_vs_single_milp_120.png")

    print(f"Wrote full-band1 plots to {SUMMARY_DIR}/")


if __name__ == "__main__":
    main()
