"""Compare fill shortfall between the balanced and unbalanced weather arms.

Reads the paired metrics recorded in data/03/ and writes

    plots/weather_fsf.pdf

Two panels, because the two readings of fill shortfall answer different
questions. The upper panel measures shortfall against maximum feasible fill
recomputed each night, which is the quantity the balance stage actually
minimizes. The lower panel measures it against the night-1 value, which is how
far a program has fallen below what looked achievable when the semester opened.

Usage:
    python tools/plot_weather_sim.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
DATA = os.path.join(PAPER, "data", "03")
PLOTS = os.path.join(PAPER, "plots")

SEED = 1
# The filler backstop has a nominal 600 hr award rather than allocated time, and
# C362 submitted no requests, so neither belongs in a shortfall comparison.
EXCLUDE = {"2026B_E475"}

ARMS = {
    "balance": dict(label="With balance stage", color="#1f77b4", marker="o"),
    "nobalance": dict(label="Without balance stage", color="#d62728", marker="s"),
}


def load(arm):
    path = os.path.join(DATA, f"metrics-seed{SEED}-{arm}.csv")
    d = pd.read_csv(path, parse_dates=["date"])
    d = d[(d.awarded_hr > 0) & (~d.program.isin(EXCLUDE))]
    grouped = d.groupby("night_index")
    return pd.DataFrame({
        "date": grouped.date.first(),
        "clear": grouped.clear.first(),
        "fsf": grouped.fsf.max(),
        "fsf_frozen": grouped.fsf_frozen.max(),
    })


def main():
    arms = {name: load(name) for name in ARMS}
    ref = arms["balance"]
    lost = ref.index[~ref.clear]

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 6.4), sharex=True)

    panels = [
        ("fsf", "Worst fill shortfall\n(vs. nightly maximum feasible fill)"),
        ("fsf_frozen", "Worst fill shortfall\n(vs. night-1 maximum feasible fill)"),
    ]
    for ax, (column, ylabel) in zip(axes, panels):
        for night in lost:
            ax.axvspan(night - 0.5, night + 0.5, color="0.88", zorder=0,
                       label="Night lost to weather" if night == lost[0] else None)
        for name, style in ARMS.items():
            ax.plot(arms[name].index, arms[name][column],
                    color=style["color"], marker=style["marker"], markersize=4,
                    linewidth=1.6, label=style["label"], zorder=3)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(alpha=0.3, linewidth=0.5)
        ax.margins(x=0.01)

    axes[0].legend(fontsize=8, loc="upper left", framealpha=0.95)
    axes[0].set_title(
        f"2026B HIRES queue under simulated weather (seed {SEED}, "
        f"{len(lost)} of {len(ref)} allocated nights lost)", fontsize=10)

    # Month labels are more legible than 29 irregular dates.
    ticks = [i for i, d in enumerate(ref.date)
             if i == 0 or d.month != ref.date.iloc[i - 1].month]
    axes[1].set_xticks(ticks)
    axes[1].set_xticklabels([ref.date.iloc[i].strftime("%b %-d") for i in ticks])
    axes[1].set_xlabel("Allocated night of the semester", fontsize=9)

    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    out = os.path.join(PLOTS, "weather_fsf.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=150)
    print(f"wrote {out}")

    for column, _ in panels:
        b, n = arms["balance"][column], arms["nobalance"][column]
        print(f"  {column:11} mean {b.mean():.3f} vs {n.mean():.3f} | "
              f"peak {b.max():.3f} vs {n.max():.3f} | "
              f"final {b.iloc[-1]:.3f} vs {n.iloc[-1]:.3f}")


if __name__ == "__main__":
    main()
