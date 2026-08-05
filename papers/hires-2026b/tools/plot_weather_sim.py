"""Compare fill shortfall between the balanced and unbalanced weather arms.

Reads the paired metrics recorded in data/03/ and writes

    plots/weather_fsf.pdf           worst shortfall over all programs
    plots/weather_fsf_programs.pdf  the same, broken out per program

The summary figure has two panels, because the two readings of fill shortfall
answer different questions. The upper panel measures shortfall against maximum
feasible fill recomputed each night, which is the quantity the balance stage
actually minimizes. The lower panel measures it against the night-1 value,
which is how far a program has fallen below what looked achievable when the
semester opened.

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


def read_arm(arm):
    """Per-(night, program) rows for the science programs of one arm."""
    path = os.path.join(DATA, f"metrics-seed{SEED}-{arm}.csv")
    d = pd.read_csv(path, parse_dates=["date"])
    return d[(d.awarded_hr > 0) & (~d.program.isin(EXCLUDE))]


def load(arm):
    grouped = read_arm(arm).groupby("night_index")
    return pd.DataFrame({
        "date": grouped.date.first(),
        "clear": grouped.clear.first(),
        "fsf": grouped.fsf.max(),
        "fsf_frozen": grouped.fsf_frozen.max(),
    })


def month_ticks(ax, dates, fmt="%b %-d"):
    """Label the first allocated night of each month.

    Spacing is by allocated night rather than calendar date, so the labels mark
    where months begin without implying uniform time between points.
    """
    ticks = [i for i, d in enumerate(dates)
             if i == 0 or d.month != dates.iloc[i - 1].month]
    ax.set_xticks(ticks)
    ax.set_xticklabels([dates.iloc[i].strftime(fmt) for i in ticks])


def plot_programs(ref, lost):
    """Per-program shortfall, one panel per program."""
    frames = {name: read_arm(name) for name in ARMS}
    programs = sorted(frames["balance"].program.unique())
    awards = (frames["balance"].groupby("program").awarded_hr.first())

    ncols = 3
    nrows = -(-len(programs) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(9.5, 2.1 * nrows),
                             sharex=True, sharey=True)
    flat = axes.ravel()

    for ax, program in zip(flat, programs):
        for night in lost:
            ax.axvspan(night - 0.5, night + 0.5, color="0.88", zorder=0,
                       label="Night lost" if night == lost[0] else None)
        for name, style in ARMS.items():
            d = frames[name]
            d = d[d.program == program].set_index("night_index")
            ax.plot(d.index, d.fsf, color=style["color"], linewidth=1.5,
                    marker=style["marker"], markersize=3,
                    label=style["label"], zorder=3)
        ax.set_title(f"{program.replace('2026B_', '')} "
                     f"({awards[program]:.0f} hr)", fontsize=9)
        ax.grid(alpha=0.3, linewidth=0.5)
        ax.margins(x=0.01)

    for ax in flat[len(programs):]:
        ax.set_visible(False)

    flat[0].legend(fontsize=7, loc="upper left", framealpha=0.95)
    for ax in axes[-1]:
        if ax.get_visible():
            # Panels are narrow, so month alone rather than month and day.
            month_ticks(ax, ref.date, fmt="%b")
            ax.set_xlabel("Allocated night", fontsize=8)
            ax.tick_params(axis="x", labelsize=8)
    for row in axes:
        row[0].set_ylabel("Fill shortfall", fontsize=8)

    fig.suptitle(
        "Per-program fill shortfall vs. nightly maximum feasible fill "
        f"(seed {SEED}, {len(lost)} of {len(ref)} nights lost)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(PLOTS, "weather_fsf_programs.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=150)
    print(f"wrote {out}")
    return frames, programs


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

    month_ticks(axes[1], ref.date)
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

    frames, programs = plot_programs(ref, lost)

    print("\nper program, mean fill shortfall over the semester and final fill")
    print(f"  {'program':8} {'mean fsf bal':>13} {'nobal':>8} "
          f"{'final fill bal':>15} {'nobal':>8}")
    for program in programs:
        b = frames["balance"]
        b = b[b.program == program].sort_values("night_index")
        n = frames["nobalance"]
        n = n[n.program == program].sort_values("night_index")
        print(f"  {program.replace('2026B_', ''):8} {b.fsf.mean():13.3f} "
              f"{n.fsf.mean():8.3f} {b.fill.iloc[-1]:15.2f} "
              f"{n.fill.iloc[-1]:8.2f}")


if __name__ == "__main__":
    main()
