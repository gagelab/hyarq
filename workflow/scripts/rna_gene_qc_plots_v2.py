import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import is_color_like
from matplotlib.ticker import LogFormatterMathtext, LogLocator, MaxNLocator, NullFormatter


def plot_color(value):
    if is_color_like(value):
        return value
    raise argparse.ArgumentTypeError("must be a valid matplotlib color")


def plot_parental_fraction_histogram(
    ax,
    fractions,
    min_total_count,
    plotted_count,
    histogram_color,
    parent1_label,
    panel_title=None,
):
    bins = [i / 20 for i in range(21)]
    if plotted_count > 0:
        ax.hist(fractions, bins=bins, color=histogram_color)
    else:
        ax.text(
            0.5,
            0.5,
            "No gene pairs meet the selected\nminimum total count",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )
    ax.set_xlim(0, 1)
    ax.axvline(0.5, linestyle="--", color="black", linewidth=1)
    ax.set_xlabel(f"{parent1_label} fraction")
    ax.set_ylabel("Gene pairs")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    subtitle = f"Minimum total count: {min_total_count}; plotted gene pairs: {plotted_count}"
    if panel_title is not None:
        subtitle = f"{panel_title}\n{subtitle}"
    ax.set_title(subtitle, fontsize=9)


def plot_parent_count_scatter(
    ax,
    plotted,
    min_total_count,
    plotted_count,
    parent1_label,
    parent2_label,
):
    both_positive = plotted.loc[(plotted["parent1_count"] > 0) & (plotted["parent2_count"] > 0)]
    both_positive_count = len(both_positive)
    if both_positive_count > 0:
        x = both_positive["parent1_count"]
        y = both_positive["parent2_count"]
        ax.scatter(x, y, s=16, alpha=0.5, edgecolors="none", rasterized=True)
        max_count = max(x.max(), y.max())
        axis_limit = 1.05 * max_count
        ax.set_xscale("log", base=10)
        ax.set_yscale("log", base=10)
        ax.xaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlim(1, axis_limit)
        ax.set_ylim(1, axis_limit)
        ax.plot([1, axis_limit], [1, axis_limit], linestyle="--", color="black", linewidth=1)
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.text(
            0.5,
            0.5,
            "No gene pairs have positive counts\nfor both parents",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    ax.set_xlabel(f"{parent1_label} count")
    ax.set_ylabel(f"{parent2_label} count")
    ax.set_title(
        f"{parent1_label} count vs {parent2_label} count\nMinimum total count: {min_total_count}; both-positive gene pairs: {both_positive_count}",
        fontsize=9,
    )


def plot_total_count_imbalance_scatter(
    ax,
    plotted,
    min_total_count,
    plotted_count,
    parent1_label,
):
    if plotted_count > 0:
        x = plotted["total_count"]
        parent1_fraction = plotted["parent1_count"] / plotted["total_count"]
        y = (parent1_fraction - 0.5).abs()
        ax.scatter(x, y, s=16, alpha=0.5, edgecolors="none", rasterized=True)
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(base=10))
        ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(0, 0.5)
    else:
        ax.text(
            0.5,
            0.5,
            "No gene pairs meet the selected\nminimum total count",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 0.5)
    ax.set_xlabel("Total count")
    ax.set_ylabel(f"Absolute parental imbalance\n|{parent1_label} fraction - 0.5|")
    ax.set_title(
        f"Total count vs absolute parental imbalance\nMinimum total count: {min_total_count}; plotted gene pairs: {plotted_count}",
        fontsize=9,
    )


def plot_gene_pair_retention(
    ax,
    pooled,
):
    total_gene_pairs = len(pooled)
    if total_gene_pairs == 0:
        ax.text(
            0.5,
            0.5,
            "No gene pairs available",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 100)
    else:
        thresholds = [1, 5, 10, 20, 50, 100, 200, 500, 1000]
        positions = range(len(thresholds))
        retained_percentages = [
            100 * (pooled["total_count"] >= threshold).sum() / total_gene_pairs
            for threshold in thresholds
        ]
        ax.plot(positions, retained_percentages, marker="o")
        ax.set_xticks(positions)
        ax.set_xticklabels(thresholds)
        ax.set_ylim(0, 100)
    ax.set_xlabel("Minimum total count")
    ax.set_ylabel("Retained gene pairs (%)")
    ax.set_title(
        f"Gene-pair retention across count thresholds\nAll pooled gene pairs: {total_gene_pairs}",
        fontsize=9,
    )


def write_rna_gene_qc_report(
    library_ids,
    tables,
    pooled,
    path,
    min_total_count,
    histogram_color,
    rows,
    columns,
    parent1_label,
    parent2_label,
):
    total_page_count = 1

    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(
            figsize=(6.5, 4.5),
            constrained_layout=True,
        )
        plotted = pooled.loc[pooled["total_count"] >= min_total_count]
        plotted_count = len(plotted)
        print(
            f"pooled plots minimum total count: {min_total_count}; plotted gene pairs: {plotted_count}",
            file=sys.stderr,
        )
        if plotted_count == 0:
            print("WARNING: no pooled gene pairs meet the selected minimum total count", file=sys.stderr)

        fractions = plotted["parent1_count"] / plotted["total_count"]
        plot_parental_fraction_histogram(
            ax,
            fractions,
            min_total_count,
            plotted_count,
            histogram_color=histogram_color,
            parent1_label=parent1_label,
            panel_title=(
                "Pooled gene-pair parental-fraction distribution\n"
                "Counts summed by gene pair across all libraries"
            ),
        )
        pdf.savefig(fig)
        plt.close(fig)
    print(f"RNA gene QC report pages written: {total_page_count}", file=sys.stderr)
