import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
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


def plot_pooled_parental_fraction_histogram(
    ax,
    pooled,
    histogram_color,
    parent1_label,
    histogram_bins,
):
    plotted = pooled.loc[pooled["total_count"] > 0]
    plotted_count = len(plotted)

    if plotted_count > 0:
        fractions = plotted["parent1_count"] / plotted["total_count"]
        bins = np.linspace(0.0, 1.0, histogram_bins + 1)
        ax.hist(
            fractions,
            bins=bins,
            color=histogram_color,
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No pooled gene pairs have\na positive total count",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )

    ax.set_xlim(0, 1)
    ax.axvline(0.5, linestyle="--", color="black", linewidth=1)
    ax.set_xlabel(f"{parent1_label} count proportion")
    ax.set_ylabel("Gene pairs")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(
        f"A. Parental count proportion\nGene pairs with total count > 0: {plotted_count:,}"
    )
    return plotted_count


def plot_parent_count_hexbin(
    ax,
    pooled,
    parent1_label,
    parent2_label,
):
    plotted = pooled.loc[pooled["total_count"] > 0]
    plotted_count = len(plotted)

    if plotted_count > 0:
        x = np.log2(plotted["parent1_count"] + 1)
        y = np.log2(plotted["parent2_count"] + 1)
        axis_limit = 1.05 * max(x.max(), y.max())
        hexbin_collection = ax.hexbin(
            x,
            y,
            gridsize=45,
            mincnt=1,
            bins="log",
            cmap="viridis",
            linewidths=0,
            extent=(0, axis_limit, 0, axis_limit),
            zorder=2,
        )
        colorbar = ax.figure.colorbar(
            hexbin_collection,
            ax=ax,
            pad=0.02,
        )
        colorbar.set_label("Gene pairs per hexagon (log scale)")
        ax.plot(
            [0, axis_limit],
            [0, axis_limit],
            linestyle="--",
            color="black",
            linewidth=1,
            zorder=3,
        )
        ax.set_xlim(0, axis_limit)
        ax.set_ylim(0, axis_limit)
    else:
        ax.text(
            0.5,
            0.5,
            "No pooled gene pairs have\na positive total count",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"log2({parent1_label} count + 1)")
    ax.set_ylabel(f"log2({parent2_label} count + 1)")
    ax.set_title(
        f"B. Parental counts per gene pair\nGene pairs with total count > 0: {plotted_count:,}"
    )
    return plotted_count


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
    histogram_bins,
):
    page_capacity = rows * columns
    sample_page_count = (len(library_ids) + page_capacity - 1) // page_capacity
    total_page_count = 1 + sample_page_count

    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(
            2,
            2,
            figsize=(10, 8),
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

        panel_a_plotted_count = plot_pooled_parental_fraction_histogram(
            axes[0, 0],
            pooled,
            histogram_color=histogram_color,
            parent1_label=parent1_label,
            histogram_bins=histogram_bins,
        )
        print(
            f"Panel A pooled positive-count gene pairs: {panel_a_plotted_count}",
            file=sys.stderr,
        )
        panel_b_plotted_count = plot_parent_count_hexbin(
            axes[0, 1],
            pooled,
            parent1_label,
            parent2_label,
        )
        print(
            f"Panel B pooled positive-count gene pairs: {panel_b_plotted_count}",
            file=sys.stderr,
        )
        plot_total_count_imbalance_scatter(
            axes[1, 0],
            plotted,
            min_total_count,
            plotted_count,
            parent1_label,
        )
        plot_gene_pair_retention(
            axes[1, 1],
            pooled,
        )
        fig.suptitle("Pooled RNA gene QC")
        pdf.savefig(fig)
        plt.close(fig)

        sample_items = list(zip(library_ids, tables))
        for page_index in range(sample_page_count):
            start = page_index * page_capacity
            page_items = sample_items[start:start + page_capacity]
            page_columns = min(columns, len(page_items))
            page_rows = (len(page_items) + page_columns - 1) // page_columns
            fig, axes = plt.subplots(
                page_rows,
                page_columns,
                figsize=(page_columns * 4, page_rows * 3),
                constrained_layout=True,
                squeeze=False,
            )
            fig.suptitle("RNA gene parental fractions")
            for ax, item in zip(axes.flat, page_items):
                library_id, table = item
                plotted = table.loc[table["total_count"] >= min_total_count]
                plotted_count = len(plotted)
                if plotted_count == 0:
                    print(
                        f"WARNING: no gene pairs meet the selected minimum total count for {library_id}",
                        file=sys.stderr,
                    )
                fractions = plotted["parent1_count"] / plotted["total_count"]
                plot_parental_fraction_histogram(
                    ax,
                    fractions,
                    min_total_count,
                    plotted_count,
                    histogram_color=histogram_color,
                    parent1_label=parent1_label,
                    panel_title=library_id,
                )
            for ax in axes.flat[len(page_items):]:
                ax.set_visible(False)
            pdf.savefig(fig)
            plt.close(fig)
    print(
        f"sample histograms minimum total count: {min_total_count}; libraries processed: {len(library_ids)}; sample pages written: {sample_page_count}",
        file=sys.stderr,
    )
    print(f"RNA gene QC report pages written: {total_page_count}", file=sys.stderr)
