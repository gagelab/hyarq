import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import is_color_like
from matplotlib.ticker import MaxNLocator


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
    ax.set_xlabel("Parent1 fraction")
    ax.set_ylabel("Gene pairs")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    subtitle = f"Minimum total count: {min_total_count}; plotted gene pairs: {plotted_count}"
    if panel_title is not None:
        subtitle = f"{panel_title}\n{subtitle}"
    ax.set_title(subtitle, fontsize=9)


def write_pooled_histogram(pooled, path, min_total_count, histogram_color):
    plotted = pooled.loc[pooled["total_count"] >= min_total_count]
    plotted_count = len(plotted)
    print(
        f"pooled histogram minimum total count: {min_total_count}; plotted gene pairs: {plotted_count}",
        file=sys.stderr,
    )
    if plotted_count == 0:
        print("WARNING: no pooled gene pairs meet the selected minimum total count", file=sys.stderr)

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    fractions = plotted["parent1_count"] / plotted["total_count"]
    plot_parental_fraction_histogram(
        ax,
        fractions,
        min_total_count,
        plotted_count,
        histogram_color=histogram_color,
    )
    fig.suptitle("Pooled RNA gene parental fractions")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="pdf")
    plt.close(fig)


def write_sample_histograms(library_ids, tables, path, min_total_count, histogram_color, rows, columns):
    page_capacity = rows * columns
    page_count = (len(library_ids) + page_capacity - 1) // page_capacity
    sample_items = list(zip(library_ids, tables))
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        for page_index in range(page_count):
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
                    panel_title=library_id,
                )
            for ax in axes.flat[len(page_items):]:
                ax.set_visible(False)
            pdf.savefig(fig)
            plt.close(fig)
    print(
        f"sample histograms minimum total count: {min_total_count}; libraries processed: {len(library_ids)}; PDF pages written: {page_count}",
        file=sys.stderr,
    )
