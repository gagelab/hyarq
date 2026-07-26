import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import is_color_like, to_rgba
from matplotlib.ticker import MaxNLocator


RNA_RETENTION_THRESHOLDS = (1, 5, 10, 20, 30, 60)


def plot_color(value):
    if is_color_like(value):
        return value
    raise argparse.ArgumentTypeError("must be a valid matplotlib color")


def draw_two_color_parental_fraction_histogram(
    ax,
    fractions,
    bins,
    histogram_colors,
):
    counts, bin_edges = np.histogram(fractions, bins=bins)
    bin_widths = np.diff(bin_edges)
    bin_midpoints = bin_edges[:-1] + bin_widths / 2
    bar_colors = [
        histogram_colors[0] if midpoint < 0.5 else histogram_colors[1]
        for midpoint in bin_midpoints
    ]
    ax.bar(
        bin_edges[:-1],
        counts,
        width=bin_widths,
        align="edge",
        color=bar_colors,
        linewidth=0,
    )


def add_parental_shift_labels(
    ax,
    histogram_colors,
    parent1_label,
    parent2_label,
):
    ax.text(
        0.02,
        0.96,
        parent2_label,
        color=histogram_colors[0],
        fontsize=8,
        ha="left",
        va="top",
        transform=ax.transAxes,
    )
    ax.text(
        0.98,
        0.96,
        parent1_label,
        color=histogram_colors[1],
        fontsize=8,
        ha="right",
        va="top",
        transform=ax.transAxes,
    )


def plot_parental_fraction_histogram(
    ax,
    fractions,
    min_total_count,
    plotted_count,
    histogram_colors,
    parent1_label,
    parent2_label,
    panel_title=None,
):
    bins = [i / 20 for i in range(21)]
    if plotted_count > 0:
        draw_two_color_parental_fraction_histogram(
            ax,
            fractions,
            bins,
            histogram_colors,
        )
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
    add_parental_shift_labels(
        ax,
        histogram_colors,
        parent1_label,
        parent2_label,
    )
    ax.set_xlim(0, 1)
    ax.axvline(0.5, linestyle="--", color="black", linewidth=1)
    ax.set_xlabel(f"{parent1_label} count proportion")
    ax.set_ylabel("Gene pairs")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    subtitle = f"Minimum total count: {min_total_count}; plotted gene pairs: {plotted_count}"
    if panel_title is not None:
        subtitle = f"{panel_title}\n{subtitle}"
    ax.set_title(subtitle, fontsize=9)


def plot_pooled_parental_fraction_histogram(
    ax,
    pooled,
    histogram_colors,
    parent1_label,
    parent2_label,
    histogram_bins,
):
    plotted = pooled.loc[pooled["total_count"] > 0]
    plotted_count = len(plotted)

    if plotted_count > 0:
        fractions = plotted["parent1_count"] / plotted["total_count"]
        bins = np.linspace(0.0, 1.0, histogram_bins + 1)
        draw_two_color_parental_fraction_histogram(
            ax,
            fractions,
            bins,
            histogram_colors,
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

    add_parental_shift_labels(
        ax,
        histogram_colors,
        parent1_label,
        parent2_label,
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


def plot_per_library_cumulative_depth_curves(
    ax,
    tables,
    depth_curve_color,
    depth_curve_alpha,
):
    sorted_counts_by_library = []
    skipped_library_count = 0
    for table in tables:
        positive_total_counts = table.loc[
            table["total_count"] > 0,
            "total_count",
        ].to_numpy()
        if len(positive_total_counts) == 0:
            skipped_library_count += 1
            continue
        sorted_counts_by_library.append(np.sort(positive_total_counts))

    included_library_count = len(sorted_counts_by_library)
    if included_library_count > 0:
        maximum_positive_count = max(
            sorted_total_counts[-1]
            for sorted_total_counts in sorted_counts_by_library
        )
        maximum_log2_depth = np.log2(maximum_positive_count)
        grid_upper_log2 = max(maximum_log2_depth, 1.0)
        grid_point_count = 512
        log2_depth_grid = np.linspace(
            0.0,
            grid_upper_log2,
            grid_point_count,
        )
        raw_depth_grid = np.exp2(log2_depth_grid)
        if maximum_positive_count > 1:
            raw_depth_grid[-1] = maximum_positive_count

        library_curves = np.asarray([
            np.searchsorted(
                sorted_total_counts,
                raw_depth_grid,
                side="right",
            )
            / len(sorted_total_counts)
            for sorted_total_counts in sorted_counts_by_library
        ])
        median_curve = np.median(library_curves, axis=0)
        for curve in library_curves:
            ax.step(
                log2_depth_grid,
                curve,
                where="post",
                color=depth_curve_color,
                alpha=depth_curve_alpha,
                linewidth=0.8,
                zorder=1,
            )
        ax.step(
            log2_depth_grid,
            median_curve,
            where="post",
            color="black",
            alpha=1,
            linewidth=2.0,
            zorder=2,
        )
        ax.set_xlim(0, grid_upper_log2)
    else:
        ax.text(
            0.5,
            0.5,
            "No libraries have gene pairs with total count > 0",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )
        ax.set_xlim(0, 1)
    ax.xaxis.set_major_locator(
        MaxNLocator(nbins=6, integer=True)
    )
    ax.set_ylim(0, 1)
    ax.set_xlabel("log2(total count per gene pair)")
    ax.set_ylabel("Cumulative fraction of counted gene pairs")
    ax.set_title(
        "C. Per-library cumulative gene-pair depth curves\n"
        "Gene pairs with total count > 0"
    )
    return included_library_count, skipped_library_count


def plot_gene_pair_retention(
    ax,
    retention,
    retention_colors,
):
    positions = np.arange(1, len(RNA_RETENTION_THRESHOLDS) + 1)
    distributions = []
    for threshold in RNA_RETENTION_THRESHOLDS:
        values = (
            retention.loc[
                retention["minimum_total_count"] == threshold,
                "retained_percent",
            ]
            .dropna()
            .astype(float)
            .to_numpy()
        )
        distributions.append(values[np.isfinite(values)])

    if not any(len(values) > 0 for values in distributions):
        ax.text(
            0.5,
            0.5,
            "No gene-pair retention values available",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
        )
    else:
        for position, values, color in zip(
            positions,
            distributions,
            retention_colors,
        ):
            if len(values) == 0:
                continue
            ax.boxplot(
                [values],
                positions=[position],
                widths=0.55,
                patch_artist=True,
                showfliers=False,
                notch=False,
                manage_ticks=False,
                boxprops={
                    "facecolor": to_rgba(color, alpha=0.2),
                    "edgecolor": color,
                },
                whiskerprops={"color": color},
                capprops={"color": color},
                medianprops={"color": "black"},
            )
            if len(values) == 1:
                jitter = np.zeros(1)
            else:
                jitter = np.linspace(-0.12, 0.12, len(values))
            ax.scatter(
                position + jitter,
                values,
                color=color,
                marker="o",
                s=16,
                alpha=0.65,
                edgecolors="none",
                zorder=3,
            )

    ax.set_xlim(0.5, len(RNA_RETENTION_THRESHOLDS) + 0.5)
    ax.set_ylim(0, 100)
    ax.set_xticks(positions)
    ax.set_xticklabels(RNA_RETENTION_THRESHOLDS)
    ax.set_xlabel("Minimum total count per gene pair")
    ax.set_ylabel("Retained gene pairs (%)")
    ax.set_title("D. Across-library gene-pair retention")


def write_rna_gene_qc_report(
    tables,
    pooled,
    path,
    histogram_colors,
    parent1_label,
    parent2_label,
    histogram_bins,
    depth_curve_color,
    depth_curve_alpha,
    retention,
    retention_colors,
):
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(
            2,
            2,
            figsize=(10, 8),
            constrained_layout=True,
        )
        panel_a_plotted_count = plot_pooled_parental_fraction_histogram(
            axes[0, 0],
            pooled,
            histogram_colors=histogram_colors,
            parent1_label=parent1_label,
            parent2_label=parent2_label,
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
        panel_c_included_count, panel_c_skipped_count = (
            plot_per_library_cumulative_depth_curves(
                axes[1, 0],
                tables,
                depth_curve_color,
                depth_curve_alpha,
            )
        )
        print(
            "Panel C libraries with gene pairs having total count > 0: "
            f"{panel_c_included_count}; skipped libraries: {panel_c_skipped_count}",
            file=sys.stderr,
        )
        plot_gene_pair_retention(
            axes[1, 1],
            retention,
            retention_colors,
        )
        retention_library_count = retention["library_id"].nunique()
        print(
            "Panel D retention libraries: "
            f"{retention_library_count}; thresholds: "
            f"{','.join(map(str, RNA_RETENTION_THRESHOLDS))}",
            file=sys.stderr,
        )
        fig.suptitle("RNA gene QC")
        pdf.savefig(fig)
        plt.close(fig)
    print("RNA gene QC report pages written: 1", file=sys.stderr)
