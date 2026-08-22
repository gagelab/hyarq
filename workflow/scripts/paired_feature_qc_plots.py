import argparse
import sys
from dataclasses import dataclass

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import is_color_like, to_rgba
from matplotlib.ticker import MaxNLocator


PAIRED_FEATURE_RETENTION_THRESHOLDS = (1, 5, 10, 20, 30, 60)
SAMPLE_HEXBIN_GRIDSIZE = 45
SAMPLE_HEXBIN_CMAP = "viridis"


@dataclass(frozen=True)
class FeatureTerminology:
    singular: str
    plural: str
    plural_capitalized: str
    hyphenated: str


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
    plotted_count,
    histogram_colors,
    parent1_label,
    parent2_label,
    terminology,
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
            f"No {terminology.plural} have\na positive total count",
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
    ax.set_ylabel(terminology.plural_capitalized)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    subtitle = (
        f"{terminology.plural_capitalized} with total count > 0: "
        f"{plotted_count:,}"
    )
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
    terminology,
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
            f"No pooled {terminology.plural} have\na positive total count",
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
    ax.set_ylabel(terminology.plural_capitalized)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(
        "A. Parental count proportion\n"
        f"{terminology.plural_capitalized} with total count > 0: "
        f"{plotted_count:,}"
    )
    return plotted_count


def plot_parent_count_hexbin(
    ax,
    pooled,
    parent1_label,
    parent2_label,
    terminology,
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
        colorbar.set_label(
            f"{terminology.plural_capitalized} per hexagon (log scale)"
        )
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
            f"No pooled {terminology.plural} have\na positive total count",
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
        f"B. Parental counts per {terminology.singular}\n"
        f"{terminology.plural_capitalized} with total count > 0: "
        f"{plotted_count:,}"
    )
    return plotted_count


def plot_per_library_cumulative_depth_curves(
    ax,
    tables,
    depth_curve_color,
    depth_curve_alpha,
    terminology,
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
            f"No libraries have {terminology.plural} with total count > 0",
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
    ax.set_xlabel(f"log2(total count per {terminology.singular})")
    ax.set_ylabel(
        f"Cumulative fraction of counted {terminology.plural}"
    )
    ax.set_title(
        f"C. Per-library cumulative {terminology.hyphenated} depth curves\n"
        f"{terminology.plural_capitalized} with total count > 0"
    )
    return included_library_count, skipped_library_count


def plot_paired_feature_retention(
    ax,
    retention,
    retention_colors,
    terminology,
):
    positions = np.arange(
        1,
        len(PAIRED_FEATURE_RETENTION_THRESHOLDS) + 1,
    )
    distributions = []
    for threshold in PAIRED_FEATURE_RETENTION_THRESHOLDS:
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
            f"No {terminology.hyphenated} retention values available",
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

    ax.set_xlim(
        0.5,
        len(PAIRED_FEATURE_RETENTION_THRESHOLDS) + 0.5,
    )
    ax.set_ylim(0, 100)
    ax.set_xticks(positions)
    ax.set_xticklabels(PAIRED_FEATURE_RETENTION_THRESHOLDS)
    ax.set_xlabel(
        f"Minimum total count per {terminology.singular}"
    )
    ax.set_ylabel(f"Retained {terminology.plural} (%)")
    ax.set_title(
        f"D. Across-library {terminology.hyphenated} retention"
    )


def write_paired_feature_qc_report(
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
    terminology,
    report_prefix,
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
            terminology=terminology,
        )
        print(
            f"Panel A pooled positive-count {terminology.plural}: "
            f"{panel_a_plotted_count}",
            file=sys.stderr,
        )
        panel_b_plotted_count = plot_parent_count_hexbin(
            axes[0, 1],
            pooled,
            parent1_label,
            parent2_label,
            terminology,
        )
        print(
            f"Panel B pooled positive-count {terminology.plural}: "
            f"{panel_b_plotted_count}",
            file=sys.stderr,
        )
        panel_c_included_count, panel_c_skipped_count = (
            plot_per_library_cumulative_depth_curves(
                axes[1, 0],
                tables,
                depth_curve_color,
                depth_curve_alpha,
                terminology,
            )
        )
        print(
            f"Panel C libraries with {terminology.plural} having "
            "total count > 0: "
            f"{panel_c_included_count}; skipped libraries: {panel_c_skipped_count}",
            file=sys.stderr,
        )
        plot_paired_feature_retention(
            axes[1, 1],
            retention,
            retention_colors,
            terminology,
        )
        retention_library_count = retention["library_id"].nunique()
        print(
            "Panel D retention libraries: "
            f"{retention_library_count}; thresholds: "
            f"{','.join(map(str, PAIRED_FEATURE_RETENTION_THRESHOLDS))}",
            file=sys.stderr,
        )
        fig.suptitle(f"{report_prefix} QC")
        pdf.savefig(fig)
        plt.close(fig)
    print(
        f"{report_prefix} QC report pages written: 1",
        file=sys.stderr,
    )


def prepare_sample_scatter_arrays(table):
    plotted = table.loc[table["total_count"] > 0]
    x = np.log2(plotted["parent1_count"] + 1).to_numpy()
    y = np.log2(plotted["parent2_count"] + 1).to_numpy()
    return x, y


def plot_sample_parent_count_hexbin(
    fig,
    ax,
    library_id,
    table,
    parent1_label,
    parent2_label,
    terminology,
):
    x, y = prepare_sample_scatter_arrays(table)
    plotted_count = len(x)
    axis_max = 1.0
    if len(x) > 0:
        axis_max = max(
            axis_max,
            float(np.max(x)),
            float(np.max(y)),
        )
        extent = (
            0,
            axis_max,
            0,
            axis_max,
        )
        hexbin_collection = ax.hexbin(
            x,
            y,
            gridsize=SAMPLE_HEXBIN_GRIDSIZE,
            mincnt=1,
            bins="log",
            cmap=SAMPLE_HEXBIN_CMAP,
            linewidths=0,
            extent=extent,
            zorder=2,
        )
        colorbar = fig.colorbar(
            hexbin_collection,
            ax=ax,
            fraction=0.05,
            pad=0.02,
        )
        colorbar.set_label(
            f"{terminology.plural_capitalized} per hexagon (log scale)"
        )
    else:
        ax.text(
            0.5,
            0.5,
            f"No {terminology.plural} have\na positive total count",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 4},
            zorder=4,
        )

    ax.plot(
        [0, axis_max],
        [0, axis_max],
        linestyle="--",
        color="black",
        linewidth=1,
        zorder=3,
    )
    ax.set_xlim(0, axis_max)
    ax.set_ylim(0, axis_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"log2({parent1_label} count + 1)")
    ax.set_ylabel(f"log2({parent2_label} count + 1)")
    ax.set_title(
        f"{library_id}\n"
        f"{terminology.plural_capitalized} with total count > 0: "
        f"{plotted_count:,}",
        fontsize=9,
    )


def plot_sample_parental_fraction_histogram(
    ax,
    library_id,
    table,
    histogram_colors,
    parent1_label,
    parent2_label,
    terminology,
):
    plotted = table.loc[table["total_count"] > 0]
    plotted_count = len(plotted)
    fractions = plotted["parent1_count"] / plotted["total_count"]
    plot_parental_fraction_histogram(
        ax,
        fractions,
        plotted_count,
        histogram_colors=histogram_colors,
        parent1_label=parent1_label,
        parent2_label=parent2_label,
        terminology=terminology,
        panel_title=library_id,
    )


def write_paired_feature_qc_sample_report(
    library_ids,
    tables,
    path,
    plot_types,
    histogram_colors,
    parent1_label,
    parent2_label,
    terminology,
    report_prefix,
):
    sample_items = list(zip(library_ids, tables))
    if plot_types == ("histogram",):
        page_capacity = 8
    elif plot_types == ("scatter",):
        page_capacity = 8
    else:
        page_capacity = 4
    page_count = (
        len(sample_items) + page_capacity - 1
    ) // page_capacity

    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        for page_index in range(page_count):
            start = page_index * page_capacity
            page_items = sample_items[
                start:start + page_capacity
            ]

            if len(plot_types) == 1:
                if plot_types == ("histogram",):
                    figure_size = (16, 7)
                else:
                    figure_size = (20, 8)
                fig, axes = plt.subplots(
                    2,
                    4,
                    figsize=figure_size,
                    constrained_layout=True,
                    squeeze=False,
                )
                for ax, (library_id, table) in zip(
                    axes.flat,
                    page_items,
                ):
                    if plot_types == ("histogram",):
                        plot_sample_parental_fraction_histogram(
                            ax,
                            library_id,
                            table,
                            histogram_colors,
                            parent1_label,
                            parent2_label,
                            terminology,
                        )
                    else:
                        plot_sample_parent_count_hexbin(
                            fig,
                            ax,
                            library_id,
                            table,
                            parent1_label,
                            parent2_label,
                            terminology,
                        )
                for ax in axes.flat[len(page_items):]:
                    ax.set_visible(False)
                if plot_types == ("histogram",):
                    fig.suptitle(
                        f"{report_prefix} parental fractions"
                    )
                else:
                    fig.suptitle(
                        f"{report_prefix} parental counts"
                    )
            else:
                fig, axes = plt.subplots(
                    2,
                    4,
                    figsize=(20, 8),
                    constrained_layout=True,
                    squeeze=False,
                )
                used_positions = set()
                for item_index, (library_id, table) in enumerate(
                    page_items
                ):
                    row_index = item_index // 2
                    histogram_column = (
                        item_index % 2
                    ) * 2
                    hexbin_column = histogram_column + 1
                    histogram_ax = axes[
                        row_index,
                        histogram_column,
                    ]
                    hexbin_ax = axes[
                        row_index,
                        hexbin_column,
                    ]
                    used_positions.update({
                        (
                            row_index,
                            histogram_column,
                        ),
                        (
                            row_index,
                            hexbin_column,
                        ),
                    })
                    plot_sample_parental_fraction_histogram(
                        histogram_ax,
                        library_id,
                        table,
                        histogram_colors,
                        parent1_label,
                        parent2_label,
                        terminology,
                    )
                    plot_sample_parent_count_hexbin(
                        fig,
                        hexbin_ax,
                        library_id,
                        table,
                        parent1_label,
                        parent2_label,
                        terminology,
                    )
                for row_index in range(2):
                    for column_index in range(4):
                        if (
                            row_index,
                            column_index,
                        ) not in used_positions:
                            axes[
                                row_index,
                                column_index,
                            ].set_visible(False)
                fig.suptitle(f"{report_prefix} sample QC")

            pdf.savefig(fig)
            plt.close(fig)

    rendered_plot_types = ",".join(plot_types)
    print(
        f"{report_prefix} QC sample report: "
        f"libraries={len(sample_items)}; "
        f"plot types={rendered_plot_types}; "
        f"pages written={page_count}",
        file=sys.stderr,
    )
