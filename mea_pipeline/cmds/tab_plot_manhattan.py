#!/usr/bin/env mea-pl

from mea_pipeline import cout, cerr, arg_parser


def init_argparser():
    p = arg_parser("Create multiple manhattan plot from columnar data")
    p.add_argument("--column", default=[], action="append")
    p.add_argument("--title", default=[], action="append")
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--dotsize", type=float, default=0.25)
    p.add_argument("--autoyscale", default=False, action="store_true")
    p.add_argument("--hline", default=False, action="store_true")
    p.add_argument("--bedfile", default=None, help="bed file to use as highlight")
    p.add_argument("--y-label", default=None)
    p.add_argument("--toptext", default=None)
    p.add_argument("--bottomtext", default=None)
    p.add_argument("-o", "--outfile", default="outplot.png")

    p.add_argument("infile")

    return p


def tab_plot_manhattan(args):

    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt

    cerr(f"Reading {args.infile}")
    if args.infile.endswith(".feather"):
        df = pd.read_feather(args.infile)
    else:
        df = pd.read_table(args.infile)

    df.sort_values(by=["CHROM", "POS"])
    # df['-log10 PValue'] = - np.log10(df['PValue'])
    # df[args.column[0]] = df[args.column[0]].abs()

    bed = None
    if args.bedfile:
        bed = pd.read_table(args.bedfile, header=None)

    ax = manhattan_plot(
        df,
        args.column,
        hline=args.hline,
        y_label=args.y_label,
        toptext=args.toptext,
        bottomtext=args.bottomtext,
        highlight=bed,
    )
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=600)


color_list_xx = [
    "#1f78b4",
    "#33a02c",
    "#e31a1c",
    "#ff7f00",
    "#6a3d9a",
]

color_list = [
    "#a6cee3",
    "#fb9a99",
    "#cab2d6",
    "#b2df8a",
    "#fdbf6f",
]

hicolor_list = [
    # hicolor
    "#1f78b4",
    "#e31a1c",
    "#6a3d9a",
    "#33a02c",
    "#ff7f00",
]

unused = [
    # unused
    "#b15928",
    "#ffff99",
]


def plot_scatter(data, column):

    from matplotlib import pyplot as plt

    paths = plt.scatter(
        x=data.index,
        y=data.loc[:, column],
        c=data["COLOR"],
        s=(data.loc[:, column] ** 2).abs().clip(lower=0.25),
    )


def manhattan_plot(
    data,
    column,
    y_label=None,
    hline: bool | float = False,
    dpi=300,
    toptext=None,
    bottomtext=None,
    highlight: list | None = None,
):
    """data=[CHROM, POS, value]"""

    import matplotlib.pyplot as plt
    import numpy as np
    from itertools import cycle
    import IPython

    N = data.shape[0]
    colors = cycle(color_list)
    hicolors = cycle(hicolor_list)

    # get unique values of first column of data
    if "CHROM_NO" in data.columns:
        chrom_no = data.CHROM_NO
    else:
        chrom_no = data.CHROM
    chroms = chrom_no.unique()
    data["COLOR"] = ""
    for chrom in chroms:
        cerr(f"Coloring {chrom}")
        c = next(colors)
        hi_c = next(hicolors)
        data.loc[chrom_no == chrom, "COLOR"] = c
        data.loc[(chrom_no == chrom) & (data.SIGNAL == 1), "COLOR"] = hi_c

    cerr("Plotting")
    fig = plt.figure(figsize=(15, 7), dpi=dpi)

    # plot background using plt.axvspan(x1, x2)
    highlight_text = []
    if highlight is not None:
        for idx, row in highlight.iterrows():
            regions = data.loc[
                (data.CHROM == row[0]) & (data.POS > row[1]) & (data.POS < row[2]), :
            ]
            if len(regions) > 0:
                plt.axvspan(regions.index[0], regions.index[1], color="lightgrey")
            highlight_text.append((row[3], regions.index[0]))

    #  do double coloring, non-signal first then signals so that signals would above non-signals
    plot_scatter(data.loc[data.SIGNAL == 0, :], column)
    plot_scatter(data.loc[data.SIGNAL == 1, :], column)

    if hline is True:
        if "UPPER_THRESHOLD" in data.columns:
            plt.axhline(data.UPPER_THRESHOLD[0], linestyle="--")
        if "LOWER_THRESHOLD" in data.columns:
            plt.axhline(data.LOWER_THRESHOLD[0], linestyle="--")
    elif type(hline) == float:
        plt.axhline(hline, linestyle="--")

    # chromosome number legend
    for chrom in chroms:
        chrom_pos = data.loc[chrom_no == chrom, :]
        min_idx = chrom_pos.index[0]
        max_idx = chrom_pos.index[-1]
        plt.text(
            (max_idx + min_idx) / 2,
            0,
            str(chrom),
            horizontalalignment="center",
            verticalalignment="center",
        )

    for text, hpos in highlight_text:
        plt.text(
            hpos,
            data[column].min() - 0.2,
            text,
            fontsize="xx-small",
            fontstretch="ultra-condensed",
            rotation="vertical",
            verticalalignment="baseline",
            horizontalalignment="right",
            color="grey",
        )

    # inner legend
    if toptext:
        plt.text(
            (data.index[0] + data.index[-1]) / 2,
            data[column].max() + 0.1,
            toptext,
            horizontalalignment="center",
            verticalalignment="center",
        )
    if bottomtext:
        plt.text(
            (data.index[0] + data.index[-1]) / 2,
            data[column].min() - 0.1,
            bottomtext,
            horizontalalignment="center",
            verticalalignment="center",
        )
    plt.ylabel(y_label or column)
    plt.xticks([])

    return

    # create a copy of the data and sort it by chrom and then position
    array = np.array(chr_pos_pvalue_array)
    if plot_threshold:
        array = array[array[:, 2] <= plot_threshold]
    else:
        plot_threshold = 1.0
    array = array[np.argsort(array[:, 1]), :]  # sort by ChrPos
    array = array[
        np.argsort(array[:, 0], kind="mergesort"), :
    ]  # Finally, sort by Chr (but keep ChrPos in case of ties)
    rle = list(_run_length_encode(array[:, 0]))

    if xaxis_unit_bp:  # compute and use cumulative basepair positions for x-axis
        if chromosome_starts is None:
            chromosome_starts = _compute_x_positions_chrom(array)
        chr_pos_list = _compute_x_positions_snps(array, chromosome_starts)
        plt.xlim([0, chromosome_starts[-1, 2] + 1])
        plt.xticks(
            chromosome_starts[:, 1:3].mean(1), chromosome_starts[:, 0].astype(int)
        )
    else:  # use rank indices for x-axis
        chr_pos_list = np.arange(array.shape[0])
        xTickMarks = [str(int(item)) for item, count in rle]
        plt.xlim([0, array.shape[0]])
        plt.xticks(list(_rel_to_midpoint(rle)), xTickMarks)
    y = -np.log10(array[:, 2])
    max_y = y.max()

    if (
        pvalue_line and vline_significant
    ):  # mark significant associations (ones that pass the pvalue_line) by a red vertical line:
        idx_significant = array[:, 2] < pvalue_line
        if np.any(idx_significant):
            y_significant = y[idx_significant]
            chr_pos_list_significant = chr_pos_list[idx_significant]
            for i in range(len(chr_pos_list_significant)):
                plt.axvline(
                    x=chr_pos_list_significant[i],
                    ymin=0.0,
                    ymax=y_significant[i],
                    color="r",
                    alpha=0.8,
                )

    plt.scatter(
        chr_pos_list,
        y,
        marker=marker,
        c=_color_list(array[:, 0], rle),
        edgecolor="none",
        s=y / max_y * 20 + 0.5,
        alpha=alpha,
    )
    plt.xlabel("chromosome")
    plt.ylabel("-log10(P value)")

    if pvalue_line:
        plt.axhline(-np.log10(pvalue_line), linestyle="--", color="gray")
    plt.ylim([-np.log10(plot_threshold), None])
    return chromosome_starts


def _compute_x_positions_chrom(positions, offset=1e5):
    chromosomes = np.unique(positions[:, 0])
    chromosomes.sort()
    chromosome_starts = np.zeros((chromosomes.shape[0], 3), dtype="float")
    chr_start_next = 0
    for i, chromosome in enumerate(chromosomes):
        pos_chr = positions[positions[:, 0] == chromosome, 1]
        chromosome_starts[i, 0] = chromosome  # the chromosome
        chromosome_starts[i, 1] = chr_start_next  # start of the chromosome
        # end of the chromosome
        chromosome_starts[i, 2] = chr_start_next + np.nanmax(pos_chr)
        chr_start_next = chromosome_starts[i, 2] + offset
    return chromosome_starts


def _compute_x_positions_snps(positions, chromosome_starts):
    cumulative_pos = np.zeros(positions.shape[0])
    for i, chromosome_start in enumerate(chromosome_starts):
        idx_chr = positions[:, 0] == chromosome_start[0]
        cumulative_pos[idx_chr] = positions[idx_chr][:, 1] + chromosome_start[1]
    return cumulative_pos


def main(args):
    tab_plot_manhattan(args)


# EOF
