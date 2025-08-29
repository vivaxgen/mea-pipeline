#!/usr/bin/env mea-pl

from mea_pipeline import cout, cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("Generate distribution plot from tabular data")
    p.add_argument("-t", "--title", default="", help="add title to plot")
    p.add_argument("-o", "--outplot", default="outdist.png", help="Output filename")
    p.add_argument("--dpi", default=600, type=int, help="resolution in dpi")
    p.add_argument("--add-hline", type=float, default=None)
    p.add_argument("--add-vline", type=float, default=None)
    p.add_argument(
        "--kind",
        default="hist",
        choices=["hist", "kde", "ecdf"],
        help="kind of distribution plot",
    )
    p.add_argument(
        "--use-y-axis",
        default=False,
        action="store_true",
        help="use y-axis for values, suitable for ecdf",
    )
    p.add_argument("--hue", default=None, help="column to use for grouping")
    p.add_argument("--add-Q1", action="store_true", default=False)
    p.add_argument("--add-Q2", action="store_true", default=False)
    p.add_argument("--add-Q3", action="store_true", default=False)
    p.add_argument("--set-ymin", default=None, type=float)
    p.add_argument("--drop-zero", default=False, action="store_true")
    p.add_argument("infile")
    return p


def add_vline(ax, x, label=None):
    ax.axvline(
        x,
        color="r",
        linestyle="--",
        linewidth=0.25,
    )

    if not label:
        return

    current_ticks = ax.get_xticks()
    current_labels = ax.get_xticklabels()
    new_ticks = list(current_ticks) + [x]
    new_labels = list(current_labels) + [label]
    ax.set_xticks(new_ticks)
    ax.set_xticklabels(new_labels)


def add_hline(ax, y, label=None):
    ax.axhline(
        y,
        color="r",
        linestyle="--",
        linewidth=0.25,
    )

    if not label:
        return

    current_ticks = ax.get_yticks()
    current_labels = ax.get_yticklabels()
    new_ticks = list(current_ticks) + [y]
    new_labels = list(current_labels) + [label]
    ax.set_yticks(new_ticks)
    ax.set_yticklabels(new_labels)


def displot(args):

    import seaborn as sns
    import pandas as pd
    from matplotlib import pyplot as plt
    from mea_pipeline.color_palettes import color_palettes

    infile, column = args.infile.split(":")
    cerr(f"Reading file: {args.infile}")
    df = pd.read_table(infile, sep="\t")

    if args.drop_zero:
        cerr(f"Dropping zero values for column {column}")
        df.drop(df.loc[df[column] == 0.0].index, inplace=True)
    cerr(f"Generating distribution plot with {df.shape[0]} data points")

    if args.hue:
        df.sort_values(by=args.hue, inplace=True)

    if args.use_y_axis:
        fg = sns.displot(
            data=df,
            y=column,
            kind=args.kind,
            hue=args.hue,
            palette=(
                sns.color_palette(color_palettes["xgfs_normal12"]) if args.hue else None
            ),
        )
    else:
        fg = sns.displot(
            data=df,
            x=column,
            kind=args.kind,
            hue=args.hue,
            palette=(
                sns.color_palette(color_palettes["xgfs_normal12"]) if args.hue else None
            ),
        )
    ax = fg.ax
    if args.kind == "ecdf":
        if args.add_Q1:
            add_vline(ax, 0.25, "Q1")
        if args.add_Q2:
            add_vline(ax, 0.50, "Q2")
        if args.add_Q3:
            add_vline(ax, 0.75, "Q3")
    if args.add_vline:
        add_vline(ax, args.add_vline, str(args.add_vline))
    if args.add_hline:
        add_hline(ax, args.add_hline, str(args.add_hline))
    ax.set_ylim(ymax=1.0)
    if args.set_ymin is not None:
        ax.set_ylim(ymin=args.set_ymin)
    if args.title:
        ax.set(title=args.title)
    if args.hue:
        sns.move_legend(fg, "upper left", bbox_to_anchor=(0.5, 0))
    plt.tight_layout()

    plt.savefig(
        args.outplot,
        metadata={"Creator": f"mea-pl plt-distribution to {args.outplot}"},
        dpi=args.dpi,
        bbox_inches="tight",
    )


def main(args):
    displot(args)


# EOF
