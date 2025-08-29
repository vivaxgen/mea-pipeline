#!/usr/bin/env mea-pl

from mea_pipeline import cout, cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("Generate categorical scatter plot from tabular data")
    p.add_argument("-t", "--title", default="", help="add title to plot")
    p.add_argument("-o", "--outplot", default="catplot.png", help="Output filename")
    p.add_argument("--dpi", default=600, type=int, help="resolution in dpi")

    p.add_argument("--add-hline", default=None, type=float)
    p.add_argument("--add-means", default=False, action="store_true")
    p.add_argument("--set-ymin", default=None, type=float)

    p.add_argument(
        "--set-label",
        default=None,
        help="set label for x-axis if the data is a single population",
    )
    p.add_argument("--hue", default=None, help="column to use for grouping")
    p.add_argument("--size", default=5, type=float)
    p.add_argument("--drop-zero", default=False, action="store_true")
    p.add_argument("--drop-one", default=False, action="store_true")
    p.add_argument("infile")

    return p


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


def tab_plot_categories(args):

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

    if args.drop_one:
        cerr(f"Dropping one values for column {column}")
        df.drop(df.loc[df[column] == 1.0].index, inplace=True)

    cerr(f"Generating distribution plot with {df.shape[0]} data points")

    if args.set_label and args.hue:
        cexit("ERROR: --set-label and --hue are mutually exclusive")

    if args.hue:
        df.sort_values(by=args.hue, inplace=True)

    fg = sns.catplot(
        data=df,
        x=args.hue,
        hue=args.hue,
        y=column,
        palette=(
            sns.color_palette(color_palettes["xgfs_normal12"]) if args.hue else None
        ),
        legend=False,
        size=args.size,
    )

    if args.add_means:
        # mean_values = df.groupby([args.hue])[column].mean().reset_index()

        # for i, hue in enumerate(mean_values[args.hue].values):

        # for ax in fg.axes.flat:
        ax = fg.ax
        if args.hue:
            for i, hue in enumerate(df[args.hue].unique()):
                mean = df.loc[df[args.hue] == hue][column].mean()
                ax.plot(
                    [i - 0.2, i + 0.2], [mean, mean], color="black", ls="--", lw=0.5
                )
        else:
            mean = df[column].mean()
            ax.axhline(mean, color="black", ls="--", lw=0.5)

    if args.add_hline:
        add_hline(fg.ax, args.add_hline, str(args.add_hline))
        fg.ax.set_ylim(ymax=1.0)

    if args.set_ymin:
        fg.ax.set_ylim(ymin=args.set_ymin)

    if args.set_label:
        fg.ax.set_xlabel(args.set_label)

    plt.xticks(rotation=60, ha="right", rotation_mode="anchor")
    plt.tight_layout()
    plt.savefig(args.outplot, dpi=600)
    plt.close()


def main(args):
    tab_plot_categories(args)


# EOF
