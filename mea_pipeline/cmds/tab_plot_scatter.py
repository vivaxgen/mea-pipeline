#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from mea_pipeline import arg_parser, cerr, cexit


def init_argparser():
    p = arg_parser()
    p.add_argument("-o", "--outplot", required=True)
    p.add_argument("--axis", default=[], action="append")
    p.add_argument("--hue", default=None)
    p.add_argument("--dpi", default=600, type=int)
    p.add_argument("infile")
    return p


def tab_plot_scatter(args):

    import pandas as pd
    import seaborn as sns
    from matplotlib import pyplot as plt

    from mea_pipeline.utils.tabutils import read_file
    from mea_pipeline.color_palettes import color_palettes

    # read file
    df = read_file(args.infile)

    if args.hue:
        df.sort_values(by=args.hue, inplace=True)

    # check the axis
    n_plots = len(args.axis)
    size = 5

    fig, axes = plt.subplots(1, n_plots, figsize=(size * n_plots, size))

    for ax, notation in zip(axes, args.axis):
        columns = notation.split(":")
        if len(columns) != 2:
            cexit(f"ERR: axis {notation} is not in format column_x:column_y")

        sns.scatterplot(
            x=columns[0],
            y=columns[1],
            data=df,
            ax=ax,
            hue=args.hue,
            size=5,
            palette=sns.color_palette(color_palettes["xgfs_normal12"]),
            legend=False,
        )

    fig.tight_layout()
    fig.savefig(args.outplot, dpi=args.dpi)


def main(args):
    tab_plot_scatter(args)


# EOF
