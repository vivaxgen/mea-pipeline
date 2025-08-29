#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr, cexit
import sys


def init_argparser():
    p = arg_parser("Generate manhattan and QQ plots from GWAS results")
    p.add_argument("--outqq", default="qq.png")
    p.add_argument("--outmht", default="mht.png")
    p.add_argument("--add-title-line", action="append", default=[])
    p.add_argument("--columns", default="Chr,ChrPos,PValue")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--use-fdr-bh", type=float, default=-1)
    g.add_argument("--use-bonferroni", type=float, default=-1)
    g.add_argument("--use-holm", type=float, default=-1)

    p.add_argument("infile")
    return p


def mhtplot(results_df, outplot, columns, title=None, pvalue_line=None):

    import fastlmm.util.util as flutil
    import matplotlib.pyplot as plt

    plt.rcParams["figure.figsize"] = (16.0, 8.0)
    # bonferroni_pvalue = 0.05 / len(results_df)
    flutil.manhattan_plot(
        results_df[columns].values,
        pvalue_line=pvalue_line,
        xaxis_unit_bp=False,
    )
    if title:
        plt.title(title, loc="left")
    plt.tight_layout()
    plt.savefig(outplot, dpi=600)
    plt.close()


def estimate_lambda(p_values):
    import numpy as np
    import scipy.stats as st

    LOD2 = np.median(st.chi2.isf(p_values, 1))
    return LOD2 / 0.456


def qqplot(p_values, outplot=None, title=""):
    import numpy as np
    import seaborn as sns
    from matplotlib import pyplot as plt

    p_values = p_values.copy()
    M = len(p_values)
    pnull = (0.5 + np.arange(M)) / M
    p_values[p_values > 1] = 1

    qnull = -np.log10(pnull)
    qemp = -np.log10(np.sort(p_values))

    # ax = sns.scatterplot(x=qnull, y=qemp, s=2)
    g = sns.JointGrid(x=qnull, y=qemp, marginal_ticks=True)
    g.plot_joint(sns.scatterplot, s=2, linewidth=0)
    g.plot_marginals(sns.histplot)
    ax = g.ax_joint
    max_value = max(max(qnull), max(qemp))
    ax.plot([0, max_value], [0, max_value], linewidth=0.5, c="#cecede")
    est_lambda = estimate_lambda(p_values)
    ax.text(0, max(qemp), f"est lambda = {est_lambda:7.6f}")
    ax.set_xlabel("expected -log10(P)")
    ax.set_ylabel("observed -log10(P)")

    plt.suptitle(title, x=0.1, ha="left", fontsize="medium")
    plt.tight_layout()

    if outplot:
        plt.savefig(outplot, dpi=600)
    else:
        plt.show()
    plt.close()


def plot_gwas_results(args):

    from mea_pipeline.utils import tabutils
    import pandas as pd
    from statsmodels.stats.multitest import multipletests

    df = tabutils.read_file(args.infile)

    # add title line
    title = "\n".join(args.add_title_line)
    columns = args.columns.split(",")

    qqplot(df[columns[2]], args.outqq, title=title)

    if args.use_fdr_bh > 0:

        pvals = df[columns[2]]
        _, pvals_fdr_bh, _, _ = multipletests(
            pvals, alpha=args.use_fdr_bh, method="fdr_bh"
        )
        df[columns[2]] = pvals_fdr_bh
        pvalue_line = args.use_fdr_bh

    elif args.use_bonferroni > 0:

        pvalue_line = args.use_bonferroni / len(df)

    elif args.use_holm > 0:

        pvals = df[columns[2]]
        _, pvals_holm, _, _ = multipletests(pvals, alpha=args.use_holm, method="holm")
        df[columns[2]] = pvals_holm
        pvalue_line = args.use_holm

    else:
        pvalue_line = None

    mhtplot(df, args.outmht, columns, pvalue_line=pvalue_line, title=title)


def main(args):
    plot_gwas_results(args)


# EOF
