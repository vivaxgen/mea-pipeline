#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser()
    p.add_argument("-o", "--outfile", required=True)
    p.add_argument("--outlog", default=None)
    p.add_argument("--metafile", default=None)
    p.add_argument("-n", "--components", type=int, default=3)
    p.add_argument("infile")
    return p


def dis_calculate_pcoa(args):

    from mea_pipeline.utils.tabutils import read_file, write_file, join_metafile
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import pandas as pd

    cerr(f"Reading distance matrix from {args.infile}")

    distm_df = read_file(args.infile)
    samples = distm_df.columns

    # prepare metadata
    if args.metafile:
        sample_df, errs = join_metafile(samples, args.metafile)
    else:
        sample_df = pd.DataFrame(dict(SAMPLE=samples))

    # calculate PCA
    scaled_distm = StandardScaler().fit_transform(distm_df.values)
    pca = PCA(n_components=args.components)
    pca_feats = pca.fit_transform(scaled_distm)
    pca_df = pd.DataFrame(
        data=pca_feats, columns=[f"PC{x + 1}" for x in range(args.components)]
    )

    # join sample_df and pca_df
    df = pd.concat([sample_df, pca_df], axis=1)
    write_file(args.outfile, df)
    cerr(f"PCoA written to {args.outfile}")

    # report on variance
    if args.outlog:
        with open(args.outlog, "w") as f_out:
            for i in range(len(pca.explained_variance_ratio_)):
                f_out.write(f"  PC{i + 1}: {pca.explained_variance_ratio_[i]}\n")


def main(args):
    dis_calculate_pcoa(args)


# EOF
