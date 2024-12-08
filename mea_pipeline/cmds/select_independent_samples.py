#!/usr/bin/env spcli

from mea_pipeline import cout, cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser(
        "Get independent samples from the least missingness samples from each cluster"
    )
    p.add_argument("-o", "--outfile")
    p.add_argument("--sample-qc", required=True)
    p.add_argument("infile")
    return p


def select_independent_samples(args):

    import pandas as pd
    import yaml

    sample_df = pd.read_table(args.sample_qc)
    containers = yaml.load(open(args.infile), Loader=yaml.SafeLoader)

    clusters = containers[-1][1]
    samples = []
    for cluster in clusters.values():
        if len(cluster) == 1:
            samples.append(cluster[0])
            continue
        cluster_sample_df = sample_df.loc[sample_df.SAMPLE.isin(cluster), :]
        sorted_sample_df = cluster_sample_df.sort_values(by="F_MISSING", ascending=True)
        samples.append(sorted_sample_df.SAMPLE.iloc[0])

    samples = sorted(samples)
    with open(args.outfile, "w") as out_f:
        out_f.write("\n".join(samples))


def main(args):
    select_independent_samples(args)


# EOF
