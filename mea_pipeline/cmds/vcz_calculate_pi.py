#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr, cexit


def init_argparser():
    p = arg_parser("generate a table of population diversity from vcz data")

    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--metafile", default=None)
    p.add_argument("--winsize", type=float, default=-1, help="window size (in kb)")
    p.add_argument(
        "--whole-genome",
        default=False,
        action="store_true",
        help="calculate diversity on whole genome",
    )
    p.add_argument("--cohort", default=None)
    p.add_argument("-o", "--outfile", required=True)
    p.add_argument("infile")

    return p


def vcz_calculate_pi(args):

    import pandas as pd
    import sgkit as sg
    import xarray as xr
    from mea_pipeline.utils import tabutils

    cerr(f"[Reading genotype data from {args.infile}]")
    ds = sg.load_dataset(args.infile)

    cohort = None
    if args.metafile:
        # make a cohort based on the metafile
        # meta filename should have colon and comma-separated fields

        sample_df, diff_err = tabutils.join_metafile(ds.sample_id.values, args.metafile)
        if diff_err is not None:
            cexit("Problem exists in joining metafile to sample ids for cohort purpose")
        cohort = pd.Categorical(sample_df[args.cohort])

        ds["sample_cohort"] = xr.DataArray(cohort.codes, dims="samples")

    if args.whole_genome:
        ds = sg.window_by_genome(ds)

    ds = sg.divergence(ds, merge=True)

    import IPython

    IPython.embed()

    # generate matrix per cohort
    keys = ["all"]
    if cohort is not None:
        keys = cohort.categories

    dfs = []
    contig = ds.contig_id.values[ds.variant_contig]
    for idx, key in enumerate(keys):
        df = pd.DataFrame(
            {
                "CONTIG": contig,
                "POS": ds.variant_position.values,
                "PI": ds.stat_divergence[:, idx, idx],
                args.cohort: key,
            }
        )
        dfs.append(df)

    out_df = pd.concat(dfs)
    if args.outfile:
        tabutils.write_file(args.outfile, out_df)


def main(args):
    vcz_calculate_pi(args)


# EOF
