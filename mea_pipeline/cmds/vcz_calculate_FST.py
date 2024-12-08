#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr, cexit


def init_argparser():
    p = arg_parser("calculate FST between 2 cohorts")
    p.add_argument("--metafile")
    p.add_argument("-o", "--outfile")
    p.add_argument("--cohort-notation", default=None)
    p.add_argument("infiles", nargs="+")

    return p


def vcz_calculate_FST(args):

    from mea_pipeline.utils import sgio
    from mea_pipeline.utils import tabutils
    import pandas as pd
    import sgkit as sg
    import xarray as xr

    match len(args.infiles):
        case 1:
            if args.cohort_notation is None:
                cexit("--cohort-notation is rquired for single input file")

            cohorts = args.cohort_notation.split(":")
            if len(cohorts) != 2:
                cexit("Invalid cohort-notation format")
            cohorts.sort()

            ds = sgio.load_dataset(args.infiles[0])

            samples = ds.sample_id.values
            cerr("INFO: joining metafile to sample ids for cohort purpose")
            sample_df, diff_err = tabutils.join_metafile(samples, args.metafile)
            if diff_err is not None:
                cexit(
                    "Problem exists in joining metafile to sample ids for cohort purpose"
                )

            categories = pd.Categorical(sample_df.iloc[:, 1], categories=cohorts)
            ds["sample_cohort"] = xr.DataArray(categories.codes, dims="samples")

        case 2:
            infile1, infile2 = args.infiles
            ds1 = sgio.load_dataset(infile1)
            ds2 = sgio.load_dataset(infile2)

            # set cohort of 0 to ds1
            ds1["sample_cohort"] = 0
            # set cohort of 1 to ds2
            ds2["sample_cohort"] = 1

            # merge ds1 and ds2
            ds = xr.concat([ds1, ds2], dim="samples")

    cerr("INFO: calculating FST")
    ds = sg.Fst(ds)
    fst = ds.stat_Fst[:, 0, 1].values

    cerr("INFO: calculating divergence")
    ds = sg.divergence(ds)
    div = ds.stat_divergence[:, 0, 1].values
    pi0 = ds.stat_divergence[:, 0, 0].values
    pi1 = ds.stat_divergence[:, 1, 1].values

    if args.outfile is not None:
        df = pd.DataFrame(
            {
                "CONTIG": ds.contig_id[ds.variant_contig.values].values,
                "POS": ds.variant_position.values,
                "DIV": div,
                "PI0": pi0,
                "PI1": pi1,
                "FST": fst,
            }
        )
        df.to_csv(args.outfile, index=False)

    import IPython

    IPython.embed()


def main(args):
    vcz_calculate_FST(args)


# EOF
