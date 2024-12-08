#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr, cexit, gzopen
from mea_pipeline.utils import tabutils
import sys


def init_argparser():
    p = arg_parser("Create colour annotation based on header of a file")
    p.add_argument(
        "--metafile",
        required=True,
        help="TSV/CSV file containing SAMPLES and group columns, "
        "with SAMPLES and group columns stated after filename, eg: "
        "a_filename.tsv:SAMPLE,Country",
    )
    p.add_argument(
        "--specfile",
        default=None,
        help="TSV/CSV file containing the group column and other columns",
    )
    p.add_argument(
        "-s",
        default=False,
        action="store_true",
        help="include sample size of each group in legend",
    )
    p.add_argument("-o", "--outprefix", default="")
    p.add_argument("--outsample", default="")
    p.add_argument("--outgroup", default="")

    p.add_argument(
        "infile",
        help="either a tab-delimited file with a SAMPLE header, or "
        "a tab-delimited file with samples as its header, such as "
        "a distance matrix file.",
    )

    return p


def tab2anno(args):

    if not (args.outprefix or args.outsample or args.outgroup):
        cexit("Please provide at least --outprefix or --outsample or --outgroup")

    if args.outprefix:
        args.outsample = args.outprefix + ".indv.tsv"
        args.outgroup = args.outprefix + ".group.tsv"

    # read infile using tabutils
    df = tabutils.read_file(args.infile)

    # use geno extension, if error then just read headers
    try:
        samples = df.geno.get_samples()
    except AttributeError:
        samples = df.columns
    cerr(f"Reading {len(samples)} samples from {args.infile}")

    # if read the metadata file
    # get filename and specs
    metafile, columns = args.metafile.split(":")
    meta_df = tabutils.read_file(metafile)
    if not columns:
        columns = meta_df.columns[[0, 1]]
    else:
        columns = columns.split(",")

    joined_df = meta_df.meta.join_to_samples(samples, columns)

    if len(joined_df) != len(samples):
        diff = set(samples) - set(joined_df["SAMPLE"])
        cerr("The following samples do not have metadata:")
        cerr(f"{diff}")
        sys.exit(1)

    if args.specfile:
        # get spec file and its specs
        specfile, column = args.specfile.split(":")
        if not column:
            raise ValueError(
                f"Please provide the column key at the end of filename, eg: "
                f"{specfile}:COLUMN_NAME"
            )
        spec_df = tabutils.read_file(specfile)
    else:
        # generate spec_df automatically
        column = columns[1]
        spec_df = tabutils.generate_spec_df(sorted(joined_df[column].unique()), column)

    joined_df = joined_df.meta.join(spec_df, column)

    if len(joined_df) != len(samples):
        diff = set(samples) - set(joined_df["SAMPLE"])
        cerr("The following samples do not have metadata:")
        cerr(f"{diff}")
        sys.exit(1)

    if args.outsample:
        tabutils.write_file(args.outsample, joined_df)
        cerr(f"Sample output written to {args.outsample}")

    if args.outgroup:
        # need to change 1st column label

        if args.s:
            sizes = joined_df[column].value_counts()

            for idx, r in spec_df.iterrows():
                spec = r[column]
                r[column] = f"{spec} ({sizes[spec]})"

        spec_df.sort_values(by=column, inplace=True)
        spec_df.rename(columns={column: "GROUP"}, inplace=True)

        tabutils.write_file(args.outgroup, spec_df)
        cerr(f"Group output written to {args.outgroup}")


def main(args):
    tab2anno(args)


# EOF
