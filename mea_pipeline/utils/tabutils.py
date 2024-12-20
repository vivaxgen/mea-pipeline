# tabutils.py

# this file is derived from https://github.com/trmznt/seqpy

import pandas as pd
import numpy as np
import pathlib
from enum import Enum
from typing import Literal
from mea_pipeline import cerr


# tabular format for genomic data
#
# Sample-based Rows
# =================
#
# headers: SAMPLE, CHR:POS, ...
#
# each row contains a single sample
#
# Variant-based Column
#


def read_file(infile):
    match exts := pathlib.Path(infile).suffixes:
        case [*_, ".feather"]:
            return pd.read_feather(infile)
        case [*_, ".txt"] | [*_, ".tsv"] | [*_, ".tsv", ".gz"]:
            return pd.read_table(infile, sep="\t")
        case [*_, ".csv"] | [*_, ".csv", ".gz"]:
            return pd.read_table(infile, sep=",")

    raise RuntimeError(f"Unknown extension type: {exts}")


def write_file(outfile, dataframe):
    match pathlib.Path(outfile).suffixes:
        case [*_, ".feather"]:
            return dataframe.to_feather(outfile)
        case [*_, ".tsv"] | [*_, ".txt"] | [*_, ".tsv", ".gz"] | [*_, ".txt", ".gz"]:
            return dataframe.to_csv(outfile, sep="\t", index=False)
        case [*_, ".csv"] | [*_, ".csv", ".gz"]:
            return dataframe.to_csv(outfile, sep=",", index=False)
    raise RuntimeError("Unknown extension type")


def read_table(infile, func, **kwargs):
    """func is function that returns the line or None"""
    import io

    lines = []
    with open(infile) as f_in:
        for line in f_in:
            if (line := func(line)) is not None:
                lines.append(line)
    buf = "\n".join(lines)
    return pd.read_table(io.StringIO(buf), **kwargs)


def join_metafile(samples, metafile, data=None, start_col=0, percenttag=False):
    tokens = metafile.split(":")
    meta_df = read_file(tokens[0])
    if len(tokens) == 1:
        # we use all columns, and set the 1st columns as key
        columns = meta_df.columns
    else:
        columns = tokens[1].split(",")

        # if only 1 column name, assumed that it will be the indexing (ie sample code) column
        # and the all of the other columns will be joined
        if len(columns) == 1:
            columns = [columns[0]] + meta_df.columns

    joined_df = meta_df.meta.join_to_samples(samples, columns, percenttag)

    diff = None
    if len(joined_df) != len(samples):
        diff = set(samples) - set(joined_df["SAMPLE"])
        cerr("The following samples do not have metadata:")
        cerr(f"{diff}")

    if data is not None:
        # perform sanity check
        if len(joined_df) != len(data):
            raise ValueError("data df does not have same length as sample size")

        # concatenate joined_df with data_df starting from column start_col
        joined_df = pd.concat([joined_df, data.iloc[:, start_col:]], axis=1)

    return joined_df, diff


@pd.api.extensions.register_dataframe_accessor("meta")
class MetaAccessor(object):

    def __init__(self, df):
        self._df = df

    def join_to_samples(
        self,
        samples,
        columns,
        percenttag=False,
        how: Literal["left", "right", "outer", "inner"] = "inner",
    ):

        actual_columns = []
        final_columns = []
        rename_dict = {}
        for c in columns:

            # if c contains '>', then the column will be renamed
            c, *c_rename = c.split("=", 1)

            # check if c exists
            if c not in self._df:
                raise ValueError(f"Column {c} does not exits")

            actual_columns.append(c)
            if any(c_rename):
                rename_dict[c] = c_rename[0]
                final_columns.append(c_rename[0])
            else:
                final_columns.append(c)

        # use pandas' join mechanism
        sample_df = pd.DataFrame(dict(SAMPLE=samples))
        joined_df = sample_df.join(
            self._df.loc[:, actual_columns].set_index(actual_columns[0]),
            on="SAMPLE",
            how=how,
        )

        # rename columns if needed
        if any(rename_dict):
            joined_df.rename(columns=rename_dict, inplace=True)

        # prepend percent signature if needed
        if percenttag:
            names = {}
            for c in final_columns[1:]:
                names[c] = f"%{c}"
            joined_df.rename(columns=names, inplace=True)

        return joined_df

    def join(self, right_df, on_column):

        # use pandas' join mechanism
        # import IPython; IPython.embed()
        return self._df.join(right_df.set_index(on_column), on=on_column, how="left")

    def to_dict(self, key_column, value_column):
        return dict(zip(self._df[key_column], self._df[value_column]))


class FormatType(Enum):
    SAMPLE = 1
    VARIANT = 2


@pd.api.extensions.register_dataframe_accessor("geno")
class GenoAccessor(object):

    def __init__(self, df):
        self._type = None
        self._df = df
        self._sampledf: pd.DataFrame | None = None
        self._metadf: pd.DataFrame | None = None
        self._validate_process(df)

    def _validate_process(self, df):
        if df.columns[0] == "SAMPLE":
            # sample-based rows
            self._type = FormatType.SAMPLE
            self._process_sample_format()

        elif df.columns[0] == "CHROM" and df.columns[1] == "POS":
            # position-based
            self._type = FormatType.VARIANT
            self._process_variant_format()

        else:
            raise AttributeError("This dataframe does not conform geno specs")

        return True

    def to_sampletype(self):

        # check that we start from position
        if self._type != "variant":
            raise ValueError("Dataframe is already position-based")

    def to_positiontype(self):

        # check that we start from sample
        if self._type != "sample":
            raise ValueError("Dataframe is already sample-based")

    def get_samplelist(self):
        if self._type == "sample":
            return self._df["SAMPLE"]
        elif self._type == "position":
            columns = self._df.columns[2:]
            sample_indexes = ~columns.str.starswith("M|")
            return columns[sample_indexes]

    def get_positionlist(self):
        """return a list of tuple (chrom, pos)"""
        positions = self.get_alleles()
        positionlist = []
        for position in positions:
            chrom, pos = position.split(":", 1)
            positionlist.append((chrom, int(pos)))
        return positionlist

    def get_positionids(self):
        pass

    def get_alleles(self):
        """return a dataframe with columns of [MARKER1, MARKER2, MARKER3, ...]"""
        return self.get_alleledf()

    def get_alleledf(self):
        return self._alleledf

    def get_samples(self):
        """return series of sample codes"""
        return self._sampledf["SAMPLE"] if self._sampledf else None

    def get_sampledf(self):
        """return dataframe containings [SAMPLE, %META1, [...], %MARKER1, ...]"""
        return self._sampledf

    def get_metarows(self):
        return self._metarows

    def get_metacolumns(self):
        return self._metacolumns

    def set_missing_alleles(self, missings=None):
        if missings:
            for allele in missings:
                self._alleledf.replace(allele, np.nan, inplace=True)
        return self._df

    def _process_sample_format(self):
        # create ** new dataframe view ** without metadata row & columns,

        metarows = self._df["SAMPLE"].str.startswith("%")
        # metarows in pandas Series

        allelecolumns = self._df.columns.str.contains(":")
        # allelecolumns is numpy array

        self._idx_samplerows = metarows.index[~metarows]
        self._idx_metarows = metarows.index[metarows]

        self._idx_allelecolumns = np.where(allelecolumns)[0]
        self._idx_metacolumns = np.where(~allelecolumns)[0]

        # self._sampledf = self._df.loc[~metarows]
        # self._metarows = self._df.loc[metarows]
        # self._metacolumns = self._df.loc[~metarows, ~allelecolumns]
        # self._alleledf = self._df.loc[~metarows, allelecolumns]

        self._sampledf = self._df.iloc[self._idx_samplerows]
        self._metarows = self._df.iloc[self._idx_metarows]
        self._metacolumns = self._df.iloc[self._idx_samplerows, self._idx_metacolumns]
        self._alleledf = self._df.iloc[self._idx_samplerows, self._idx_allelecolumns]

    def _process_variant_format(self):
        # create new dataframe without metadata columns
        raise NotImplementedError("the functionality is still being worked on")
        metacolumns = self._df.columns.str.startswith("%")
        self._metacolumns = self._df.loc[:, metacolumns]
        samplecolumns = self._df.columns[~metacolumns]

    def to_mhaps(self, mhap_specs):
        """return a new DF containing microhaplotype strings"""
        pass


@pd.api.extensions.register_dataframe_accessor("mhap")
class MicroHaplotypeAccessor(object):

    def __init__(self, df):
        self._df = df

    def to_numerics(self):
        """return a new DF containing microhaplotype numerics"""
        # for each mhap:
        #   perform value_count(), sort and convert to dict
        #   create new series based on the dict
        pass


def dataframe_from_variants(dataset, variants):

    from mea_pipeline.utils import sgio

    position_ids = sgio.get_position_ids(dataset)

    # import IPython; IPython.embed()

    cerr("[Creating new dataframe for alleles...]")
    df = pd.DataFrame(data=variants, columns=position_ids)
    df.insert(0, "SAMPLE", dataset.sample_id)

    return df


def generate_spec_df(a_list, column_name):
    """this function generate a spec df (currently color only, symbol/marker for future)"""

    import itertools

    # 12 colours from ColorBrewer2
    colour12 = [
        "#1f78b4",
        "#33a02c",
        "#e31a1c",
        "#ff7f00",
        "#6a3d9a",
        "#b15928",
        "#a6cee3",
        "#b2df8a",
        "#fb9a99",
        "#fdbf6f",
        "#cab2d6",
        "#ffff99",
    ]

    # 13 colours from https://gist.github.com/xgfs/37436865b6616eebd09146007fea6c09
    xgfs_normal13 = [
        "#ebac23",
        "#b80058",
        "#008cf9",
        "#006e00",
        "#00bbad",
        "#d163e6",
        "#b24502",
        "#ff9287",
        "#5954d6",
        "#00c6f8",
        "#878500",
        "#00a76c",
        "#bdbdbd",
    ]

    # 20 colours from Vega
    colour20 = [
        "#1f77b4",
        "#aec7e8",
        "#ff7f0e",
        "#ffbb78",
        "#2ca02c",
        "#98df8a",
        "#d62728",
        "#ff9896",
        "#9467bd",
        "#c5b0d5",
        "#8c564b",
        "#c49c94",
        "#e377c2",
        "#f7b6d2",
        "#7f7f7f",
        "#c7c7c7",
        "#bcbd22",
        "#dbdb8d",
        "#17becf",
        "#9edae5",
    ]

    # 20 colours from IWantHue
    colour20 = [
        "#7ebbc7",
        "#6d38c0",
        "#72dc59",
        "#cc4dc4",
        "#cfd94a",
        "#6f68c9",
        "#649a40",
        "#c2477b",
        "#68d5a9",
        "#cd3f41",
        "#637dac",
        "#dc6931",
        "#4d285e",
        "#bf953c",
        "#cc9acd",
        "#536840",
        "#74372c",
        "#c9d19d",
        "#363638",
        "#c69085",
    ]

    colour31 = """
173, 216, 230
  0, 191, 255
 30, 144, 255
  0,   0, 255
  0,   0, 139
 72,  61, 139
123, 104, 238
138,  43, 226
128,   0, 128
218, 112, 214
255,   0, 255
255,  20, 147
176,  48,  96
220,  20,  60
240, 128, 128
255,  69,   0
255, 165,   0
244, 164,  96
240, 230, 140
128, 128,   0
139,  69,  19
255, 255,   0
154, 205,  50
124, 252,   0
144, 238, 144
143, 188, 143
 34, 139,  34
  0, 255, 127
  0, 255, 255
  0, 139, 139
128, 128, 128
"""

    # if grouping <= 12, use colour12 otherwise use colour20
    colour_cycle = itertools.cycle(xgfs_normal13 if len(a_list) <= 13 else colour20)
    if len(a_list) > 20:
        import colorcet as cc

        # another option is https://github.com/taketwo/glasbey or
        # manually from https://mokole.com/palette.html
        colour_cycle = itertools.cycle(cc.b_glasbey_bw)
    colour_set = [next(colour_cycle) for c in range(len(a_list))]

    spec_df = pd.DataFrame({column_name: a_list, "COLOUR": colour_set})

    return spec_df


# EOF
