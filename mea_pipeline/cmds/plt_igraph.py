#!/usr/bin/env mea-pl

from mea_pipeline import cout, cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("Generate iGraph plot from tabular data")
    p.add_argument(
        "-s",
        "--samplefile",
        default=None,
        help="optional sample list to create vertex/node",
    )
    p.add_argument("--metafile", default=None)
    p.add_argument("--specfile", default=None)
    p.add_argument(
        "-t",
        "--threshold",
        default="0.0078125,0.015625,0.03125,0.0625,0.125,0.25,0.5,1",
    )
    p.add_argument(
        "--variable-name",
        default="IBD",
        help="variable name to be written in the title",
    )
    p.add_argument("--label-size", type=float, default=0)
    p.add_argument(
        "-o", "--outfile", default=None, help="output plot file (in PNG or PDF)"
    )
    p.add_argument(
        "--outcluster",
        default=None,
        help="output for cluster information on the last threshold",
    )
    p.add_argument("infile")
    return p


def graph_from_dataframe(df, vertices):

    import igraph as ig
    import IPython

    # generate numeric ids for each vertices
    v_keys = {}
    for idx, vert in enumerate(vertices):
        v_keys[vert] = idx

    # IPython.embed()

    # generate edges
    edges = []
    for idx, row in df.iterrows():
        idx1 = v_keys[row.iloc[0]]
        idx2 = v_keys[row.iloc[1]]
        edges.append((idx1, idx2))

    g = ig.Graph(len(vertices), edges)
    g.vs["label"] = vertices

    return g


def plt_igraph(args):

    import itertools
    import math
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt, lines
    import igraph as ig
    from mea_pipeline.utils import tabutils
    import IPython

    # read file

    infile, columns = args.infile.split(":")
    cerr(f"Reading file: {args.infile}")
    df = pd.read_table(infile, sep="\t")
    columns = columns.split(",")

    # 1st and 2nd columns are ids, 3rd are ibd fractions
    igraph_df = pd.DataFrame(
        {"id1": df[columns[0]], "id2": df[columns[1]], "weight": df[columns[2]]}
    )

    vertices = sorted(pd.concat([igraph_df.id1, igraph_df.id2]).unique())
    thresholds = [float(x) * 0.975 for x in args.threshold.split(",")]

    cerr("Generating graph(s)...")
    graphs = []
    for t in thresholds:
        curr_df = igraph_df[igraph_df.weight > t]
        g = graph_from_dataframe(curr_df, vertices=vertices)
        graphs.append(g)

    if args.outfile:

        vertex_color = None
        spec_df = None
        if args.metafile:
            # set up colours
            samples = vertices
            metafile, column = args.metafile.split(":")
            meta_df = tabutils.read_file(metafile)
            joined_df = meta_df.meta.join_to_samples(samples, ["SAMPLE", column])

            # fail guard, fill nan with empty string
            joined_df.loc[joined_df[column].isnull()] = ""

            # import IPython; IPython.embed()

            if len(joined_df) != len(samples):
                diff = set(samples) - set(joined_df["SAMPLE"])
                cerr("The following samples do not have metadata:")
                cerr(f"{diff}")
                cexit("Please recheck the sample and metadata file again")

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
                # column = column
                spec_df = tabutils.generate_spec_df(
                    sorted(joined_df[column].unique()), column
                )

            joined_df = joined_df.meta.join(spec_df, column)

            if len(joined_df) != len(samples):
                diff = set(samples) - set(joined_df["SAMPLE"])
                cerr("The following samples do not have metadata:")
                cerr(f"{diff}")
                cexit("Please recheck the sample. metadata and spec file again")

            vertex_color = joined_df["COLOUR"]

            # IPython.embed()

        # IPython.embed()
        # set layout
        cols = 3
        rows = math.ceil(len(thresholds) / cols)
        fig, axes = plt.subplots(rows, cols)
        if rows == 1:
            axes = np.array([axes])
        idx = 0

        cerr("Generating plot(s)...")
        for i, j in itertools.product(range(rows), range(cols)):
            if idx >= len(graphs):
                axes[i, j].axis("off")
                continue
            ax = axes[i, j]
            ig.plot(
                graphs[idx],
                target=ax,
                vertex_label=None if args.label_size <= 0 else False,
                vertex_label_size=args.label_size,
                vertex_size=3,
                vertex_color=vertex_color,
                vertex_frame_width=0.0675,
                edge_width=0.125,
                edge_color="#ababab",
                layout="fr",
            )
            ax.set_title(
                f"{args.variable_name} >= {thresholds[idx]*100:2.1f}%", fontsize=7
            )
            idx += 1

        # put legend on last grid/panel if spec_df is defined
        if spec_df is not None:
            ax = axes[-1, -1]
            elements = [
                lines.Line2D(
                    [],
                    [],
                    marker="o",
                    markerfacecolor=r.iloc[1],
                    label=r.iloc[0],
                    linestyle="none",
                    color="black",
                    markeredgewidth=0.25,
                )
                for i, r in spec_df.iterrows()
            ]
            # import IPython; IPython.embed()
            ax.legend(handles=elements, loc="center", fontsize="x-small")
            ax.axis("off")

        # save plots
        fig.savefig(args.outfile)
        cerr(f"Plot written to {args.outfile}")

    if args.outcluster:
        # create a cluster file from all graphs
        containers = []
        for t, curr_graph in zip(thresholds, graphs):
            sub_graphs = curr_graph.components().subgraphs()
            clusters = {}
            for idx, sg in enumerate(sub_graphs):
                clusters[idx] = sg.vs["label"]
            containers.append([t, clusters])

        import yaml

        with open(args.outcluster, mode="w") as outfile:
            yaml.dump(containers, outfile)

        cerr(f"Cluster info written to {args.outcluster}")


def main(args):
    plt_igraph(args)


# EOF
