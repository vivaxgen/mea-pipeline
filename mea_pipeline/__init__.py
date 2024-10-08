# __init__.py
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from . import subcommands


cerr = subcommands._cerr
cout = subcommands._cout
cexit = subcommands._cexit

arg_parser = subcommands.arg_parser


def setup_config(d={}):
    d["MEA_PIPELINE_RULEDIR"] = __path__[0] + "/rules"
    return d


def gzopen(filename, options="rt"):

    if filename.endswith(".gz"):
        import gzip

        return gzip.open(filename, options)
    else:
        return open(filename, options)


# EOF
