#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import os
import sys
import pathlib
from mea_pipeline import arg_parser, cerr, snakeutils
import mea_pipeline


def init_argparser(desc=None):
    p = snakeutils.init_argparser(desc=desc or "run arbitrary snakefile")
    return p


def run_snakefile(args, config: dict = {}, show_status: bool = True):

    executor = snakeutils.SnakeExecutor(
        args,
        setup_config_func=setup_config,
        env_basedir=pathlib.Path.cwd(),
        from_module=mea_pipeline,
        default_config_file=pathlib.Path(os.environ["MEA_PIPELINE_BASE"])
        / "etc"
        / "default-config.yaml",
    )

    status, elapsed_time = executor.run(
        snakefile=args.snakefile,
        config=config,
        force=True,
        no_config_cascade=True,
    )

    if show_status:
        if not status:
            cerr(
                f"[WARNING: snakefile {args.snakefile} did not successfully "
                "complete]"
            )
        cerr(f"[Finish running snakefile {args.snakefile} (time: {elapsed_time})]")

    return status, elapsed_time


def setup_config(config):
    config["MEA_PIPELINE_BASE"] = os.environ["MEA_PIPELINE_BASE"]
    return config


def main(args):
    run_snakefile(args)


# EOF
