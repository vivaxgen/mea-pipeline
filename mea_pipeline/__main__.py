#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

# [https://github.com/vivaxgen/mea-pipeline]
__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
import os
import logging

# check that we have MEA_PIPELINE_BASE environemt
if "MEA_PIPELINE_BASE" not in os.environ:
    print(
        "ERROR: please set the environment by executing MEA-Pipeline activation script",
        file=sys.stderr,
    )
    sys.exit(1)

# prepare logging ahead of everything else
if (LOGLEVEL := int(os.environ.get("MEA_PIPELINE_LOGLEVEL", 0))) > 0:
    if LOGFILE := os.environ.get("MEA_PIPELINE_LOGFILE", ""):
        logging.basicConfig(filename=LOGFILE, level=LOGLEVEL)
    else:
        logging.basicConfig(level=LOGLEVEL)


from mea_pipeline.subcommands import _cerr as cerr, _cexit as cexit, SubCommands
import platform


def greet():
    cerr(
        f'{sys.argv[0].split("/")[-1]} - MEA-Pipeline command line interface\n'
        f"[https://github.com/vivaxgen/mea-pipeline]"
    )
    cerr(f"Host: {platform.uname().node}")


def usage():
    cexit("  usage:\n      mea-pl CMD [ARGS]\n  try: mea-pl showcmds")


def main():

    cmds = SubCommands(
        modules=["mea_pipeline.cmds"],
        module_env="MEA_PIPELINE_CMD_MODS",
        env_takes_precedence=True,
        allow_any_script=True,
        allow_shell=True,
        greet_func=greet,
        usage_func=usage,
    )

    cmds.main()


if __name__ == "__main__":
    main()

# EOF
