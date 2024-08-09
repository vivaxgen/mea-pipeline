# snakeutils.py
# [https://github.com/trmznt/py-snakeutils]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"
__version__ = "2024.08.09.01"

# this module provides wrapper to execute Snakemake file from Python code


import os
import sys
import pathlib
import time
import datetime
import logging
import argparse
import types
import typing


L = logging.getLogger(__name__)


__DEFAULT_RULE_PATH__ = None


def _cout(msg: str):
    print(msg, file=sys.stdout)


def _cerr(msg: str):
    print(msg, file=sys.stderr)
    sys.stderr.flush()


def _cexit(msg: str, exit_code: int = 1):
    _cerr(msg)
    sys.exit(exit_code)


def set_default_rule_path(module: types.ModuleType):
    global __DEFAULT_RULE_PATH__
    __DEFAULT_RULE_PATH__ = pathlib.Path(module.__path__[0]) / "rules"


def init_argparser(desc: str = "", p: argparse.ArgumentParser | None = None):
    """provide common arguments for snakemake-based cli"""

    p = p or argparse.ArgumentParser(description=str)
    p.arg_dict = {}

    # snakemake arguments
    p.add_argument(
        "-j", type=int, default=32, help="number of jobs to be executed in parallel"
    )
    p.add_argument(
        "--dryrun", default=False, action="store_true", help="run in dry-mode"
    )
    p.add_argument(
        "--showcmds",
        default=False,
        action="store_true",
        help="show shell commands, similar to --printshellcmds",
    )
    p.add_argument(
        "-p",
        "--printshellcmds",
        default=False,
        action="store_true",
        help="show shell commands",
    )
    p.add_argument(
        "--show-config-files",
        default=False,
        action="store_true",
        help="show the order of config files to parse",
    )
    p.add_argument(
        "--unlock",
        default=False,
        action="store_true",
        help="unlock the directory from unfinished/failed snakemake run",
    )
    p.add_argument(
        "--rerun",
        default=False,
        action="store_true",
        help="continue running the snakemake workflow from previoius point",
    )
    p.add_argument(
        "--touch",
        default=False,
        action="store_true",
        help="touch all output files, to avoid re-running the ngs-pipeline "
        "such as after modifying/debugging snakemake file",
    )
    p.add_argument(
        "--profile",
        default=None,
        help="snakemake profile to be used, also can be set from SNAKEMAKE_PROFILE env",
    )
    p.add_argument(
        "--nocluster",
        default=False,
        action="store_true",
        help="run without cluster support (eg. only on local node), useful for debugging",
    )

    # general options
    p.arg_dict["target"] = p.add_argument(
        "-t",
        "--target",
        default=[],
        action="append",
        help="target rule(s) in the snakefile [all]",
    )
    p.arg_dict["snakefile"] = p.add_argument(
        "--snakefile", default=None, help="snakemake file to be called"
    )

    # configuration files
    p.add_argument(
        "--base-config",
        default=None,
        help="path for base configuration file, relative to "
        "base environment directory",
    )
    p.add_argument(
        "--panel",
        default=None,
        help="panel to be used (eg. PANEL -> configs/PANEL.yaml as base config)",
    )
    p.add_argument(
        "-c", "--config", default=[], action="append", help="config file(s) to append"
    )
    p.add_argument(
        "-f",
        "--force",
        default=False,
        action="store_true",
        help="force the processing even if the working directory is not "
        "under current pipeline environment base directory",
    )
    p.add_argument(
        "--no-config-cascade",
        default=False,
        action="store_true",
        help="prevent from reading cascading configuration file",
    )

    return p


def check_env(env_name: str):
    if env_name in os.environ:
        return True
    return False


def setup_config(config):
    # dummy setup configurator
    return config


class SnakeExecutor(object):

    def __init__(
        self,
        # arguments parsed from the init_argparser above
        args,
        *,
        # basic configuration
        setup_config_func: typing.Callable[
            [
                dict,
            ],
            dict,
        ] = setup_config,
        # working directory
        workdir: str | pathlib.Path | None = None,
        # show configuration files in the parsing order to stderr
        show_config_files: bool = False,
        # environment base dir as root for cascading configuration
        env_basedir: str | pathlib.Path | None = None,
        # module to get the snakemake file from
        from_module: types.ModuleType | None = None,
        # default configuration file (eg from software installation)
        # if exists, this would be added as the first configuration file to parse
        default_config_file: str | pathlib.Path | None = None,
    ):

        from snakemake import cli

        self.args = args
        self.setup_config_func = setup_config_func
        self.workdir = workdir
        self.show_config_files = show_config_files or args.show_config_files
        self.env_basedir = env_basedir
        self.from_module = from_module
        self.default_config_file = default_config_file

        set_default_rule_path(self.from_module)

        # monkey-patch snakemake cli class
        cli_parse_config = cli.parse_config

        def parse_config(entries):
            if type(entries) == dict:
                return entries
            return cli_parse_config(entries)

        cli.parse_config = parse_config
        # end of monkey patching

    def run(
        self,
        # snakefile to run
        snakefile: str | pathlib.Path | None = None,
        # configuration to append/update
        config: dict = {},
        *,
        from_module: types.ModuleType | None = None,
        # allow to run the snakefile even if not inside environment directory
        force: bool = False,
        # prevent from performing cascading configuration parsing
        no_config_cascade: bool = False,
    ):

        from snakemake import cli

        cwd = self.workdir or pathlib.Path.cwd()
        if not (force or self.args.force) and not cwd.is_relative_to(self.env_basedir):
            _cexit(
                f"ERROR: current directory {cwd} is not relative to {self.env_basedir}"
            )

        snakefile = snakefile or self.args.snakefile

        # check sanity
        if not snakefile:
            _cexit(
                "ERR: Please provide snakefile to execute using --snakefile argument."
            )

        # process config files

        configfiles = [pathlib.Path(cf) for cf in reversed(self.args.config)]

        if no_config_cascade or self.args.no_config_cascade:
            if (configfile := self.env_basedir / "config.yaml").is_file():
                configfiles.append(configfile)
            config_dirs = [cwd]
        else:
            L.debug("processing cascading configuration")
            # for each config directory, check config file existence
            config_dirs = []
            config_path = cwd
            while config_path.is_relative_to(self.env_basedir):
                config_dirs.append(config_path)
                configfile = config_path / "config.yaml"
                if configfile.is_file():
                    configfiles.append(configfile)
                config_path = config_path.parent

        # get panel configuration and set as base configuration from configs/
        if self.args.panel:
            if self.args.base_config:
                _cexit(f"ERROR: cannot use both --panel and --base-config")
            self.args.base_config = "configs/" + self.args.panel + ".yaml"

        if self.args.base_config:
            configfiles.append(self.env_basedir / self.args.base_config)

        # get config file at root of environment base directory
        # this config file can be overridden by base/panel/custom config files
        if (root_configfile := self.env_basedir / "configs" / "config.yaml").is_file():
            configfiles.append(root_configfile)

        # if provided, this will be the default config file that comes with
        # the software/tool installation
        if self.default_config_file:
            configfiles.append(self.default_config_file)

        if not any(configfiles):
            _cexit(f"ERROR: cannot find any config.yaml in {config_dirs}")

        configfiles.reverse()
        if self.show_config_files:
            _cerr(f"config_files: {configfiles}")

        # setting up profile

        if "SNAKEMAKE_PROFILE" in os.environ:
            if self.args.profile is None:
                self.args.profile = os.environ["SNAKEMAKE_PROFILE"]

        if self.args.nocluster:
            # nocluster means prevent from running using batch/job scheduler,
            # then use special profile "none"
            self.args.profile = "none"

        # set targets
        if type(self.args.target) != list:
            targets = [self.args.target]
        else:
            targets = self.args.target if any(self.args.target) else ["all"]

        try:
            # set with --profile first
            if self.args.profile:
                argv = ["--profile", self.args.profile]
            else:
                argv = []

            # XXX: need to modify to use snakemake API
            L.debug("parsing snakemake arguments")
            parser, args = cli.parse_args(argv)

            args.snakefile = get_snakefile_path(
                snakefile, from_module=from_module or self.from_module
            )

            # set args further from self.args
            args.configfile = configfiles
            # cargs.config = [f'{k}={v}' for k, v in setup_config(config).items()]
            args.config = self.setup_config_func(config)
            args.targets = targets

            # running mode
            args.dryrun = self.args.dryrun
            args.touch = self.args.touch
            args.rerun_incomplete = self.args.rerun
            args.unlock = self.args.unlock

            # running parameters
            args.cores = self.args.j
            args.printshellcmds = self.args.showcmds or self.args.printshellcmds

            L.debug("invoking snakemake client")
            start_time = time.monotonic()
            status = cli.args_to_api(args, parser)
            finish_time = time.monotonic()

            return (status, datetime.timedelta(seconds=int(finish_time - start_time)))

        except Exception as e:
            cli.print_exception(e)
            sys.exit(1)

        raise RuntimeError("FATAL ERROR: should not execute this part of code")


def run_snakefile(args, config={}, workdir=None, show_configfiles=False):
    """execute snakefile based on args"""

    L.debug("importing snakemake libraries")

    import snakemake
    from snakemake import cli
    import yaml

    # monkey-patch snakemake
    cli_parse_config = cli.parse_config

    def parse_config(entries):
        if type(entries) == dict:
            return entries
        return cli_parse_config(entries)

    cli.parse_config = parse_config
    # end of monkey patching

    # getting values from environment
    L.debug("processing environment variables")
    INSTALL_BASEDIR = check_NGS_PIPELINE_BASE()
    ENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    if __envname_FORCE__ in os.environ:
        args.force = True
    if __envname_NO_CASCADE__ in os.environ:
        args.no_config_cascade = True

    # check sanity

    cwd = workdir or pathlib.Path.cwd()
    if not args.force and not cwd.is_relative_to(ENV_BASEDIR):
        _cexit(f"ERROR: current directory {cwd} is not relative to {ENV_BASEDIR}")

    configfiles = [pathlib.Path(cf) for cf in reversed(args.config)]

    if args.no_config_cascade:
        if (configfile := ENV_BASEDIR / "config.yaml").is_file():
            configfiles.append(configfile)
    else:
        L.debug("processing cascading configuration")
        # for each config directory, check config file existence
        config_dirs = []
        config_path = cwd
        while config_path.is_relative_to(SENV_BASEDIR):
            config_dirs.append(config_path)
            configfile = config_path / "config.yaml"
            if configfile.is_file():
                configfiles.append(configfile)
            config_path = config_path.parent

    # get panel configuration and set as base configuration from configs/
    if args.panel:
        if args.base_config:
            _cexit(f"ERROR: cannot use both --panel and --base-config")
        args.base_config = "configs/" + args.panel + ".yaml"

    if args.base_config:
        configfiles.append(ENV_BASEDIR / args.base_config)

    # get config file at root of NGSENV_BASEDIR
    # this config file can be overridden by base/panel/custom config files
    if (root_configfile := ENV_BASEDIR / "configs" / "config.yaml").is_file():
        configfiles.append(root_configfile)

    if not any(configfiles):
        _cexit(f"ERROR: cannot find any config.yaml in {config_dirs}")

    configfiles.reverse()

    if args.showconfigfiles:
        _cerr("Config files to read:")
        _cerr(yaml.dump([cf.as_posix() for cf in configfiles]))
        _cexit("\n")

    # for profile purposes

    if "SNAKEMAKE_PROFILE" in os.environ:
        if args.profile is None:
            args.profile = os.environ["SNAKEMAKE_PROFILE"]

    if args.nocluster:
        # nocluster means prevent from running using batch/job scheduler
        args.profile = None

    # set targets
    if type(args.target) != list:
        targets = [args.target]
    else:
        targets = args.target if any(args.target) else ["all"]

    try:
        # set with --profile first
        if args.profile:
            argv = ["--profile", args.profile]
        else:
            argv = []

        # XXX: need to modify to use snakemake API
        L.debug("parsing snakemake arguments")
        parser, cargs = cli.parse_args(argv)

        # set cargs further from args
        cargs.snakefile = get_snakefile_path(args.snakefile, from_module=ngs_pipeline)
        cargs.configfile = configfiles
        # cargs.config = [f'{k}={v}' for k, v in setup_config(config).items()]
        cargs.config = setup_config(config)
        cargs.targets = targets

        # running mode
        cargs.dryrun = args.dryrun
        cargs.touch = args.touch
        cargs.rerun_incomplete = args.rerun
        cargs.unlock = args.unlock

        # running parameters
        cargs.cores = args.j
        cargs.printshellcmds = args.showcmds

        L.debug("invoking snakemake client")
        start_time = datetime.datetime.now()
        status = cli.args_to_api(cargs, parser)
        finish_time = datetime.datetime.now()

        return (status, finish_time - start_time)

    except Exception as e:
        cli.print_exception(e)
        sys.exit(1)

    raise RuntimeError("FATAL ERROR: should not execute this part of code")


def get_snakefile_path(
    filepath: str | pathlib.Path,
    snakefile_root: pathlib.Path | None = None,
    from_module: types.ModuleType | None = None,
):
    """return real path of  snakefile"""

    if is_abs_or_rel_path(filepath):
        if type(filepath) == str:
            return pathlib.Path(filepath)
        return filepath
    if from_module is not None:
        snakefile_root = pathlib.Path(from_module.__path__[0]) / "rules"
    if snakefile_root is not None:
        return snakefile_root / filepath
    return __DEFAULT_RULE_PATH__ / filepath


def is_abs_or_rel_path(filepath: str | pathlib.Path):
    filepath = filepath.as_posix() if isinstance(filepath, pathlib.Path) else filepath
    if (
        filepath.startswith("/")
        or filepath.startswith("./")
        or filepath.startswith("../")
    ):
        return True
    return False


def scan_for_config_keywords(path):
    """return a list of keywords used as config keys in any of the
    snakemake rules
    """

    import re

    mo_bracket = re.compile(r"""config\s*\[\s*['"]([^]]*)['"]\s*\]""")
    mo_get = re.compile(r"""config\s*\.\s*get\s*\(\s*['"]([^'"]*)""")

    keywords = []

    path = pathlib.Path(path)

    for rule_file in path.glob("*.smk"):
        with open(rule_file) as f_in:
            for line in f_in:
                keys = mo_bracket.findall(line) + mo_get.findall(line)
                if any(keys):
                    keywords += keys

    return sorted(set(keywords))


# EOF
