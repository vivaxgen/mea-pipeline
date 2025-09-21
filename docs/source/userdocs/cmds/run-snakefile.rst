mea-pl run-snakefile
====================

Synopsis
--------

**mea-pl** run-snakefile [*OPTIONS*] [--debug]

Description
-----------

:program:`mea-pl run-snakefile` runs a snakefile.

Options
-------

.. program:: mea-pl run-snakefile

.. option::  --snakefile SNAKEFILE

   Path to the snakefile to be executed. (required)
   If path can be absolute or relative to the current working directory,
   (eg. ./my_snakefile.smk).
   If only provided the name of the snakefile (eg. my_snakefile.smk), it will
   be searched from the rules/ directory of the MEA-Pipeline installation.

.. option:: -c, --configfile CONFIGFILE

   Path to a configuration file in YAML or JSON format.
   If path can be absolute or relative to the current working directory,
   (eg. ./my_config.yaml).
   If only provided the name of the config file (eg. my_config.yaml), it will
   be searched from the config/ directory of the MEA-Pipeline installation.
   Multiple config files can be provided by repeating this option.

.. option:: -t, --targets TARGETS

   If not provided, the default target specified in the snakefile will be built.
   Multiple targets can be provided by repeating this option.

Environment
-----------

.. envvar:: SNAKEFILE_PROFILE

   Directory for Snakemake profile