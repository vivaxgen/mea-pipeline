mea-pl
======

Synopsis
--------

**mea-pl** [-h] [-i] [-l]

**mea-pl** <*COMMDAN*> [*OPTIONS] [--debug]

Description
-----------

:program:`mea-pl` is the front end to execute all MEA-Pipeliene commands.

Options
-------

.. program:: mea-pl

.. options::  -h, --help

   Show this help message and exit.
    
.. options:: -i

   Run in interactive mode, using IPython shell.
    
.. options:: -l, --list

   List all available commands and exit.
    
.. option:: --debug

   Enable debug mode for detailed logging.

Environment
-----------

.. envvar:: MEA_PIPLINE_BASE

   Path to the MEA-Pipeline installation directory.

.. envvar:: MEA_PIPLINE_LOGLEVEL

   Set the logging level, with numeric values.

.. envvar:: MEA_PIPELINE_CMD_MODS

   Colon-separated list of additional directories to search for command
   modules.
