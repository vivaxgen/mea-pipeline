mea-pl
======

Synopsis
--------

**mea-pl** [-h] [-i] [-l]

**mea-pl** *COMMAND* [*OPTIONS*] [--debug]

Description
-----------

:program:`mea-pl` is the front end to execute all MEA-Pipeline commands.

Options
-------

.. program:: mea-pl

.. option::  -h, --help

   Show this help message and exit.
    
.. option:: -i

   Run in interactive mode, using IPython shell.
    
.. option:: -l, --list

   List all available commands and exit.
    
.. option:: --debug

   Enable debug mode for detailed logging.

Environment
-----------

.. envvar:: MEA_PIPELINE_BASE

   Path to the MEA-Pipeline installation directory.

.. envvar:: MEA_PIPELINE_LOGLEVEL

   Set the logging level, with numeric values.

.. envvar:: MEA_PIPELINE_CMD_MODS

   Colon-separated list of additional directories to search for command
   modules.
