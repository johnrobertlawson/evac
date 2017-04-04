evac package
===========

Tutorial
--------

This will lead you through an example for automating WRF runs, 
creating statistics from an ensemble, and plotting data!

Installation
------------

Ensure you have `git` installed on your system or server. Then execute ``git
clone https://github.com/johnrobertlawson/evac.git``. The example scripts are
located in ``evac/evac/examples/``. You can copy a `.py` file from there into your
own personal scripts folder. `evac` works best when you don't interact directly
with the codebase, but only change the top-level script.

The following allows you to use evac without installing anything:

.. code-block:: python

    import sys
    sys.append('path/to/evac')

Make sure you change the path to where you have downloaded the `evac` codebase
from GitHub. Next:

.. code-block:: python

    from evac.datafiles.ensemble import Ensemble


Examples
--------

This creates an instance of the Ensemble class. An example of how to load a directory
of WRF out files is as follows:

.. code-block:: python

    outdir = '/absolute/path/to/figures/'
    ncdir = '/absolute/path/to/data'

More to follow...

