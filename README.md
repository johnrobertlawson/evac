# **evac**

Ensemble Verification and Creation. An evolution of my WEM package. Written
and tested in Python 3.5.x on various Linux distros and Mac OS, using only
the Anaconda distributions of Python 3.

**NOT FULLY TESTED YET - still refactoring WEM.**
Documentation will eventually be here:

https://johnrobertlawson.github.io/evac/

If you're a slave to Python 2.x, it's easy to create a new Anaconda environment in
Python 3.x (replace `yourcustomname` with the name of the new environment). 

`conda create -n yourcustomname python=3.5`

It's important to use no higher than 3.5, as pygrib will not work with 3.6 via conda install.

As a sidenote, I encourage all users of Python to cite packages and acknowledge fellow
users in publications and presentations. There is often little incentive for
researchers to open source their efforts, yet it is a process that needs
to be encouraged.

## Dependencies

Thanks to weirdness with pygrib, it is recommended you install the following in this order with this command:

`conda install -c conda-forge <package>`

* pygrib 
* netCDF4
* numpy
* scipy
* matplotlib
* basemap
* xarray

## Where do I start?

Documentation incoming. For now, check out the READMEs and in-line
commenting.

## TO DO
* Folder naming - should it clash with scripts within the directories?
* Create stable master branch, release version 1.0, and fork for edits 
* Consistency with API (utc, doms, get, etc)
* Consistency with convention (datetimes, custom integers, pathlib)
* Do tutorial for running ensembles with `LazyEnsemble`, and
    verifying (stats/plot) using `Verif`.

## Contributors & Attributions

Some files or methods contain attributions to other programmers whose
code has been refactored to fit this project (or is/will become a
prerequisite). In summary, thanks to:

### SHARPpy

* Patrick Marsh
* John Hart

### HootPy/PulsatrixWx project

URL: http://www.atmos.washington.edu/~lmadaus/pyscripts.html

* David-John Gagne
* Tim Supinie
* Luke Madaus

### PyWRF project (Monash)

URL: http://code.google.com/p/pywrf/
URL: https://github.com/scaine1/pyWRF/

### PyWRFPlot project

URL: https://code.google.com/p/pywrfplot/

* Geir Arne Waagbø
