# `plot`

Create figures and wrap matplotlib commands.

* `Clicker()` creates interactive and pop-up plots.
* `colourtables.py` has look up tables
* `scales.py` sets default limits and contouring for plots

| Class | Inherits | Description |
| --- | --- | --- |
| `Figure` | | Parent class for figures |
| | | |
| `BirdsEye()` | `Figure()` | Top-down plots |
| `Clicker()` | `Figure()` | Interactive plots |
| `HeatMap()` | `Figure()` | Heatmap/matrix plots |
| `Map()` | `Figure()` | Generic maps, like domains and labelled figures |
| `SkewT()` | `Figure()` | Skew-T sounding plots |
| `TimeSeries()` | `Figure()` | Time series plots |
| `XSection()` | `Figure()` | Cross-section plots |

## TO DOs

* Merge `map.py` and `maps.py`
* Merge `lookuptable.py` and `colourtables.py`
* Move `rucplot.py` into datafiles and other files here
* Maybe change `scales.py` from class to methods.
* Merge `skewt.py` and `skewt2.py` scripts
