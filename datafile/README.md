# datafiles

## Helper classes

* The `Ensemble()` class creates a dictionary of ensemble data file members.
* The `Obs()` class collates observational data of a given type so that e.g. precipitation data can be inspected at one across numerous times.
* ?

## Data file formats

### Family tree:

| Class name | Inherits | Description |
| --- | --- | --- |
| `DataFile()`  | | Top level class  |
||||
| `GribFile()` | `DataFile()` | Grib2 format |
| `NCFile()` | `DataFile()` | NetCDF4 format |
| `PNGFile()` | `DataFile()` | PNG format |
| `CSVFile()` | `DataFile()` | CSV format |
| `BinaryFile()` | `DataFile()` | Any old binary file|
||||
| `ECMWFFile()`| `NCFile()` | ECMWF data (NetCDF)|
| `WRFOut()` | `NCFile()` | Output from WRF-ARW (NetCDF) |
| `HRRR()` | `GribFile()` | HRRR model (Grib)|
| `StageIV()` | `GribFile()` | STAGEIV accumulated precip (Grib) |
| `Radar()` | `PNGFile()` | Radar data (PNG) |
| `SPCLSR()` | `CSVFile()` | SPC storm reports (CSV) |
| `StormReports()` | `CSVFile()` | Storm reports from where? |
| `MRMS()` | `BinaryFile()` | MRMS rainfall etc data (binary)
||||
| `RUC()` | `WRFOut()` | RUC or RAP model (?)|

## TO DO

* `Radar()` should allow different inheritance depending on its format.
* `DataFile()` could wrap plotting scripts for data.
