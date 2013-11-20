nrml_converters
===============

Collection of standalone python scripts for converting OpenQuake outputs from NRML to other formats

Current Functions
=================

Of November 2013 the repository contains functions for converting the following:

1) Hazard Curves from Nrml to i) Csv, ii) Matlab Binary
2) Hazard Maps from Nrml to i) Csv, ii) .xyz (for GMT format), iii) shapefile
3) Uniform hazard spectra from Nrml to i) Csv, ii) Matlab Binary
4) Ground motion field set from Nrml to i) Csv, ii) Matlab Binary


Installation
===============

No installation is needed; however, the scripts require the following 
dependences:

* Numerical/Scientific Python (numpy/scipy)
* lxml
* Pyshp

If working in an environment where OpenQuake is already installed then the
only extra dependency is "PyShp", which can be installed using the standard
python package installer

>> sudo pip install pyshp

In other environments it is recommended to install Numpy/Scipy and lxml
using the standard packages (dependent on the OS).



