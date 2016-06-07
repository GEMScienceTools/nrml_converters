nrml_converters
===============

Collection of standalone python scripts for converting NRML files (as generated
or consumed by OpenQuake) to other formats.

Current Functions
=================

The ``oq_input`` folder contains a script for converting OpenQuake NRML input
source model and rupture model files to ESRI shapefile, and vice versa.

The ``oq_output`` folder contains scripts for converting OpenQuake NRML output
files (hazard curves, hazard maps, uniform hazard spectra, stochastic event
sets, ground motion fields, disaggregation matrices) to .csv (comma separated)
or .txt (tab delimited) files.

Each script can be executed from shell by invoking the python command followed
by the script name, and by one or more flag arguments (depending on the script).

For each script, an ``help`` flag is available providing instructions for use
(just type: python SCRIPT_NAME.py --help).

Installation
===============

No installation is needed; however, the scripts require the following 
dependences:

* numpy
* oq-hazardlib (https://github.com/gem/oq-hazardlib) - only
    for eventset_converter.py
* oq-nrmllib (https://github.com/gem/oq-nrmllib) - for source_model_converter.py,
    rupture_model_converter.py and hazard_curve_converter.py
* shapely - only for source_model_converter.py
* pyshp - only for source_model_converter.py
* GMT (http://gmt.soest.hawaii.edu) - only for disaggregation_converter.py

If working in an environment where OpenQuake is already installed then the first
five python dependencies are already available. The only missing one is ``pyshp``
which can be installed using the following command:

>> sudo pip install pyshp

In other environments it is recommended to install numpy and shapely
using the standard packages (dependent on the OS).



