#!/usr/bin/env python
# LICENSE
#
# Copyright (c) 2013, GEM Foundation, G. Weatherill, M. Pagani, D. Monelli.
#
# The nrml_convertes is free software: you can redistribute
# it and/or modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>
#
# DISCLAIMER
#
# The software nrml_convertes provided herein is released as a prototype
# implementation on behalf of scientists and engineers working within the GEM
# Foundation (Global Earthquake Model).
#
# It is distributed for the purpose of open collaboration and in the
# hope that it will be useful to the scientific, engineering, disaster
# risk and software design communities.
#
# The software is NOT distributed as part of GEM's OpenQuake suite
# (http://www.globalquakemodel.org/openquake) and must be considered as a
# separate entity. The software provided herein is designed and implemented
# by scientific staff. It is not developed to the design standards, nor
# subject to same level of critical review by professional software
# developers, as GEM's OpenQuake software suite.
#
# Feedback and contribution to the software is welcome, and can be
# directed to the hazard scientific staff of the GEM Model Facility
# (hazard@globalquakemodel.org).
#
# The nrml_convertes is therefore distributed WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details.
#
# The GEM Foundation, and the authors of the software, assume no liability for
# use of the software.
'''
Convert NRML hazard map file to .csv file.
'''

import re
import os
import argparse
import numpy
from openquake.commonlib.nrml import read

NRML = '{http://openquake.org/xmlns/nrml/0.5}'

AK2007 = {"PGA": {"C1": 2.65, "C2": 1.39, "C3": -1.91, "C4": 4.09,
                  "logy15": 1.69, "sigma1": 1.01, "cfact": 980.665},
          "SA(0.3)": {"C1": 2.40, "C2": 1.36, "C3": -1.83, "C4": 3.56,
                      "logy15": 1.92, "sigma1": 0.88, "cfact": 980.665},
          "SA(1.0)": {"C1": 3.23, "C2": 1.18, "C3": 0.57, "C4": 2.95,
                      "logy15": 1.50, "sigma1": 0.84, "cfact": 980.665},
          "SA(2.0)": {"C1": 3.72, "C2": 1.29, "C3": 1.99, "C4": 3.00,
                      "logy15": 1.00, "sigma1": 0.86, "cfact": 980.665},
          "PGV": {"C1": 4.37, "C2": 1.32, "C3": 3.54, "C4": 3.03,
                  "logy15": 0.48, "sigma1": 0.8, "cfact": 1.}}

OPTIONAL_PATHS = [("statistics", "statistics"),
                  ("gsimlt_path", "gsimTreePath"),
                  ("smlt_path", "sourceModelTreePath"),
                  ("quantile_value", "quantileValue")]


def atkinson_kaka_2007_rsa2mmi(imt, values):
    """
    Implements the spectral acceleration to PGA conversion model of
    Atkinson & Kaka (2007)
    """
    if imt not in AK2007.keys():
        print("IMT %s not convertable to MMI" % imt)
        return values, None
    values = values * AK2007[imt]["cfact"]
    mmi = AK2007[imt]["C1"] + AK2007[imt]["C2"] * numpy.log10(values)
    idx = numpy.log10(values) > AK2007[imt]["logy15"]
    mmi[idx] = AK2007[imt]["C3"] + AK2007[imt]["C4"] * numpy.log10(values[idx])
    mmi[mmi > 10.0] = 10.0
    return mmi, AK2007[imt]["sigma1"]


def parse_nrml_hazard_map(nrml_hazard_map):
    """
    Reads the NRML file and returns the metadata as a dictionary and the value
    as a numpy array of [lon, lat, IML]
    """
    node_set = read(nrml_hazard_map, stop="node").hazardMap
    metadata = {
        "imt": node_set.attrib["IMT"],
        "investigation_time": float(node_set.attrib["investigationTime"])}
    for option, name in OPTIONAL_PATHS:
        if name in node_set.attrib:
            metadata[option] = node_set.attrib[name]
        else:
            metadata[option] = None
    if "SA" in metadata["imt"]:
        imt_str = node_set.attrib['IMT']
        m = re.search('.*\((\d*\.\d*)\).*', imt_str)
        period = m.group(1)
        metadata["sa_period"] = period
        # TODO need to fix this after the damping will be added again to
        # nrml
        # metadata['sa_damping'] = node_set.attrib['saDamping']
        metadata['sa_damping'] = '5'
    values = []
    for node in node_set.nodes:
        values.append([float(node.attrib["lon"]),
                       float(node.attrib["lat"]),
                       float(node.attrib["iml"])])
    values = numpy.array(values)
    return metadata, values


def save_hazard_map_to_csv(nrml__hazard_map_file, file_name_root,
                           to_mmi=False):
    """
    Read hazard map in `nrml__hazard_map_file` and save to .csv file
    with root name `file_name_root`
    """
    output_file = '%s.csv' % file_name_root
    if os.path.isfile(output_file):
        raise ValueError('Output file already exists.'
                         ' Please specify different name or remove old file')

    metadata, values = parse_nrml_hazard_map(nrml__hazard_map_file)
    if to_mmi:
        mmi, _ = atkinson_kaka_2007_rsa2mmi(metadata["imt"], values[:, 2])
    header = ','.join(
        ['%s=%s' % (k, v) for k, v in metadata.items() if v is not None]
    )

    header = '# ' + header + '\nlon,lat,iml'

    f = open(output_file, 'w')
    f.write(header+'\n')
    numpy.savetxt(f, values, fmt='%g', delimiter=',')
    f.close()


def save_hazard_map_to_netcdf(nrml__hazard_map_file, file_name_root,
                              to_mmi=False, resolution="10k"):
    """
    Reads the hazard map in nrml format exports to xyz file then uses GMT
    to convert to netcdf
    """
    output_file = '%s.xyz' % file_name_root
    if os.path.isfile(output_file):
        raise ValueError('Output file already exists.'
                         ' Please specify different name or remove old file')
    metadata, values = parse_nrml_hazard_map(nrml__hazard_map_file)
    values = numpy.array(values)

    if to_mmi:
        values[:, 2], _ = atkinson_kaka_2007_rsa2mmi(metadata["imt"],
                                                     values[:, 2])
    f = open(output_file, "w")
    numpy.savetxt(f, values, fmt="%g", delimiter=" ")
    f.close()

    # Convert to netcdf
    llon = numpy.min(values[:, 0])
    ulon = numpy.max(values[:, 0])
    llat = numpy.min(values[:, 1])
    ulat = numpy.max(values[:, 1])
    istring = " -I%s/%s" % (resolution, resolution)
    rstring = " -R" + "{:.5f}/".format(llon) + "{:.5f}/".format(ulon) +\
        "{:.5f}/".format(llat) + "{:.5f}".format(ulat)
    command_string = ("xyz2grd %s -G%s" % (output_file,
                      file_name_root + ".NC")) + istring + rstring
    os.system(command_string)
    os.remove(output_file)


def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert Hazard Map file from Nrml to .csv file.'
                    'To run just type: python hazard_map_converter.py '
                    '--input-file=/PATH/TO/INPUT_FILE '
                    '--output-file=/PATH/TO/OUTPUT_FILE', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
                       help='path to NRML hazard map file (Required)',
                       default=None,
                       required=True)
    flags.add_argument('--output-file',
                       help='path to output file, without file extension '
                            ' (Optional, default is root of NRML file)',
                       default=None)
    flags.add_argument('--to-mmi',
                       help="Convert the ground motion values to MMI using "
                       "the model of Atkinson & Kaka (2007) - Note this will"
                       "only apply to PGA, PGV SA(0.3), SA(1.0), SA(2.0)",
                       default=False,
                       required=False)
    flags.add_argument('--to-netcdf',
                       help="Convert the output to netcdf for use with "
                       "GMT and/or QGIS",
                       default=False,
                       required=False)
    flags.add_argument('--spacing',
                       help="Approximate spacing (km) of points: ##k",
                       default="10k",
                       required=False)

    return parser

if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        output_file = \
            os.path.splitext(parser.parse_args().input_file)[0] \
            if args.output_file is None else args.output_file
        if args.to_netcdf:
            save_hazard_map_to_netcdf(args.input_file, args.output_file,
                                      args.to_mmi, args.spacing)
        else:
            save_hazard_map_to_csv(args.input_file, output_file)
    else:
        parser.print_usage()
