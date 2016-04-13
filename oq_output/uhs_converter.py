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
Convert NRML uniform hazard spectra file to .csv file.
'''

import os
import argparse
import numpy
#from lxml import etree
from openquake.commonlib.nrml import read_lazy
import matplotlib.pyplot as plt

NRML='{http://openquake.org/xmlns/nrml/0.5}'
GML='{http://www.opengis.net/gml}'


#def parse_nrml_uhs_curves(nrml_uhs_map):
#    """
#    Parse NRML uhs file.
#    """
#    metadata = {}
#    periods = None
#    values = []
#
#    parse_args = dict(source=nrml_uhs_map)
#    for _, element in etree.iterparse(**parse_args):
#        if element.tag == '%suniformHazardSpectra' % NRML:
#            a = element.attrib
#            metadata['statistics'] = a.get('statistics')
#            metadata['quantile_value'] = a.get('quantileValue')
#            metadata['smlt_path'] = a.get('sourceModelTreePath')
#            metadata['gsimlt_path'] = a.get('gsimTreePath')
#            metadata['investigation_time'] = a['investigationTime']
#            metadata['poe'] = a.get('poE')
#        elif element.tag == '%speriods' % NRML:
#            periods = map(float, element.text.split())
#        elif element.tag == '%suhs' % NRML:
#            lon, lat = map(
#                float, element.find('%sPoint/%spos' % (GML, GML)).text.split()
#            )
#            imls = map(float, element.find('%sIMLs' % NRML).text.split())
#
#            uhs = [lon, lat]
#            uhs.extend(imls)
#            values.append(uhs)
#
#    return metadata, periods, numpy.array(values)
OPTIONAL_PATHS = [("statistics", "statistics"),
                  ("gsimlt_path", "gsimTreePath"),
                  ("smlt_path", "sourceModelTreePath"),
                  ("quantile_value", "quantileValue")]


def parse_nrml_uhs_curves(nrml_uhs_map):
    """
    Reads the NRML file and returns the metadata (as a dictionary), the
    periods (as a numpy array) and the uhs values as an array of
    [lon, lat, uhs]
    """
    node_set = read_lazy(nrml_uhs_map, "IMLs")[0]
    # Read metadata

    metadata = {
        "smlt_path": node_set.attrib["sourceModelTreePath"],
        "investigation_time": float(node_set.attrib["investigationTime"]),
        "poe": float(node_set.attrib["poE"]),
        "gsimlt_path": node_set["gsimTreePath"]}
    for option, name in OPTIONAL_PATHS:
        if name in node_set.attrib:
            metadata[option] = node_set.attrib[name]
        else:
            metadata[option] = None

    periods = numpy.array(map(float, node_set.nodes[0].text.split()))
    values = []
    for node in node_set.nodes[1:]:
        subnodes = list(node.nodes)
        lon, lat = map(float, subnodes[0].nodes[0].text.split())
        uhs = [lon, lat]
        uhs.extend(map(float, subnodes[1].text.split()))
        values.append(uhs)
    return metadata, periods, numpy.array(values)


def plot_uhs(file_name_root, uhs, periods, metadata):
    """
    Takes the UHS data and produces a set of curves as pdf images in the
    output folder
    """
    os.mkdir(file_name_root)
    num_curves = numpy.shape(uhs)[0]
    if not metadata["statistics"]:
        metadata["statistics"] = ""
    for iloc, row in enumerate(uhs):
        print "Plotting curve %d of %d" % (iloc + 1, num_curves)
        fig = plt.figure(figsize=(7, 5))
        fig.set_tight_layout(True)
        plt.plot(periods, row[2:], 'bo-', linewidth=2.0)
        plt.xlabel("Period (s)", fontsize=14)
        plt.ylabel("Spectral Acceleration (g)", fontsize=14)
        plt.grid(b=True, color='0.66', linestyle="--")
        if row[0] < 0.0:
            long_ind = "W"
        else:
            long_ind = "E"
        if row[1] < 0.0:
            lat_ind = "S"
        else:
            lat_ind = "N"
        title_string_upper = "{:s} UHS with a {:s} PoE in {:s} Years\n".format(
             metadata["statistics"],
             str(metadata["poe"]),
             str(metadata["investigation_time"]))
        title_string_lower = "Location: {:.6f}{:s}, {:.6f}{:s}".format(
            numpy.abs(row[0]), long_ind, numpy.abs(row[1]), lat_ind)
        plt.title(title_string_upper + title_string_lower, fontsize=16)
        output_file = os.path.join(file_name_root,
            "UHS_{:.5f}{:s}_{:.5f}{:s}.pdf".format(row[0], long_ind,
            row[1], lat_ind))
        plt.savefig(output_file, dpi=300, format="pdf", papertype="a4")

def save_uhs_to_csv(nrml_uhs_file, file_name_root, plot_spectra=False):
    """
    Read hazard map in `nrml__hazard_map_file` and save to .csv file
    with root name `file_name_root`
    """
    output_file = '%s.csv' % file_name_root
    if os.path.isfile(output_file):
        raise ValueError('Output file already exists [%s].'
                         ' Please specify different name or '
                         'remove old file' % output_file)

    metadata, periods, values = parse_nrml_uhs_curves(nrml_uhs_file)

    header = ','.join(
        ['%s=%s' % (k, v) for k, v in metadata.items() if v is not None]
    )
    header = '# ' + header
    header += '\nlon,lat,'+','.join([str(p) for p in periods])

    f = open(output_file, 'w')
    f.write(header+'\n')
    numpy.savetxt(f, values, fmt='%g', delimiter=',')
    f.close()
    if plot_spectra:
        plot_uhs(file_name_root, values, periods, metadata)

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert UHS file from Nrml to .csv file.'
            'To run just type: python uhs_converter.py ' 
            '--input-file=/PATH/TO/INPUT_FILE '
            '--output-file=/PATH/TO/OUTPUT_FILE', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
                        help='path to NRML uhs file (Required)',
                        default=None,
                        required=True)
    flags.add_argument('--output-file',
                        help='path to output file, without file extension '
                             ' (Optional, default is root of NRML file)',
                        default=None
                       )
    flags.add_argument('--plot-spectra',
                       help="Plot the uniform hazard spectra to pdf (True) " 
                       "or not (False) - may take time for many hazard curves",
                       default=False)

    return parser


if __name__ == "__main__":
    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        output_file = \
            os.path.splitext(parser.parse_args().input_file)[0] \
            if args.output_file is None else args.output_file

        save_uhs_to_csv(args.input_file, output_file, args.plot_spectra)
    else:
        parser.print_usage()
