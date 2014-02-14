#!/usr/bin/env python
# LICENSE
#
# Copyright (c) 2014, GEM Foundation, G. Weatherill, M. Pagani,
# D. Monelli.
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
# The software nrml_converters provided herein
# is released as a prototype implementation on behalf of
# scientists and engineers working within the GEM Foundation (Global
# Earthquake Model).
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
# The nrml_converters is therefore distributed WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# The GEM Foundation, and the authors of the software, assume no
# liability for use of the software.

'''
Convert (scenario based) GMF NRML file to .csv files.
'''
import os
import csv
import argparse
import numpy
from collections import OrderedDict
from lxml import etree

NRML='{http://openquake.org/xmlns/nrml/0.4}'


def parse_gmfs_file(file_name):
    """
    Parse NRML 0.4 GMF set file.
    """
    parse_args = dict(source=file_name)

    gmfs = OrderedDict()
    for _, element in etree.iterparse(**parse_args):
        if element.tag == '%sgmf' % NRML:
            imt = element.attrib['IMT']
            if imt == 'SA':
                damping = element.attrib['saDamping']
                period = element.attrib['saPeriod']
                imt = '%s(%s)' % (imt, period)
            gmf = []
            for node in element:
                gmv = float(node.attrib['gmv'])
                lon = float(node.attrib['lon'])
                lat = float(node.attrib['lat'])
                gmf.append([lon, lat, gmv])
            if imt in gmfs.keys():
                gmfs[imt].append(gmf)
            else:
                gmfs[imt] = [gmf]

    return gmfs

def save_gmfs_to_csv(gmfs, out_dir):
    """
    Save GMFs to .csv files
    """
    for imt in gmfs.keys():
        dir_name = '%s/GMFS_%s' % (out_dir, imt)
        os.makedirs(dir_name)
        for i, gmf in enumerate(gmfs[imt]):
            header = 'lon,lat,value'
            fname = '%s/gmf_%s.csv' % (dir_name, (i + 1))
            f = open(fname, 'w')
            f.write(header+'\n')
            numpy.savetxt(f, numpy.array(gmf), fmt='%g', delimiter=',')
            f.close()

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML scenario-based ground motion fields files'
            ' to .csv files. Inside the specified output directory, create '
            'subdirectories for each intensity measure type. A .csv '
            'file is created for each ground motion field.'
            'To run just type: python gmfset_converter.py '
            '--input-file=PATH_TO_GMF_NRML_FILE'
            '--output-dir=PATH_TO_OUTPUT_DIRECTORY', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
        help='path to gmf NRML file (Required)',
        default=None,
        required=True)
    flags.add_argument('--output-dir',
        help='path to output directory (Required, raise an error if it '
             'already exists)',
        default=None,
        required=True)

    return parser


if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        # create the output directory immediately. Raise an error if
        # it already exists
        os.makedirs(args.output_dir)

        gmfc = parse_gmfs_file(args.input_file)
        save_gmfs_to_csv(gmfc, args.output_dir)
    else:
        parser.print_usage()