#!/usr/bin/env/python
# LICENSE
#
# Copyright (c) 2010-2013, GEM Foundation, G. Weatherill, M. Pagani,
# D. Monelli.
#
# The Hazard Modeller's Toolkit is free software: you can redistribute
# it and/or modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>
#
# DISCLAIMER
# 
# The software Hazard Modeller's Toolkit (hmtk) provided herein
# is released as a prototype implementation on behalf of
# scientists and engineers working within the GEM Foundation (Global
# Earthquake Model).
#
# It is distributed for the purpose of open collaboration and in the
# hope that it will be useful to the scientific, engineering, disaster
# risk and software design communities.
#
# The software is NOT distributed as part of GEM’s OpenQuake suite
# (http://www.globalquakemodel.org/openquake) and must be considered as a
# separate entity. The software provided herein is designed and implemented
# by scientific staff. It is not developed to the design standards, nor
# subject to same level of critical review by professional software
# developers, as GEM’s OpenQuake software suite.
#
# Feedback and contribution to the software is welcome, and can be
# directed to the hazard scientific staff of the GEM Model Facility
# (hazard@globalquakemodel.org).
#
# The Hazard Modeller's Toolkit (hmtk) is therefore distributed WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# The GEM Foundation, and the authors of the software, assume no
# liability for use of the software.
'''
Post-process uniform hazard data to convert OpenQuake uniform hazard
spectra from XML into different formats
'''

import argparse
import numpy as np
from lxml import etree
from scipy.io import savemat


xmlNRML = '{http://openquake.org/xmlns/nrml/0.4}'
xmlGML = '{http://www.opengis.net/gml}'


class UHSParser(object):
    """
    Base class for parsing the uniform hazard spectra nrml files
    """
    def __init__(self, filename):
        """
        """
        self.input_file = filename
        self.investigation_time = None
        self.poe = None
        self.statistics = None
        self.periods = None
        self.longitude = []
        self.latitude = []
        self.uhs = []

    def read_file(self):
        """
        Reads the input file
        """
        parse_args = dict(source=self.input_file)

        for _, element in etree.iterparse(**parse_args):
            if element.tag == "%suniformHazardSpectra" % xmlNRML:
                # UHS header info
                self.poe = float(element.attrib.get('poE'))
                self.statistics = element.attrib.get('statistics')
                self.investigation_time = float(
                    element.attrib.get('investigationTime'))
            elif element.tag == "%speriods" % xmlNRML:
                self.periods = map(float, (element.text).split())
            elif element.tag == "%spos" % xmlGML:
                location = map(float, (element.text).split())
                self.longitude.append(location[0])
                self.latitude.append(location[1])
            elif element.tag == "%sIMLs" % xmlNRML:
                self.uhs.append(map(float, (element.text).split()))
            else:
                pass

        self.longitude = np.array(self.longitude)
        self.latitude = np.array(self.latitude)
        self.uhs = np.array(self.uhs)
        self.periods = np.array(self.periods)


class UHStoCSV(UHSParser):
    """
    Returns the UHS as a simple csv file
    """
    def write_file(self, output_file):
        """
        Writes to a simple csv
        """
        npts = len(self.longitude)
        fid = open(output_file, 'wt')
        # Write header info
        print >> fid, 'Investigation Time, %6.2f, PoE, %.8e, Statistics, %s'\
            % (self.investigation_time, self.poe, self.statistics)
        # Print column headers
        period_string = ", ".join(['%7.3f' % per for per in self.periods])
        print >> fid, 'Longitude, Latitude, %s' % period_string
        for iloc in range(0, npts):
            uhs_string = ", ".join(['%.6e' % iml for iml in self.uhs[iloc, :]])
            print >> fid, '%10.5f, %10.5f, %s' % (self.longitude[iloc],
                                                  self.latitude[iloc],
                                                  uhs_string)
        fid.close()


class UHStoMatlabBinary(UHSParser):
    """
    Returns the UHS data as a data structure in a Matlab Binary
    """
    def write_file(self, output_file):
        """
        Writes to Matlab Binary
        """
        # Creates the Matlab Dictionary
        mat_dict = {'longitude': self.longitude,
                    'latitude': self.latitude,
                    'periods': self.periods,
                    'poe': self.poe,
                    'uhs': self.uhs,
                    'statistics': self.statistics,
                    'investigation_time': self.investigation_time}
        savemat(output_file, mat_dict, oned_as='column')


def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert Hazard Map file from Nrml to Something Readable'
                    'To run just type: python uhs_converter.py '
                    '--input-file=/PATH/INPUT_FILE_NAME '
                    '--output-file=/PATH/OUTPUT_FILE_NAME')
    parser.add_argument(
        '--input-file',
        help='path to NRML UHS file',
        default=None)
    parser.add_argument(
        '--output-file',
        help='path to output file root (without file extension)',
        default=None)
    parser.add_argument(
        '--file-type',
        help='File type as extension (i.e. "csv" or "mat")',
        default='csv')
    return parser

FILEMAP = {'csv': UHStoCSV,
           'mat': UHStoMatlabBinary}

if __name__ == "__main__":
    parser = set_up_arg_parser()
    args = parser.parse_args()
    # Read in file
    converter = FILEMAP[args.file_type](args.input_file)
    converter.read_file()
    # Write to file
    converter.write_file(args.output_file)
