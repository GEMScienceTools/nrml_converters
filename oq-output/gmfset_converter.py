#!/usr/bin/env/python
# LICENSE
#
# Copyright (c) 2010-2013, GEM Foundation, G. Weatherill, M. Pagani,
# D. Monelli.
#
# The Hazard Modeller's Toolkit is free software: you can redistribute
# it and/or modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>
#
# DISCLAIMER
#
# The software Hazard Modeller's Toolkit (hmtk) provided herein
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
# The Hazard Modeller's Toolkit (hmtk) is therefore distributed WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# The GEM Foundation, and the authors of the software, assume no
# liability for use of the software.

'''
Post-process hazard calculation data process a set of ground motion fields
from xml format to other formats
'''

import os
import csv
import argparse
import shapefile
import numpy as np
from lxml import etree
from collections import OrderedDict
from scipy.io import savemat

xmlNRML='{http://openquake.org/xmlns/nrml/0.4}'


def parse_gmf_node(element):
    """
    Parse Hazard Map Node element.
    """
    return float(element.attrib.get('lon')), float(element.attrib.get('lat')),\
        float(element.attrib.get('gmv'))

class GMFSetParser(object):
    """
    Base class for converting ground motion field sets from xml
    """
    def __init__(self, input_file):
        """
        """
        self.input_file = input_file
        self.longitude = []
        self.latitude = []
        self.gmfset = OrderedDict()

    def read_gmfset_file(self):
        """
        Loads the input file into memory
        """
        parse_args = dict(source=self.input_file)
        for _, element in etree.iterparse(**parse_args):
            if element.tag == '%sgmf' % xmlNRML:
                # Gets the event set
                self.process_gmf_element(element)
            else:
                pass

        for key in self.gmfset.keys():
            self.gmfset[key] = np.array(self.gmfset[key])
            self.gmfset[key] = self.gmfset[key].T

    def process_gmf_element(self, element):
        """
        Processes the element to get the ground motion field
        """
        imt = element.attrib.get("IMT")
        if imt == "SA":
            # Concatenate into a single key
            imt_key = "SA" + "_" + element.attrib.get("saPeriod") + "_" +\
                element.attrib.get("saDamping")
        else:
            imt_key = imt
        # Get the field 
        lons = []
        lats = []
        gmvs = []
        for gmf_element in element.getiterator():
            if gmf_element.tag == ("%snode" % xmlNRML):
                lons.append(float(gmf_element.attrib.get("lon")))
                lats.append(float(gmf_element.attrib.get("lat")))
                gmvs.append(float(gmf_element.attrib.get("gmv")))
        if len(self.longitude) == 0:
            self.longitude = np.array(lons)
            self.latitude = np.array(lats)

        if not imt_key in self.gmfset.keys():
            # IMT is not in dictionary
            self.gmfset[imt_key] = [gmvs]
        else:
            # IMT is in dictionary - so parse data
            self.gmfset[imt_key].append(gmvs)



class GMFSetToCsv(GMFSetParser):
    """
    Writes the GMFSet to a Csv File
    """
    def write_file(self, output_file):
        """
        Writes to the output file
        """
        fid = open(output_file, 'a+')
        for key in self.gmfset.keys():
            print >> fid, "Longitude, Latitude, %s" % key
            
            data = np.column_stack([self.longitude, 
                                    self.latitude, 
                                    self.gmfset[key]])
            np.savetxt(fid, data, fmt='%.8e', delimiter=',')
        fid.close()



class GMFSetToMultipleCsv(GMFSetParser):
    """
    Writes the GMFSet to a set of Muliple csv files
    """
    def write_file(self, file_stem):
        """
        Writes to a set of files with a common stem. Each file is the stem
        that is appended with the GMF ID
        """

        for key in self.gmfset.keys():
            filename = file_stem + "_" + key + ".csv"
            data = np.column_stack([self.longitude, 
                                    self.latitude, 
                                    self.gmfset[key]])
            fid = open(filename, 'wt')
            print >> fid, "Longitude, Latitude, %s" % key
            np.savetxt(fid, data, fmt="%.8e", delimiter=",")
            fid.close()
           
class GMFSetToMatlabBinary(GMFSetParser):
    """
    Write the GMFSet to a Matlab Binary
    """
    def write_file(self, output_file):
        """
        """
        mat_dict = {'longitude': self.longitude,
                    'latitude': self.latitude}
        for key in self.gmfset.keys():
            mat_dict[key] = self.gmfset[key]
        savemat(output_file, mat_dict, oned_as='column')


FILE_MAP = {'csv': GMFSetToCsv,
            'multi-csv': GMFSetToMultipleCsv,
            'mat': GMFSetToMatlabBinary}

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert GMF Set file from Nrml to Something Readable'
            'To run just type: python gmfset_converter.py ' 
            '--input-file=/PATH/INPUT_FILE_NAME '
            '--output-file=/PATH/OUTPUT_FILE_NAME.xml'
            '--file-type=csv | multi-csv | mat')
    parser.add_argument('--input-file',
        help='path to gmfset file',
        default=None)
    parser.add_argument('--output-file',
        help='path to output file (without file extension for multi-csv)',
        default=None)
    parser.add_argument('--file-type',
        help='File type as extension (i.e. "csv", "multi-csv" or "mat")',
        default='csv')
    return parser


if __name__ == "__main__":
    parser = set_up_arg_parser()
    args = parser.parse_args()
    # Read input file
    converter = FILE_MAP[args.file_type](args.input_file)
    converter.read_gmfset_file()
    converter.write_file(args.output_file)
