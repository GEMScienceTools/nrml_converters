#!/usr/bin/env/python
# LICENSE
#
# Copyright (c) 2013, GEM Foundation, G. Weatherill, M. Pagani,
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
# The software nrml_convertes provided herein
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
# The nrml_convertes is therefore distributed WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# The GEM Foundation, and the authors of the software, assume no
# liability for use of the software.
'''
Post-process hazard calculation data to convert hazard curves into different
formats
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
xmlGML = '{http://www.opengis.net/gml}'

def parse_hazard_curve_node(element):
    '''
    Reads the hazard curve element to return the longitude, latitude and 
    poes
    '''
    for e in element.iter():        
        if e.tag == '%spos' % xmlGML:
            coords = str(e.text).split()
            lon = float(coords[0])
            lat = float(coords[1])
        elif e.tag == '%spoEs' % xmlNRML:
            poes = str(e.text).split()
            poes = map(float, poes)
        else:
            continue
    return lon, lat, poes


def parse_metadata(element):
    '''
    Returns the statistics, IMT and investigation time
    '''
    meta_info = {}
    meta_info['imt'] = element.attrib.get('IMT')
    if 'SA' in meta_info['imt']:
        meta_info['period'] = float(element.attrib.get('saPeriod'))
        meta_info['damping'] = float(element.attrib.get('saDamping'))
    else:
        meta_info['period'] = None
        meta_info['damping'] = None
    meta_info['statistics'] = element.attrib.get('statistics')
    investigation_time = float(element.attrib.get('investigationTime'))
    return meta_info, investigation_time

def get_imls(element):
    '''
    Returns the IML list
    '''
    imls = element.attrib.get('IMLs')
    return map(float, imls.split(','))


class HazardCurveParser(object):
    '''
    General class to parse the hazard curves
    '''
    def __init__(self, input_file):
        '''
        '''
        self.input_file = input_file
        self.longitude = []
        self.latitude = []
        self.curves = []
        self.meta_info = {}
        self.imls = None
        self.investigation_time = None

    def read_hazard_curve_file(self):
        '''
        Reads the data from the file
        '''
        parse_args = dict(source=self.input_file)

        for _, element in etree.iterparse(**parse_args):
            if element.tag == '%shazardCurves' % xmlNRML:
                self.meta_info, self.investigation_time = \
                    parse_metadata(element)
            elif element.tag == '%sIMLs' % xmlNRML:
                #imls = str(element.text).split()
                #print imls
                self.imls = map(float, str(element.text).split())
            elif element.tag == '%shazardCurve' % xmlNRML:
                lon, lat, poes = parse_hazard_curve_node(element)
                self.longitude.append(lon)
                self.latitude.append(lat)
                self.curves.append(poes)
            else:
                continue
        self.curves = np.array(self.curves)
        self.longitude = np.array(self.longitude)
        self.latitude = np.array(self.latitude)
   
    def adjust_curves_to_different_poe(self, new_it):
        '''
        Adjusts the hazard curves to a new probability of exceedence using the
        poisson model
        '''
        lamda = -np.log(1. - self.curves) / self.investigation_time
        self.curves = 1. - np.exp(-lamda * new_it)
        self.investigation_time = new_it
    
    def build_file_string(self, input_string):
        '''
        Builds a filename string including the IMT and investigation time
        '''
        if self.meta_info['imt'] == 'SA':
            marker = self.meta_info['imt'] + \
                ('%.2f' % self.meta_info['period']) + '_' + \
                ('damp%.1f' % self.meta_info['damping']) + '_'
        else:
            marker = self.meta_info['imt'] + '_'
           
        filename = input_string + '_' + marker + self.meta_info['statistics'] +\
            '_' + ('%.1f' %self.investigation_time) + '.'
        return filename


class HazardCurves2Csv(HazardCurveParser):
    '''
    Writes the hazard curve set to csv
    '''
    def write_file(self, ofile, append=False):
        '''
        Writes the output to csv
        '''
        fieldnames = ['Longitude', 'Latitude']
        iml_headers = ['%s' % iml for iml in self.imls]
        fieldnames.extend(iml_headers)
        if append:
            fout = open(ofile, 'a+')
        else:
            fout = open(ofile, 'wt')

        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        if not append:
            # If not appending then write the header row
            writer.writerow(OrderedDict((fn, fn) for fn in fieldnames))
        
        for iloc, curve in enumerate(self.curves):
            row_dict = {'Longitude': str(self.longitude[iloc]),
                        'Latitude': str(self.latitude[iloc])}
            for jloc, hdr in enumerate(iml_headers):
                row_dict[hdr] = '%.8e' % curve[jloc]
            writer.writerow(row_dict)
        fout.close()


class HazardCurves2MatlabBinary(HazardCurveParser):
    '''
    Writes the hazard curve to a matlab binary file
    '''
    def write_file(self, ofile):
        '''
        Writes the file
        '''
        # Creates the Matlab dictionary
        mat_dict = {'longitude': self.longitude,
                    'latitude': self.latitude,
                    'imls': self.imls,
                    'imt': self.meta_info['imt'],
                    'period': None,
                    'damping': None,
                    'curves': self.curves,
                    'statistics': self.meta_info['statistics'],
                    'investigation_time': self.investigation_time}
        if self.meta_info['imt'] == 'SA':
            mat_dict['period'] = self.meta_info['period']
            mat_dict['damping'] = self.meta_info['damping']
        elif self.meta_info['imt'] == 'PGA':
            mat_dict['period'] = 0.
            
            mat_dict['damping'] = np.nan
        else:
            pass
        # Save to binary
        savemat(ofile, mat_dict, oned_as='row')
        

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert Hazard Map file from Nrml to Something Readable'
            'To run just type: python hazard_map_converter.py ' 
            '--input-file=/PATH/INPUT_FILE_NAME '
            '--output-file=/PATH/OUTPUT_FILE_NAME.xml')
    parser.add_argument('--input-file',
                        help='path to NRML hazard curve file',
                        default=None)
    parser.add_argument('--output-file',
                        help='path to output file root (without file extension)',
                        default=None)
    parser.add_argument('--investigation-time',
                        help='Investigation time (years) for output',
                        default=None)
    parser.add_argument('--file-type',
                        help='File type as extension (i.e. "csv" or "mat")',
                        default='csv')
    return parser
    


FILE_MAP = {'csv': HazardCurves2Csv,
            'mat': HazardCurves2MatlabBinary}


if __name__ == "__main__":
    parser = set_up_arg_parser()
    args = parser.parse_args()
    # Read input file
    converter = FILE_MAP[args.file_type](args.input_file)
    converter.read_hazard_curve_file()
    if args.investigation_time:
        converter.adjust_curves_to_different_poe(
            float(args.investigation_time))
    filename = converter.build_file_string(args.output_file) + args.file_type
    converter.write_file(filename)





