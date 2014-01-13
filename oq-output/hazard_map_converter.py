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
Post-process hazard calculation data to join map data from tiles into
a single csv file
'''

import os
import csv
import argparse
import shapefile
import numpy as np
from lxml import etree

xmlNRML='{http://openquake.org/xmlns/nrml/0.4}'

def parse_hazard_map_node(element):
    """
    Parse Hazard Map Node element.
    """
    return float(element.get('lon')), float(element.get('lat')), \
        float(element.get('iml'))


    for e in element.iter():        
        if e.tag == '%spos' % xmlGML:
            coords = str(e.text).split()
            lon = float(coords[0])
            lat = float(coords[1])
            if e.tag == '%sIML' % xmlNRML:
                value = float(e.text)
    return lon,lat,value


class HazardMapParser(object):
    '''
    Base Class for converting hazard maps from nrml to other formats
    '''
    def __init__(self, input_file):
        '''
        '''
        self.input_file = input_file
        self.longitude = []
        self.latitude = []
        self.data = []
        self.meta_info = {}

    def read_hazard_map_file(self):
        """
        Parse NRML hazard map file.
        """
        parse_args = dict(source=self.input_file)

        self.longitude = []
        self.latitude = []
        self.data = []
        self.meta_info = {}
        for _, element in etree.iterparse(**parse_args):
            if element.tag == '%snode' % xmlNRML:
                # Get node
                lon,lat,value = parse_hazard_map_node(element)
                self.longitude.append(lon)
                self.latitude.append(lat)
                self.data.append(value)
            elif element.tag == '%shazardMap' % xmlNRML:
                # Get header info
                self.meta_info['imt'] = element.attrib.get('IMT')
                if 'SA' in self.meta_info['imt']:
                    # Get period and damping
                    self.meta_info['period'] = \
                        float(element.attrib.get('saPeriod'))
                    self.meta_info['damping'] = \
                        float(element.attrib.get('saDamping'))
                self.meta_info['inv_time'] = float(element.attrib.get(
                    'investigationTime'))
                self.meta_info['poe'] = float(element.attrib.get('poE'))
                self.meta_info['stats'] = element.attrib.get('statistics')
            else:
                pass
        
    def number_data_points(self):
        '''
        Gets the number of data points
        '''
        assert (len(self.longitude) == len(self.latitude)) and (
            len(self.longitude) == len(self.data))
        return len(self.longitude)

    def generate_meta_string(self):
        '''
        '''
        if len(self.meta_info.keys()) == 0:
            raise ValueError('Cannot generate info string when empty')
        
        if 'SA' in self.meta_info['imt']:
            imt = self.meta_info['imt'] + ('%.2f' % self.meta_info['period'])
        else:
            imt = self.meta_info['imt']

        return imt + '-poe' + str(self.meta_info['poe']) + \
            ('in%syr' % str(self.meta_info['inv_time'])) + \
            ('-%s' % self.meta_info['stats'])
        


class HazardMap2XYZ(HazardMapParser):
    '''
    Writes the Hazard Map into a simple xyz format - useful for CMT!
    '''
    def write_file(self, output_filename, append=False):
        '''
        Writes to file
        :param str output_filename
            Name of output file
        :param bool append:
            Allows for appending (True) or not (False)
        '''


        if len(self.data) == 0:
            raise ValueError('Cannot write empty data to file!')
        
        if append:
            # Open a file for appending
            f = open(output_filename, 'a+')
        else:
            # Create new file - add properties to name
            filestring = output_filename + '_' + \
                self.generate_meta_string() + '.xyz'
            f = open(filestring, 'wt')
    
        for i in range(0, self.number_data_points()):
            print >> f, '%s  %s  %s' % (str(self.longitude[i]), 
                                        str(self.latitude[i]),
                                        str(self.data[i]))
        f.close()


class HazardMap2csv(HazardMapParser):
    '''
    Writes to csv - similar to xyz but simply comma separated
    '''
    def write_file(self, output_filename, append=False):
        '''
        Write to file
        :param str output_filename
            Name of output file
        '''

        if len(self.data) == 0:
            raise ValueError('Cannot write empty data to file!')
        
        if append:
            # Open a file for appending
            f = open(output_filename, 'a+')
        else:
            # Create new file - add properties to name 
            filestring = output_filename + '_' + \
                self.generate_meta_string() + '.csv'
            f = open(filestring, 'wt')

        writer = csv.writer(f)
        # Write header!
        writer.writerow(('Longitude', 'Latitude', self.meta_info['imt']))
        for i in range(0, self.number_data_points()):
            writer.writerow((self.longitude[i], 
                             self.latitude[i], 
                             self.data[i]))
        f.close()


class HazardMap2Shapefile(HazardMapParser):
    '''
    Writes the hazard map to shapefile (as a set of points)
    '''

    def write_file(self, output_filename):
        '''
        Write to file
        :param str output_filename
            Name of output file
        '''
        if len(self.data) == 0:
            raise ValueError('Cannot write empty data to file!')
       
        output_filename = output_filename + '_' + self.generate_meta_string()
        
        w = shapefile.Writer(shapefile.POINT)
        # Set up file headers
        if 'SA' in self.meta_info['imt']:
            sa_period = str(round(self.meta_info['period'], 3))
            sa_period.replace('.', 'p')
            sa_string = 'SA_' + sa_period
            w.field(sa_string, 'N', 15, 10)
        else:
            w.field(self.meta_info['imt'], 'N', 15, 10)
        # Generate shapefile database
        for i in range(0, self.number_data_points()):
            w.point(self.longitude[i], self.latitude[i], 0, 0)
            w.record(round(self.data[i], 10))
        # Write to shapefile
        w.save(output_filename)


MAP_TYPE = {'xyz': HazardMap2XYZ,
            'csv': HazardMap2csv,
            'shp': HazardMap2Shapefile}


def maps_to_single_file(map_list, output_filename, fmt='csv'):
    '''
    Parses a list of maps to a single output file
    '''
    if not fmt in ['csv', 'xyz']:
        raise IOError('Output file type %s not supported!' % fmt)

    for haz_map in map_list:
        converter = MAP_TYPE[fmt](haz_map)
        converter.read_hazard_map_file()
        converter.writer_file(output_filename, append=True)
    

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
                        help='path to NRML hazard map file',
                        default=None)
    parser.add_argument('--output-file',
                        help='path to output csv file',
                        default=None)

    return parser


if __name__=="__main__":
    # As a default process, just convert to csv
    parser = set_up_arg_parser()
    args = parser.parse_args()
    converter = HazardMap2csv(args.input_file)
    converter.read_hazard_map_file()
    converter.writer_file(args.output_file)

