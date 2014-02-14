#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
# LICENSE
#
# Copyright (c) 2010-2014, GEM Foundation, V. Silva
#
# The Risk Modeller's Toolkit is free software: you can redistribute
# it and/or modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>
#
# DISCLAIMER
# 
# The software Risk Modeller's Toolkit (rmtk) provided herein
# is released as a prototype implementation on behalf of
# scientists and engineers working within the GEM Foundation (Global
# Earthquake Model).
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
# directed to the risk scientific staff of the GEM Model Facility
# (risk@globalquakemodel.org).
#
# The Risk Modeller's Toolkit (rmtk) is therefore distributed WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# The GEM Foundation, and the authors of the software, assume no
# liability for use of the software.
# -*- coding: utf-8 -*-
'''
Post-process risk calculation data to convert loss curves into different
formats
'''

import os
import csv
import argparse
import numpy as np
from lxml import etree
from collections import OrderedDict

xmlNRML='{http://openquake.org/xmlns/nrml/0.4}'
xmlGML = '{http://www.opengis.net/gml}'

def parse_single_loss_curve(element):
    '''
    Reads the loss curve element to return the longitude, latitude and 
    poes and losses
    '''
    for e in element.iter(): 
        ref = element.attrib.get('assetRef')
        if e.tag == '%spos' % xmlGML:
            coords = str(e.text).split()
            lon = float(coords[0])
            lat = float(coords[1])
        elif e.tag == '%spoEs' % xmlNRML:
            poes = str(e.text).split()
            poes = map(float, poes)
        elif e.tag == '%slosses' % xmlNRML:
            losses = str(e.text).split()
            losses = map(float, losses)
        else:
            continue
    return lon, lat, ref, poes, losses


def parse_metadata(element):
    '''
    Returns the statistics
    '''
    meta_info = {}
    meta_info['statistics'] = element.attrib.get('statistics')
    meta_info['quantile_value'] = element.attrib.get('quantileValue')
    meta_info['investigationTime'] = element.attrib.get('investigationTime')
    meta_info['unit'] = element.attrib.get('unit')
    
    return meta_info

def LossCurveParser(input_file):

    refs = []
    longitude = []
    latitude = []
    losses = []
    poes = []
    meta_info = {}

    for _, element in etree.iterparse(input_file):
        if element.tag == '%slossCurves' % xmlNRML:
            meta_info = parse_metadata(element)
        elif element.tag == '%slossCurve' % xmlNRML:
            lon, lat, ref, poe, loss = parse_single_loss_curve(element)
            longitude.append(lon)
            latitude.append(lat)
            refs.append(ref)
            poes.append(poe)
            losses.append(loss)
        else:
            continue
    longitude = np.array(longitude)
    latitude = np.array(latitude)
    
    return refs, longitude, latitude, poes, losses
   
def LossCurves2Csv(nrml_loss_curves):
    '''
    Writes the Loss curve set to csv
    '''
    refs, longitude, latitude, poes, losses = LossCurveParser(nrml_loss_curves)
    output_file = open(nrml_loss_curves.replace('xml','csv'),'w')
    for iloc in range(len(refs)):
        print len(poes)
        poes_list = ','.join(map(str, poes[iloc]))
        losses_list = ','.join(map(str, losses[iloc]))
        output_file.write(str(refs[iloc])+','+str(longitude[iloc])+','+str(latitude[iloc])+','+poes_list+','+losses_list+'\n')
    output_file.close()


def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML loss curves file to tab delimited '
            ' .txt files. Inside the specified output directory, create a .txt '
            'file for each stochastic event set.'
            'To run just type: python parse_loss_curves.py '
            '--input-file=PATH_TO_LOSS_CURVE_NRML_FILE ', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
        help='path to loss curves NRML file (Required)',
        default=None,
        required=True)

    return parser

if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        # create the output directory immediately. Raise an error if
        # it already exists
        LossCurves2Csv(args.input_file)