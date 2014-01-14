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
Convert GMF NRML file to .csv or .mat files.
'''

import os
import csv
import argparse
import numpy
from lxml import etree
from scipy.io import savemat

NRML='{http://openquake.org/xmlns/nrml/0.4}'


class GMF(object):
    """
    GMF
    """
    def __init__(self, IMT, ruptureId, values, saDamping, saPeriod):
        self.IMT = IMT
        self.ruptureId = ruptureId
        self.values = values
        self.saDamping = saDamping
        self.saPeriod = saPeriod


class GmfSet(object):
    """
    GMF Set.
    """
    def __init__(self, stochasticEventSetId, investigationTime, gmfs):
        self.stochasticEventSetId = stochasticEventSetId
        self.investigationTime = investigationTime
        self.gmfs = gmfs


class GmfCollection(object):
    """
    GMF collection.
    """
    def __init__(self, sm_tp, gsim_tp, gmfss):
        self.sm_tp = sm_tp
        self.gsim_tp = gsim_tp
        self.gmfss = gmfss


def parse_gmfc_file(file_name):
    """
    Parse NRML 0.4 GMF collection file.
    """
    parse_args = dict(source=file_name)

    for _, element in etree.iterparse(**parse_args):
        if element.tag == '%sgmfCollection' % NRML:
            gmfc = parse_gmfc_collection(element)
            element.clear()

    return gmfc

def parse_gmfc_collection(element):
    """
    Parse NRML 0.4 GMF collection element.
    """
    sm_tp = element.attrib['sourceModelTreePath']
    gsim_tp = element.attrib['gsimTreePath']

    gmfss = []
    for e in element.iterchildren():
        gmfss.append(parse_gmf_set(e))
        e.clear()

    return GmfCollection(sm_tp, gsim_tp, gmfss)

def parse_gmf_set(element):
    """
    Parse NRML 0.4 GMF set element.
    """
    stochasticEventSetId = element.attrib['stochasticEventSetId']
    investigationTime = element.attrib['investigationTime']
    print 'Processing gmf set %s' % stochasticEventSetId

    gmfs = {}
    for e in element.iterchildren():
        gmf = parse_gmf(e)
        IMT = gmf.IMT if gmf.IMT != 'SA' else '%s(%s)' % (gmf.IMT, gmf.saPeriod)

        if IMT in gmfs.keys():
            gmfs[IMT].append((gmf.ruptureId, gmf.values))
        else:
            gmfs[IMT] = [(gmf.ruptureId, gmf.values)]

        e.clear()

    return GmfSet(stochasticEventSetId, investigationTime, gmfs)

def parse_gmf(element):
    """
    Parse NRML 0.4 GMF element.
    """
    ruptureId = element.attrib['ruptureId']
    IMT = element.attrib['IMT']
    saDamping = None
    saPeriod = None
    if IMT == 'SA':
        saDamping = element.attrib['saDamping']
        saPeriod = element.attrib['saPeriod']

    print 'Processing gmf %s' % ruptureId
    lons = []
    lats = []
    values = []
    for e in element.iterchildren():
        lons.append(float(e.attrib['lon']))
        lats.append(float(e.attrib['lat']))
        values.append(float(e.attrib['gmv']))
        e.clear()
    values = numpy.array([lons, lats, values]).T

    return GMF(IMT, ruptureId, values, saDamping, saPeriod)

def save_gmfs_to_csv(gmf_collection, out_dir):
    """
    Save GMFs to .csv files
    """
    for gmf_set in gmf_collection.gmfss:
        for imt in gmf_set.gmfs.keys():
            dir_name = '%s/StochasticEventSet_%s_%s' % \
                (out_dir, gmf_set.stochasticEventSetId, imt)
            os.makedirs(dir_name)

            for rup_id, values in gmf_set.gmfs[imt]:
                print 'saving file for rupture %s' % rup_id
                header = 'smltp=%s, gsimltp=%s' % \
                    (gmf_collection.sm_tp, gmf_collection.gsim_tp)
                header += '\nIMT=%s' % imt
                fname = '%s/%s.csv' % (dir_name, rup_id)
                numpy.savetxt(fname, values, header=header, fmt='%5.2f,%5.2f,%g')

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML ground motion fields file to .csv files. '
            'Inside the specified output directory, create subdirectories '
            'for each intensity measure type and stochastic event set. A .csv '
            'file is created for each ground motion field.'
            'To run just type: python gmfset_converter.py '
            '--input-file=PATH_TO_GMF_NRML_FILE'
            '--output-dir=PATH_TO_OUTPUT_DIRECTORY')
    parser.add_argument('--input-file',
        help='path to gmf NRML file',
        default=None,
        required=True)
    parser.add_argument('--output-dir',
        help='path to output directory (default is the current working directory)',
        default=os.getcwd())

    return parser


if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        gmfc = parse_gmfc_file(args.input_file)
        save_gmfs_to_csv(gmfc, args.output_dir)
        #converter = FILE_MAP[args.file_type](args.input_file)
        #converter.read_gmfset_file()
        #converter.write_file(args.output_file)
    else:
        parser.print_usage()
