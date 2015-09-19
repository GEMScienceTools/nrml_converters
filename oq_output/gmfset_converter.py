#!/usr/bin/env python
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
Convert GMF NRML file to .csv files.
'''

import os
import csv
import argparse
import numpy
from collections import OrderedDict
#from lxml import etree
from openquake.commonlib.nrml import read_lazy
from hazard_map_converter import atkinson_kaka_2007_rsa2mmi, AK2007

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
    node_set = read_lazy(file_name, "node")
    gmfss = []
    for element in node_set:
        if "gmfSet" in element.tag:
            gmfss.append(parse_gmf_set(element))
        elif "gmfCollection" in element.tag:
            gmfc = GmfCollection(element.attrib["sourceModelTreePath"],
                                 element.attrib["gsimTreePath"],
                                 None)
        else:
            pass
    gmfc.gmfss = gmfss
    return gmfc



def parse_gmf_set(gmfset_node):
    """
    Parse NRML 0.4 GMF set element.
    """
    stochasticEventSetId = gmfset_node.attrib['stochasticEventSetId']
    investigationTime = float(gmfset_node.attrib['investigationTime'])
    print 'Processing gmf set %s' % stochasticEventSetId
    gmfs = OrderedDict()
    for node in gmfset_node.nodes:
        gmf = parse_gmf(node)
        IMT = gmf.IMT if gmf.IMT != 'SA' else '%s(%s)' % (gmf.IMT, gmf.saPeriod)

        if IMT in gmfs.keys():
            gmfs[IMT].append((gmf.ruptureId, gmf.values))
        else:
            gmfs[IMT] = [(gmf.ruptureId, gmf.values)]

    return GmfSet(stochasticEventSetId, investigationTime, gmfs)    
    

def parse_gmf(gmf_node):
    """
    Parse NRML 0.4 GMF element.
    """
    ruptureId = gmf_node.attrib["ruptureId"]
    IMT = gmf_node.attrib["IMT"]
    saDamping = None
    saPeriod = None
    if "IMT" == "SA":
        saDamping = gmf_node.attrib["saDamping"]
        saPeriod = gmf_node.attrib["saPeriod"]
    values = []
    for node in gmf_node.nodes:
        values.append([float(node.attrib["lon"]),
                       float(node.attrib["lat"]),
                       float(node.attrib["gmv"])])
    return GMF(IMT, ruptureId, numpy.array(values), saDamping, saPeriod)


def save_gmfs_to_csv(gmf_collection, out_dir):
    """
    Save GMFs to .csv files
    """
    for gmf_set in gmf_collection.gmfss:
        for imt in gmf_set.gmfs.keys():
            dir_name = '%s/StochasticEventSet_%s_%s' % \
                (out_dir, gmf_set.stochasticEventSetId, imt)
            os.makedirs(dir_name)

            map_values = []
            for rup_id, values in gmf_set.gmfs[imt]:
                print 'saving file for rupture %s' % rup_id
                header = '# smltp=%s, gsimltp=%s' % \
                    (gmf_collection.sm_tp, gmf_collection.gsim_tp)
                header += '\n# IMT=%s' % imt
                header += '\nlon,lat,gmf_value'
                fname = '%s/%s.csv' % (dir_name, rup_id)
                f = open(fname, 'w')
                f.write(header+'\n')
                numpy.savetxt(f, values, fmt='%5.2f,%5.2f,%g')
                f.close()

def save_gmfs_to_netcdf(gmf_collection, out_dir, to_mmi=False, spacing="10k",
        cleanup=False):
    """
    Exports the ground motion fields to NetCDF format
    """
    for gmf_set in gmf_collection.gmfss:
        for imt in gmf_set.gmfs.keys():
            dir_name = '%s/StochasticEventSet_%s_%s' % \
                (out_dir, gmf_set.stochasticEventSetId, imt)
            os.makedirs(dir_name)
            if to_mmi:
                if imt in AK2007.keys():
                    get_mmi = True
                else:
                    print "Cannot convert %s to MMI" % imt
                    get_mmi = False
            else:
                get_mmi = False

            map_values = []
            for rup_id, values in gmf_set.gmfs[imt]:
                rup_id = rup_id.replace("|","_")
                rup_id = rup_id.replace("=","")

                print 'saving file for rupture %s' % rup_id
                if get_mmi:
                    values[:, 2], sigma = atkinson_kaka_2007_rsa2mmi(
                        imt, 
                        values[:, 2])
                f_stem = '%s/%s' % (dir_name, rup_id)

                f = open(f_stem + ".xyz", 'w')
                numpy.savetxt(f, values, fmt='%5.2f %5.2f %g')
                f.close()
                # Generate GMT xyz2grd string
                values = numpy.array(values)
                llon = numpy.min(values[:, 0])
                ulon = numpy.max(values[:, 0])
                llat = numpy.min(values[:, 1])
                ulat = numpy.max(values[:, 1])
                istring = " -I%s/%s" % (spacing, spacing)
                rstring = " -R" + "{:.5f}/".format(llon) + \
                    "{:.5f}/".format(ulon) + "{:.5f}/".format(llat) +\
                    "{:.5f}".format(ulat)
                if "SA" in imt:
                    mod1 = f_stem.replace("(", "\(")
                    f_stem = mod1.replace(")", "\)")
                
                command_string = ("xyz2grd %s -G%s" % (f_stem + ".xyz",
                                  f_stem + ".NC")) + istring + rstring
                os.system(command_string)
    if cleanup:
        os.system("find . -name '*.xyz' -type f -delete")



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
            '--output-dir=PATH_TO_OUTPUT_DIRECTORY', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
        help='path to gmf NRML file (Required)',
        default=None,
        required=True)
    flags.add_argument('--output-dir',
        help='path to output directory (Required, raise an error if it already exists)',
        default=None,
        required=True)
    flags.add_argument('--to-netcdf',
        help='Converts files to netcdf format for use with GMT/QGis',
        default=False,
        required=False)
    flags.add_argument('--spacing',
        help="Approximate spacing (km) of interpolated mesh",
        default="10k",
        required=False)
    flags.add_argument('--to-mmi',
        help="Convert the ground motion values to MMI using the model of "
             "Atkinson & Kaka (2007) - Note this will only apply to PGA, PGV "
             "SA(0.3), SA(1.0), SA(2.0)",
        default=False,
        required=False)
    flags.add_argument('--cleanup',
        help="Removes .xyz files after creation of NetCDF",
        default=False,
        required=False)
    return parser


if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        # create the output directory immediately. Raise an error if
        # it already exists
        os.makedirs(args.output_dir)

        gmfc = parse_gmfc_file(args.input_file)
        if args.to_netcdf:
            save_gmfs_to_netcdf(gmfc, args.output_dir, args.to_mmi,
                args.spacing, args.cleanup)
        else:
            save_gmfs_to_csv(gmfc, args.output_dir)
    else:
        parser.print_usage()
