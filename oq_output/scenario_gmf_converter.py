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
from openquake.commonlib.nrml import read
from hazard_map_converter import atkinson_kaka_2007_rsa2mmi, AK2007

NRML='{http://openquake.org/xmlns/nrml/0.5}'


def parse_gmfs_file(file_name):
    """
    Parses the NRML 0.4 GMF set file
    """

    node_set = read(file_name, stop="node")[0]
    gmf_set = []
    gmfs = []
    imt_list = []
    rup_id_list = []

    for gmfset_node in node_set:
        scenario_keys = []
        for gmf_node in gmfset_node:    
            # Parse field
            imt = gmf_node["IMT"]
            if "SA" in imt:
                imt = "SA({:s})".format(str(gmf_node["saPeriod"]))
            if not imt in imt_list:
                if len(rup_id_list):
                    gmfs.append((imt_list[-1], OrderedDict(rup_id_list)))
                    rup_id_list = []
                imt_list.append(imt)

            rup_id = gmf_node["ruptureId"]
            if not rup_id in scenario_keys:
                scenario_keys.append(rup_id)
            gmf = [[val["lon"], val["lat"], val["gmv"]] for val in gmf_node]
            rup_id_list.append((rup_id, numpy.array(gmf)))
        gmfs.append((imt, OrderedDict(rup_id_list)))
        gmfs = OrderedDict(gmfs)
        # Re-order to give in terms of scenario ids
        scenarios = []
        for rup_id in scenario_keys:
            # Determine the IMTs available for the scenario

            scenario_gmfs = {
                "IMT": ["lon", "lat"] + imt_list,
                "GMFs": numpy.zeros([len(gmf), len(imt_list) + 2])}
            for i, imt in enumerate(imt_list):
                scenario_gmfs["GMFs"][:, i + 2] = gmfs[imt][rup_id][:, 2]
            scenario_gmfs["GMFs"][:, :2] = gmfs[imt][rup_id][:, :2]
            scenarios.append((rup_id, scenario_gmfs))
        gmf_set.append((gmfset_node["stochasticEventSetId"],
                        OrderedDict(scenarios)))
        gmfs = []
        imt_list = []
    return OrderedDict(gmf_set)


def save_gmfs_to_csv(gmfs, out_dir):
    """
    Save GMFs to .csv files
    """
    for ses_id in gmfs.keys():
        dir_name = os.path.join(out_dir, "SES{:s}".format(ses_id))
        os.makedirs(dir_name)
        # One file per rupture scenario
        for rup_id in gmfs[ses_id].keys():
            fname = os.path.join(dir_name,
                "ses{:s}_{:s}.csv".format(ses_id,
                                          rup_id.replace("-", "_")))
            header = ",".join(gmfs[ses_id][rup_id]["IMT"])
            f = open(fname, "w")
            f.write(header+"\n")
            numpy.savetxt(f, gmfs[ses_id][rup_id]["GMFs"], fmt="%g",
                          delimiter=",")
            f.close()

#    for imt in gmfs.keys():
#        dir_name = '%s/GMFS_%s' % (out_dir, imt)
#        os.makedirs(dir_name)
#        for i, gmf in enumerate(gmfs[imt]):
#            header = 'lon,lat,value'
#            fname = '%s/gmf_%s.csv' % (dir_name, (i + 1))
#            f = open(fname, 'w')
#            f.write(header+'\n')
#            numpy.savetxt(f, numpy.array(gmf), fmt='%g', delimiter=',')
#            f.close()

def magic_flipud(values, llat, ulat):
    """

    """
    dlat = ulat - llat
    lat_diff = values[:, 1] - llat
    new_lat = ulat - lat_diff
    values[:, 1] = new_lat
    return values


def save_gmfs_to_netcdf(gmfs, out_dir, to_mmi=False, resolution="10k",
        cleanup=False, magic=False):
    """
    Save ground motion fields to netCDF files for use in GMT and/or QGis
    """
    for ses_id in gmfs.keys():
        dir_name = os.path.join(out_dir, "SES{:s}".format(ses_id))
        os.makedirs(dir_name)
        for rup_id in gmfs[ses_id].keys():
            #lon = gmfs[ses_id][rup_id]["GMFs"][:, 0]
            lat = gmfs[ses_id][rup_id]["GMFs"][:, 1]
            #gmvs = gmfs[ses_id][rup_id]["GMFs"][:, 2:]
            ogmf = gmfs[ses_id][rup_id]["GMFs"]
           
            for i, imt in enumerate(gmfs[ses_id][rup_id]["IMT"][2:]):
                if to_mmi:
                    if imt in AK2007.keys():
                        get_mmi = True
                    else:
                        print "Cannot convert %s to MMI" % imt
                        get_mmi = False
                else:
                    get_mmi = False
                f_stem = os.path.join(dir_name, 
                    "ses{:s}_{:s}_{:s}".format(ses_id,
                                               rup_id.replace("-", "_"),
                                               imt))
                if magic:
                    ogmf = magic_flipud(ogmf, numpy.min(lat), numpy.max(lat))
                if get_mmi:
                    ogmf[:, 2], sigma = atkinson_kaka_2007_rsa2mmi(imt,
                                                                   ogmf[:, 2])
                    f_stem = f_stem.replace(imt, "{:s}2MMI".format(imt))
                f = open(f_stem + ".xyz", "w")
                numpy.savetxt(f, ogmf, fmt="%g", delimiter=" ")
                f.close()
                # Generate GMT xyz2grd string
                llon = numpy.min(ogmf[:, 0])
                ulon = numpy.max(ogmf[:, 0])
                llat = numpy.min(ogmf[:, 1])
                ulat = numpy.max(ogmf[:, 1])
                istring = " -I%s/%s" % (resolution, resolution)
                rstring = " -R" + "{:.5f}/".format(llon) + \
                    "{:.5f}/".format(ulon) +\
                    "{:.5f}/".format(llat) + "{:.5f}".format(ulat)
                if "SA" in imt:
                    mod1 = f_stem.replace("(", "\(")
                    f_stem = mod1.replace(")", "\)")
                    
                command_string = ("xyz2grd %s -G%s" % (f_stem + ".xyz",
                                  f_stem + ".NC")) + istring + rstring
                os.system(command_string)
        if cleanup:
            os.system("find . -name '*.xyz' -type f -delete")

#    for imt in gmfs.keys():
#        dir_name = '%s/GMFS_%s' % (out_dir, imt)
#        os.makedirs(dir_name)
#        if to_mmi:
#            if imt in AK2007.keys():
#                get_mmi = True
#            else:
#                print "Cannot convert %s to MMI" % imt
#                get_mmi = False
#        else:
#            get_mmi = False
#            
#        for i, gmf in enumerate(gmfs[imt]):
#            f_stem = '%s/gmf_%s' % (dir_name, (i + 1))
#            ogmf = numpy.array(gmf)
#            if magic:
#                ogmf = magic_flipud(ogmf,
#                                    numpy.min(ogmf[:, 1]),
#                                    numpy.max(ogmf[:, 1]))
#            if get_mmi:
#                ogmf[:, 2], sigma = atkinson_kaka_2007_rsa2mmi(imt, ogmf[:, 2])
#            f = open(f_stem + ".xyz", "w")
#            numpy.savetxt(f, ogmf, fmt="%g", delimiter=" ")
#            f.close()
#            # Generate GMT xyz2grd string
#            llon = numpy.min(ogmf[:, 0])
#            ulon = numpy.max(ogmf[:, 0])
#            llat = numpy.min(ogmf[:, 1])
#            ulat = numpy.max(ogmf[:, 1])
#            istring = " -I%s/%s" % (resolution, resolution)
#            rstring = " -R" + "{:.5f}/".format(llon) + "{:.5f}/".format(ulon) +\
#                "{:.5f}/".format(llat) + "{:.5f}".format(ulat)
#            if "SA" in imt:
#                mod1 = f_stem.replace("(", "\(")
#                f_stem = mod1.replace(")", "\)")
#                
#            command_string = ("xyz2grd %s -G%s" % (f_stem + ".xyz",
#                              f_stem + ".NC")) + istring + rstring
#            os.system(command_string)
#    if cleanup:
#        os.system("find . -name '*.xyz' -type f -delete")
        

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

    flags.add_argument('--magic',
        help="Flips raster in vertical",
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
        gmfc = parse_gmfs_file(args.input_file)
        if args.to_netcdf:
            save_gmfs_to_netcdf(gmfc, args.output_dir, args.to_mmi,
                                args.spacing, args.cleanup, args.magic)    
        else:
            save_gmfs_to_csv(gmfc, args.output_dir)
    else:
        parser.print_usage()
