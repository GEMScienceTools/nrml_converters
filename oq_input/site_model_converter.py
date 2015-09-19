#!/usr/bin/env python
# LICENSE
#
# Copyright (c) 2014, GEM Foundation, G. Weatherill, M. Pagani, D. Monelli.
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

"""
Convert a site model from csv to Site_Model.xml and vice-versa

"""
import ast
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
from collections import OrderedDict

from openquake.commonlib.node import read_nodes, LiteralNode
from openquake.commonlib import nrml
from openquake.commonlib.valid import site_param


def xml_to_csv(input_xml, output_csv):
    """
    Parses the site model from an input xml file to a headed csv file
    """
    # Read in from XML
    sites = read_nodes(input_xml, lambda el: el.tag.endswith("site"),
                       nrml.nodefactory["siteModel"])
    fid = open(output_csv, "w")
    print >> fid, "%s" % "longitude,latitude,vs30,vs30Type,z1pt0,z2pt5,backarc"
    for site in sites:
        if "backarc" in site.attrib:
            if ast.literal_eval(site.attrib["backarc"]):
                site.attrib["backarc"] = 1
            else:
                site.attrib["backarc"] = 0

        else:
            site.attrib["backarc"] = 0

        if site["vs30Type"] == "measured":
            vs30_type = 1
        else:
            vs30_type = 0

        print >> fid, "%s" % ",".join([
            site["lon"], site["lat"], site["vs30"],
            str(vs30_type), site["z1pt0"], site["z2pt5"],
            str(site["backarc"])])
    fid.close()

def csv_to_xml(input_csv, output_xml):
    """
    Parses the site model from an input (headed) csv file to an output xml
    """
    data = np.genfromtxt(input_csv, delimiter=",", names=True)
    site_nodes = []
    for i in range(0, len(data)):
        site_attrib = [("lon", str(data["longitude"][i])),
                       ("lat", str(data["latitude"][i])),
                       ("vs30", str(data["vs30"][i])),
                       ("vs30Type", str(bool(data["vs30Type"][i]))),
                       ("z1pt0", str(data["z1pt0"][i])),
                       ("z2pt5", str(data["z2pt5"][i]))]
        if "backarc" in data:
            site_attrib.append(("backarc", str(bool(data["backarc"][i]))))
        else:
            site_attrib.append(("backarc", "False"))
        site_nodes.append(LiteralNode("site",
                                      OrderedDict(site_attrib),
                                      nodes=None))
    site_model = LiteralNode("siteModel", nodes=site_nodes)
    with open(output_xml, "w") as fid:
        nrml.write([site_model], fid, "%s")


def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """

    description = ('Convert Site Model from csv to xml and '
                   'vice versa.\n\nTo convert from csv to xml: '
                   '\npython site_model_converter.py '
                   '--input-csv-file PATH_TO_SITE_MODEL_CSV_FILE. '
                   '--output-xml-file PATH_TO_OUTPUT_FILE'
                   '\n\nTo convert from XML to CSV type: '
                   '\npython source_model_converter.py '
                   '--input-xml-file PATH_TO_SITE_MODEL_XML_FILE. '
                   '--output-csv-file PATH_TO_OUTPUT_CSV_FILE')
    parser = argparse.ArgumentParser(description=description,
                                     add_help=False,
                                     formatter_class=RawTextHelpFormatter)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    group_output = flags.add_mutually_exclusive_group()

    group_output.add_argument('--output-xml-file',
                              help='path to output xml file',
                              default=None,
                              required=False)
    group_output.add_argument('--output-csv-file',
                              help='path to output xml file',
                              default=None,
                              required=False)
    group_input = flags.add_mutually_exclusive_group()
    group_input.add_argument('--input-xml-file',
                       help='path to site model NRML file',
                       default=None)
    group_input.add_argument('--input-csv-file',
                       help='path to site model csv file',
                       default=None)

    return parser

if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()
    if args.input_csv_file:
        csv_to_xml(args.input_csv_file, args.output_xml_file)
    elif args.input_xml_file:
        xml_to_csv(args.input_xml_file, args.output_csv_file)
    else:
        parser.print_usage()
