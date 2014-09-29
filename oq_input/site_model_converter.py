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
import argparse
import openquake.nrmllib
from openquake.nrmllib.hazard.parsers import SiteModelParser
from lxml import etree
from argparse import RawTextHelpFormatter
import numpy as np
from shapely import wkt

def csv_to_xml(input_csv, output_xml, geo_precision="%12.6f",
        site_precision="%8.2f"):
    """
    Converts the csv file to the xml format
    """
    data = np.genfromtxt(input_csv, delimiter=",", skip_header=1)
    nsites = np.shape(data)[0]
    with open(output_xml, "w") as fh:
        root = etree.Element('nrml',
                              nsmap=openquake.nrmllib.SERIALIZE_NS_MAP)
        site_model = etree.SubElement(root, "siteModel")
        for iloc in range(0, nsites):
            site_value = etree.SubElement(site_model, "site")
            site_value.set("lon", geo_precision % data[iloc, 0])
            site_value.set("lat", geo_precision % data[iloc, 1])
            site_value.set("vs30", site_precision % data[iloc, 2])
            if data[iloc, 3]:
                site_value.set("vs30Type", "measured")
            else:
                site_value.set("vs30Type", "inferred")
            site_value.set("z1pt0", site_precision % data[iloc, 4])
            site_value.set("z2pt5", site_precision % data[iloc, 5])
        fh.write(etree.tostring(root, pretty_print=True,
                                xml_declaration=True, encoding='UTF-8'))


def xml_to_csv(input_xml, output_csv, geo_precision="%12.6f",
        site_precision="%8.2f"):
    """
    Converts the xml to csv format
    """
    parser = SiteModelParser(input_xml)
    sites = list(parser.parse())
    f = open(output_csv, "w")

    target_string = ",".join([geo_precision, geo_precision, site_precision,
                              "%s", site_precision, site_precision])
    for site in sites:
        locn = wkt.loads(site.wkt)
        print >> f, target_string % (locn.x, locn.y, site.vs30, site.vs30_type,
                                     site.z1pt0, site.z2pt5)
    f.close()

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
