#!/usr/bin/env python
# LICENSE
#
# Copyright (c) 2013, GEM Foundation, G. Weatherill, M. Pagani, D. Monelli.
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
'''
Convert NRML hazard curves file to .csv file.
'''
import os
import argparse
import numpy
from lxml import etree
import matplotlib.pyplot as plt

from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser


def _set_curves_matrix(hcm):
    """
    Store locations and poes in :class:`openquake.nrml.models.HazardCurveModel`
    in numpy array.
    """
    curves = []
    for loc, poes in hcm:
        row = [loc.x, loc.y]
        row.extend(poes)
        curves.append(row)

    return numpy.array(curves)

def _set_header(hcm):
    """
    Save metadata in :class:`openquake.nrml.models.HazardCurveModel`
    in a string to be used as header
    """
    header = ','.join(
        ['%s=%s' % (k,v) for k,v in hcm.metadata.items()
        if v is not None and k != 'imls']
    )
    header = '# ' + header
    header += '\nlon,lat,'+','.join([str(iml) for iml in hcm.metadata['imls']])

    return header

def plot_hazard_curve(filename_root, curves, hcm):
    """
    Exports the hazard curves to a set of pdf files
    """
    os.mkdir(filename_root)
    if ("PGA" in hcm.metadata["imt"]) or ("SA" in hcm.metadata["imt"]):
        imt_units = "g"
    else:
        imt_units = "cm/s"
    num_curves = numpy.shape(curves)[0]
    for iloc, row in enumerate(curves):
        print "Plotting curve %d of %d" % (iloc, num_curves)
        fig = plt.figure(figsize=(7, 5))
        fig.set_tight_layout(True)
        plt.loglog(hcm.metadata["imls"], row[2:], 'bo-', linewidth=2.0)
        plt.xlabel("%s (%s)" %(hcm.metadata["imt"], imt_units), fontsize=14)
        plt.ylabel("Probability of Being Exceeded in %s years" %
                   hcm.metadata["investigation_time"], fontsize=14)
        if row[0] < 0.0:
            long_ind = "W"
        else:
            long_ind = "E"
        if row[1] < 0.0:
            lat_ind = "S"
        else:
            lat_ind = "N"
        plt.title("Location: %12.6f %s, %12.6f %s" %(row[0], long_ind,
                                                     row[1], lat_ind))
        output_file = os.path.join(filename_root,
            "HazCurve_{:.5f}{:s}_{:.5f}{:s}.pdf".format(row[0], long_ind,
            row[1], lat_ind))
        plt.savefig(output_file, dpi=300, format="pdf", papertype="a4")

def save_hazard_curves_to_csv(nrml__hazard_curves_file, file_name_root,
    plot_curves=False):
    """
    Read hazard curves in `nrml__hazard_curves_file` and save to .csv file
    with root name `file_name_root`
    """
    output_file = '%s.csv' % file_name_root
    if os.path.isfile(output_file):
        raise ValueError('Output file already exists.'
                         ' Please specify different name or remove old file')

    hcm = HazardCurveXMLParser(nrml__hazard_curves_file).parse()

    curves = _set_curves_matrix(hcm)
    header = _set_header(hcm)
    f = open(output_file, 'w')
    f.write(header+'\n')
    numpy.savetxt(f, curves, fmt='%g', delimiter=',')
    f.close()
    if plot_curves:
        plot_hazard_curve(file_name_root, curves, hcm)



def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert Hazard Curves file from Nrml to .csv file.'
            'To run just type: python hazard_curve_converter.py ' 
            '--input-file=/PATH/TO/INPUT_FILE '
            '--output-file=/PATH/TO/OUTPUT_FILE', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
                        help='path to NRML hazard curve file (Required)',
                        default=None,
                        required=True)
    flags.add_argument('--output-file',
                        help='path to output file, without file extension '
                             ' (Optional, default is root of NRML file)',
                        default=None
                       )
    flags.add_argument('--plot-curves',
                       help="Plot the hazard curves to pdf (True) or not "
                       "(False) - may take time for many hazard curves",
                       default=False)
    return parser


if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        output_file = \
            os.path.splitext(parser.parse_args().input_file)[0] \
            if args.output_file is None else args.output_file

        save_hazard_curves_to_csv(args.input_file,
                                  output_file,
                                  args.plot_curves)
    else:
        parser.print_usage()
