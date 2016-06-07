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
Convert SES NRML file to rupture xmls.
'''
import abc
import os
import csv
import argparse
import numpy as np
from openquake.commonlib.node import striptag, LiteralNode
from openquake.commonlib import nrml
from openquake.hazardlib.geo.mesh import RectangularMesh
from openquake.hazardlib.geo import Point, Line, PlanarSurface


def parse_planar_surface(rupture):
    """

    """
    mag = LiteralNode("magnitude", text=rupture["magnitude"])
    rake = LiteralNode("rake", text=rupture["rake"])
    top_left = Point(float(rupture.planarSurface.topLeft["lon"]),
                     float(rupture.planarSurface.topLeft["lat"]),
                     float(rupture.planarSurface.topLeft["depth"]))
    top_right = Point(float(rupture.planarSurface.topRight["lon"]),
                     float(rupture.planarSurface.topRight["lat"]),
                     float(rupture.planarSurface.topRight["depth"]))
    bottom_left = Point(float(rupture.planarSurface.bottomLeft["lon"]),
                     float(rupture.planarSurface.bottomLeft["lat"]),
                     float(rupture.planarSurface.bottomLeft["depth"]))
    bottom_right = Point(float(rupture.planarSurface.bottomRight["lon"]),
                     float(rupture.planarSurface.bottomRight["lat"]),
                     float(rupture.planarSurface.bottomRight["depth"]))
    planar_surface = PlanarSurface(1.0,
        float(rupture["strike"]),
        float(rupture["dip"]),
        top_left,
        top_right,
        bottom_right,
        bottom_left)
    hypocentre = planar_surface.get_middle_point()
    hypo = LiteralNode("hypocenter", {
        "lon": str(hypocentre.longitude),
        "lat": str(hypocentre.latitude),
        "depth": str(hypocentre.depth)})
    surface = LiteralNode("planarSurface",
        {"strike": rupture["strike"], "dip": rupture["dip"]},
        nodes=rupture.planarSurface.nodes)
    return LiteralNode("singlePlaneRupture", nodes=[mag, rake, hypo, surface])

def _row_to_linestring(row, loc):
    """

    """
    ncol = row.shape[1]
    geom_str = ["%s %s %s" % (row[loc, i, 0], row[loc, i, 1], row[loc, i, 2])
                for i in range(ncol)]
    poslist_node = LiteralNode("gml:posList", text=geom_str)
    return LiteralNode("gml:LineString", nodes=[poslist_node])

def parse_mesh_surface(rupture):
    """
    Parse the rectangular mesh to a complex fault surface
    """
    nrow, ncol = int(rupture.mesh["rows"]), int(rupture.mesh["cols"])
    gridpoints = np.zeros([nrow, ncol, 3])
    for node in rupture.mesh:
        i, j = int(node["row"]), int(node["col"])
        gridpoints[i, j, 0] = float(node["lon"])
        gridpoints[i, j, 1] = float(node["lat"])
        gridpoints[i, j, 2] = float(node["depth"])
    edges = []
    for i in range(nrow):
        if i == 0:
            # Top edge
            edges.append(LiteralNode("faultTopEdge",
                nodes=[_row_to_linestring(gridpoints, i)]))
        elif i == (nrow - 1):
            edges.append(LiteralNode("faultBottomEdge",
                nodes=[_row_to_linestring(gridpoints, i)]))
        else:
            edges.append(LiteralNode("intermediateEdge",
                nodes=[_row_to_linestring(gridpoints, i)]))
    m_r, m_c = (nrow / 2), (ncol / 2)
    hypo = LiteralNode("hypocenter", {
        "lon": str(gridpoints[m_r, m_c, 0]),
        "lat": str(gridpoints[m_r, m_c, 1]),
        "depth": str(gridpoints[m_r, m_c, 2])})

    geom = LiteralNode("complexFaultGeometry", nodes=edges)
    mag = LiteralNode("magnitude", text=rupture["magnitude"])
    rake = LiteralNode("rake", text=rupture["rake"])
    return LiteralNode("complexFaultRupture", nodes=[mag, rake, hypo, geom])        


def event_set_to_rupture_xmls(input_ses, output_dir):
    """
    Parses the entire event set to a set of files
    """
    if os.path.exists(output_dir):
        raise IOError("Output directory %s already exists" % output_dir)
    else:
        os.mkdir(output_dir)
    nodeset = nrml.read(input_ses, chatty=False)
    for sesc in nodeset:
        sesc_dir = os.path.join(
            output_dir, 
            "smltp_{:s}".format(sesc["sourceModelTreePath"]))
        os.mkdir(sesc_dir)
        for i, ses in enumerate(sesc):
            ses_dir = os.path.join(sesc_dir, "ses_{:s}".format(str(ses["id"])))
            os.mkdir(ses_dir)
            for rupture in ses:
                print "Parsing event %s" % rupture["id"]
                if hasattr(rupture, "planarSurface"):
                    rupture_node = parse_planar_surface(rupture)
                elif hasattr(rupture, "mesh"):
                    rupture_node = parse_mesh_surface(rupture)
                rup_id = rupture["id"].replace("=", "_")
                filename = os.path.join(ses_dir,
                                        rup_id.replace("|", "_") + ".xml")
                with open(filename, "w") as f:
                    nrml.write([rupture_node], f, "%s")

def selected_event_set_to_rupture_xmls(input_ses, output_dir, selected_ids):
    """
    Parse only ruptures with the selected IDs to the output dir
    """
    if os.path.exists(output_dir):
        raise IOError("Output directory %s already exists" % output_dir)
    else:
        os.mkdir(output_dir)
    nodeset = nrml.read(input_ses, chatty=False)
    for sesc in nodeset:
        for ses in sesc:
            for rupture in ses:
                if rupture["id"] in selected_ids:
                    print "Parsing event %s" % rupture["id"]
                    if hasattr(rupture, "planarSurface"):
                        rupture_node = parse_planar_surface(rupture)
                    elif hasattr(rupture, "mesh"):
                        rupture_node = parse_mesh_surface(rupture)
                    rup_id = rupture["id"].replace("=", "_")
                    filename = os.path.join(output_dir,
                                            rup_id.replace("|", "_") + ".xml")
                    with open(filename, "w") as f:
                        nrml.write([rupture_node], f, "%s")


def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML stochastic event set file to tab delimited '
            ' .txt files. Inside the specified output directory, create a .txt '
            'file for each stochastic event set.'
            'To run just type: python eventset_converter.py '
            '--input-file=PATH_TO_GMF_NRML_FILE '
            '--output-dir=PATH_TO_OUTPUT_DIRECTORY', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
        help='path to ses NRML file (Required)',
        default=None,
        required=True)
    flags.add_argument('--output-dir',
        help='path to output directory (Required, raise an error if it already exists)',
        default=None,
        required=True)

    flags.add_argument('--selected-ids',
        help="Selected rupture IDs for parsing",
        default=[],
        nargs="+")

    return parser


if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        # create the output directory immediately. Raise an error if
        # it already exists

        if len(args.selected_ids):
            selected_event_set_to_rupture_xmls(args.input_file,
                                               args.output_dir,
                                               args.selected_ids)
        else:
            event_set_to_rupture_xmls(args.input_file, args.output_dir)
    else:
        parser.print_usage()
     
