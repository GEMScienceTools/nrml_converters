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
Convert SES NRML file to .txt files.
'''
import abc
import os
import csv
import argparse
import numpy
from lxml import etree
from openquake.commonlib.node import striptag
from openquake.commonlib.nrml import read
from openquake.hazardlib.geo.mesh import RectangularMesh

NRML='{http://openquake.org/xmlns/nrml/0.5}'

PLANAR_TAGS = ["{:s}".format(tag)
               for tag in ["topLeft", "topRight", "bottomLeft", "bottomRight"]]


def parse_sesc_file(file_name):
    """
    Parse NRML 0.4 SES collection file.
    """
    
    element = read(file_name, chatty=False)[0]
    #sesc = parse_ses_collection(element)
    return parse_ses_collection(element)
#    parse_args = dict(source=file_name)
#  
#    for _, element in etree.iterparse(**parse_args):
#        if element.tag == '%sstochasticEventSetCollection' % NRML:
#            sesc = parse_ses_collection(element)
#
#    return sesc

def parse_ses_collection(element):
    """
    Parse NRML 0.4 'stochasticEventSetCollection' and return
    class StochasticEventSetCollection.
    """
    #smtp = element.attrib['sourceModelTreePath']

    sess = []
    num_ses = 0
    for node in element.nodes:
        sess.append(parse_ses(node))

    print 'number of stochastic event sets: %s' % len(sess)

    #return StochasticEventSetCollection(smtp, sess)
    return StochasticEventSetCollection(sess)


def parse_ses(element):
    """
    Parse NRML 0.4 'stochasticEventSet' and return
    class StochasticEventSet.
    """
    ID = element.attrib['id']
    time_span = element.attrib['investigationTime']

    print 'Processing ses %s' % ID

    rups = []
    for node in element.nodes:
        rups.append(parse_rup(node))

    return StochasticEventSet(ID, time_span, rups)

def parse_rup(element):
    """
    Parse NRML 0.4 'rupture' and return class Rupture.
    """
    ID = element.attrib['id']
    magnitude = float(element.attrib['magnitude'])
    strike = float(element.attrib['strike'])
    dip = float(element.attrib['dip'])
    rake = float(element.attrib['rake'])
    tectonic_region = element.attrib['tectonicRegion']

    print 'Processing rupture %s' % ID
    for node in element.nodes:
        if "planarSurface" in node.tag:
            surf = parse_planar_surf(node)
        else:
            surf = parse_mesh(node)
    return Rupture(ID, magnitude, strike, dip, rake, tectonic_region, surf)

#    
#    
#    planar_surfs = element.findall('%splanarSurface' % NRML)
#    mesh = element.find('%smesh' % NRML)
#
#    if planar_surfs:
#        if len(planar_surfs) == 1:
#            surf = parse_planar_surf(planar_surfs[0])
#        else:
#            surf = []
#            for e in planar_surfs:
#                surf.append(parse_planar_surf(e))
#            surf = MultiPlanarSurface(surf)
#    else:
#        assert(mesh is not None)
#        surf = parse_mesh(mesh)
#
#    return Rupture(ID, magnitude, strike, dip, rake, tectonic_region, surf)

def parse_planar_surf(element):
    """
    Parse NRML 0.4 'planarSurface' and return class PlanarSurface.
    """
    tag_set = [striptag(node.tag) for node in element.nodes]
    (top_left, top_right, bottom_left, bottom_right) = tuple(
        [element.nodes[tag_set.index(tag)] for tag in PLANAR_TAGS])
#    for node in element.nodes:
#        if "topLeft" in node.tag:
#            top_left = copy(node)
#        elif "topRight" in node.tag:
#            top_right = copy(node)
#        elif "bottomRight" in node.tag:
#          
#
#    top_left = element.find('%stopLeft' % NRML)
#    top_right = element.find('%stopRight' % NRML)
#    bottom_left = element.find('%sbottomLeft' % NRML)
#    bottom_right = element.find('%sbottomRight' % NRML)

    corners_lons = numpy.array(
        [float(top_left.attrib['lon']), float(top_right.attrib['lon']),
        float(bottom_right.attrib['lon']), float(bottom_left.attrib['lon'])]
    )

    corners_lats = numpy.array(
        [float(top_left.attrib['lat']), float(top_right.attrib['lat']),
        float(bottom_right.attrib['lat']), float(bottom_left.attrib['lat'])]
    )

    corners_depths = numpy.array(
        [float(top_left.attrib['depth']), float(top_right.attrib['depth']),
        float(bottom_right.attrib['depth']), float(bottom_left.attrib['depth'])]
    )

    return PlanarSurface(corners_lons, corners_lats, corners_depths)


def parse_mesh(element):
    """
    Parse NRML 0.4 'mesh' and return class MeshSurface.
    """
    nrows = int(element.attrib['rows'])
    ncols = int(element.attrib['cols'])

    lons = numpy.ndarray((nrows, ncols))
    lats = numpy.ndarray((nrows, ncols))
    depths = numpy.ndarray((nrows, ncols))

    for e in element.nodes:
        row = int(e.attrib['row'])
        col = int(e.attrib['col'])
        lons[row, col] = float(e.attrib['lon'])
        lats[row, col] = float(e.attrib['lat'])
        depths[row, col] = float(e.attrib['depth'])

    return MeshSurface(lons, lats, depths)


#def parse_mesh(element):
#    """
#    Parse NRML 0.4 'mesh' and return class MeshSurface.
#    """
#    nrows = int(element.attrib['rows'])
#    ncols = int(element.attrib['cols'])
#
#    lons = numpy.ndarray((nrows, ncols))
#    lats = numpy.ndarray((nrows, ncols))
#    depths = numpy.ndarray((nrows, ncols))
#
#    for e in element.iterchildren():
#        row = int(e.attrib['row'])
#        col = int(e.attrib['col'])
#        lons[row, col] = float(e.attrib['lon'])
#        lats[row, col] = float(e.attrib['lat'])
#        depths[row, col] = float(e.attrib['depth'])
#
#    return MeshSurface(lons, lats, depths)


class BaseSurface(object):
    """
    Class representing base surface.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_surface_boundaries(self):
        """
        Return surface boundaries.
        """

    @abc.abstractmethod
    def get_middle_point(self):
        """
        Return surface middle point
        """

class PlanarSurface(BaseSurface):
    """
    Class representing a planar surface in terms of its four corner
    coordinates.
    """
    def __init__(self, corners_lons, corners_lats, corners_depths):
        self.corners_lons = corners_lons
        self.corners_lats = corners_lats
        self.corners_depths = corners_depths

    def get_surface_boundaries(self):
        return [self.corners_lons.take([0, 1, 2, 3, 0])], \
               [self.corners_lats.take([0, 1, 2, 3, 0])]

    def get_middle_point(self):
        lons = numpy.array([self.corners_lons.take([0, 1]),
                            self.corners_lons.take([3, 2])])
        lats = numpy.array([self.corners_lats.take([0, 1]),
                            self.corners_lats.take([3, 2])])
        depths = numpy.array([self.corners_depths.take([0, 1]),
                              self.corners_depths.take([3, 2])])                                  
        mesh = RectangularMesh(lons, lats, depths)

        middle_point = mesh.get_middle_point()

        return (middle_point.longitude, middle_point.latitude, middle_point.depth)

class MultiPlanarSurface(BaseSurface):
    """
    Class representing surface as collection of planar surfaces.
    """
    def __init__(self, surfaces):
        self.surfaces = surfaces

    def get_surface_boundaries(self):
        lons = []
        lats = []
        for surf in self.surfaces:
            assert(isinstance(surf, PlanarSurface))
            lons_surf, lats_surf = surf.get_surface_boundaries()
            lons.append(lons_surf[0])
            lats.append(lats_surf[0])

        return lons, lats

    def get_middle_point(self):
        surf = self.surfaces[0]
        assert(isinstance(surf, PlanarSurface))
        return surf.get_middle_point()

class MeshSurface(BaseSurface):
    """
    Class representing surface as mesh of points.
    """
    def __init__(self, lons, lats, depths):
        self.lons = lons
        self.lats = lats
        self.depths = depths

    def get_surface_boundaries(self):
        lons = numpy.concatenate((self.lons[0,:],
                                  self.lons[1:,-1],
                                  self.lons[-1,:-1][::-1],
                                  self.lons[:-1,0][::-1]))
        lats = numpy.concatenate((self.lats[0,:],
                                  self.lats[1:,-1],
                                  self.lats[-1,:-1][::-1],
                                  self.lats[:-1,0][::-1]))
        return [lons], [lats]

    def get_middle_point(self):
        mesh = RectangularMesh(self.lons, self.lats, self.depths)
        middle_point = mesh.get_middle_point()
        return (middle_point.longitude, middle_point.latitude, middle_point.depth)

class Rupture(object):
    """
    Class representing a single rupture, in terms of
    ID, magnitude, strike, dip, rake, tectonic region type,
    and a surface.
    """
    def __init__(self, ID, magnitude, strike, dip, rake, tect_reg, surf):
        self.ID = ID
        self.magnitude = magnitude
        self.strike = strike
        self.dip = dip
        self.rake = rake
        self.tect_reg = tect_reg
        self.surf = surf

class StochasticEventSet(object):
    """
    Class representing a single SES, identified by an ID
    and associated to a time span, and constructed from
    a list of ruptures.
    """
    def __init__(self, ID, time_span, rups):
        self.ID = ID
        self.time_span = time_span
        self.rups = rups

class StochasticEventSetCollection(object):
    """
    Class representing collection of SESs associated to
    given source model and GSIM logic tree paths.
    """
    def __init__(self, sess):
#    def __init__(self, smtp, sess):
        data = []
        for ses in sess:
            for rup in ses.rups:
                lon, lat, depth = rup.surf.get_middle_point()
                multi_lons, multi_lats = rup.surf.get_surface_boundaries()
                boundary = 'MULTIPOLYGON(%s)' % \
                    ','.join(
                            '((%s))' % ','.join('%s %s' % \
                                (lon, lat) for lon, lat in zip(lons, lats)) \
                                for lons, lats in zip(multi_lons, multi_lats)
                        )
                data.append([ses.ID, rup.ID, rup.magnitude,
                             lon, lat, depth, rup.tect_reg, rup.strike,
                             rup.dip, rup.rake, boundary])
#                data.append([smtp, ses.ID, rup.ID, rup.magnitude,
#                             lon, lat, depth, rup.tect_reg, rup.strike,
#                             rup.dip, rup.rake, boundary])

        self.data = numpy.array(data, dtype=object)


def save_sess_to_txt(sesc, output_dir):
    """
    Save stochastic event sets collection to .txt files.
    """
    ses_ids = numpy.unique(sesc.data[:, 1])
    for ID in ses_ids:
        idx = sesc.data[:, 1] == ID
        fname = '%s/ses_%s.txt' % (output_dir, ID)
        header = 'id\tmag\tcentroid_lon\tcentroid_lat\tcentroid_depth\ttrt\tstrike\tdip\trake\tboundary'
        f = open(fname, 'w')
        f.write(header+'\n')
        numpy.savetxt(f, sesc.data[idx, 2 :],
        fmt='t%2.1f\t%5.2f\t%5.2f\t%5.2f\t%s\t%5.2f\t%5.2f\t%5.2f\t%s')
#            fmt='%s\t%2.1f\t%5.2f\t%5.2f\t%5.2f\t%s\t%5.2f\t%5.2f\t%5.2f\t%s')
        f.close()


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

    return parser


if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        # create the output directory immediately. Raise an error if
        # it already exists
        os.makedirs(args.output_dir)

        sesc = parse_sesc_file(args.input_file)
        save_sess_to_txt(sesc, args.output_dir)
    else:
        parser.print_usage()
