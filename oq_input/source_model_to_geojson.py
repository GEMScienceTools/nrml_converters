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
Convert NRML source model file to geojson (http://geojson.org) file.
"""
import os
import argparse
import re
import geojson
from shapely import wkt

from openquake.nrmllib.hazard.parsers import SourceModelParser
from openquake.nrmllib.models import (
    PointSource,
    AreaSource,
    SimpleFaultSource,
    ComplexFaultSource,
    CharacteristicSource,
    IncrementalMFD,
    TGRMFD
)

from openquake.hazardlib.mfd import TruncatedGRMFD

AREA_GEO = re.compile('^POLYGON\(\(.+\)\)')
COMPLEX_GEO = re.compile('^LINESTRING\(.+\)\_LINESTRING\(.+\)')

def set_up_arg_parser():
    """
    Define command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML source model file into'
                    'geojson file'
    )

    parser.add_argument(
        '--nrml_source_model_file',
        help='NRML source model file to be converted',
        required=True
    )

    return parser

def _get_geometry(src):
    """
    Return source geometry
    """
    if isinstance(src, ComplexFaultSource):
        # for complex faults we neglect intermediate edges (sources that
        # have the same top and bottom edges will have anyhow the same
        # surface projection)
        return '_'.join([
            src.geometry.top_edge_wkt, src.geometry.bottom_edge_wkt
        ])

    elif isinstance(src, CharacteristicSource):
        raise ValueError(
            'Not yet able to extract geometry for CharacteristicSource'
        )
    else:
        return src.geometry.wkt

def _extract_sources(nrml_source_model_file):
    """
    Extract sources from NRML source model file.

    Sources are returned by typology, tectonic region type and 'geometry'
    coordinates by means of three nested dictionaries
    """
    srcm = SourceModelParser(nrml_source_model_file).parse()

    srcs = dict()

    for src in srcm:
        typology = src.__class__.__name__
        trt = src.trt
        geo = _get_geometry(src)

        if typology in srcs:
            if trt in srcs[typology]:
                if geo in srcs[typology][trt]:
                    srcs[typology][trt][geo].append(src)
                else:
                    srcs[typology][trt][geo] = [src]
            else:
                srcs[typology][trt] = {geo: [src]}
        else:
            srcs[typology] = {trt: {geo: [src]}}

    return srcs

def _get_geojson_geometry(geo):
    """
    Convert 'geo' string to geojson geometry
    """
    if AREA_GEO.match(geo):
        coords = wkt.loads(geo).exterior.coords
        coords = [[lon, lat] for lon, lat in coords]
        return geojson.Polygon([coords])

    elif COMPLEX_GEO.match(geo):
        coords = geo.split('_')
        top = [[lon, lat] for lon, lat, _ in wkt.loads(coords[0]).coords]
        bottom = [[lon, lat] for lon, lat, _ in wkt.loads(coords[1]).coords]

        # define boundary of complex fault
        top.extend(bottom[::-1])
        return geojson.Polygon([top])
    else:
        raise ValueError('Geometry %s not recognized' % geo)

def _get_mfds(srcs):
    """
    Return sources' magnitude frequency distributions.
    """
    mfds = {}
    tot_occur_rate = 0
    for src in srcs:
        mfd = {}

        if isinstance(src.mfd, IncrementalMFD):
            mfd['min_mag'] = float(src.mfd.min_mag)
            mfd['bin_width'] = float(src.mfd.bin_width)
            mfd['occur_rates'] = map(float, src.mfd.occur_rates)

        elif isinstance(src.mfd, TGRMFD):
            bin_width = 0.1
            mfd = TruncatedGRMFD(
                src.mfd.min_mag,
                src.mfd.max_mag,
                bin_width,
                src.mfd.a_val,
                src.mfd.b_val
            )

            mags_rates = mfd.get_annual_occurrence_rates()
            mags = [m for m, _ in mags_rates]
            occur_rates = [r for _, r in mags_rates]

            mfd['min_mag'] = mags[0]
            mfd['bin_width'] = bin_width
            mfd['occur_rates'] = occur_rates

        else:
            raise ValueError(
                'MFD %s not recognized' % src.mfd.__class__.__name__
            )

        tot_occur_rate += sum(mfd['occur_rates'])

        assert src.id not in mfds
        mfds[src.id] = mfd

    return mfds, tot_occur_rate

def _get_source_properties(srcs):
    """
    Extract source properties depending on source typology.
    """
    props = {}
    for src in srcs:

        data = {}
        # this is valid also for Area Sources
        if isinstance(src, PointSource):
            data['upper_seismo_depth'] = float(src.geometry.upper_seismo_depth)
            data['lower_seismo_depth'] = float(src.geometry.lower_seismo_depth)
            data['mag_scale_rel'] = src.mag_scale_rel
            data['rupt_aspect_ratio'] = float(src.rupt_aspect_ratio)

            nodal_plane_dist = [
                map(float, [np.probability, np.strike, np.dip, np.rake]) \
                for np in src.nodal_plane_dist
            ]
            data['nodal_plane_dist'] = nodal_plane_dist

            hypo_depth_dist = [
                map(float, [hd.probability, hd.depth]) \
                for hd in src.hypo_depth_dist
            ]
            data['hypo_depth_dist'] = hypo_depth_dist

        elif isinstance(src, ComplexFaultSource):
            data['mag_scale_rel'] = src.mag_scale_rel
            data['rupt_aspect_ratio'] = float(src.rupt_aspect_ratio)
            data['rake'] = float(src.rake)

        elif isinstance(src, SimpleFaultSource):
            data['upper_seismo_depth'] = float(src.geometry.upper_seismo_depth)
            data['lower_seismo_depth'] = float(src.geometry.lower_seismo_depth)
            data['dip'] = float(src.geometry.dip)
            data['mag_scale_rel'] = src.mag_scale_rel
            data['rupt_aspect_ratio'] = float(src.rupt_aspect_ratio)
            data['rake'] = float(src.rake)

        else:
            raise ValueError(
                'Source typology %s not recognized' % src.__class__.__name__
            )

        props[src.id] = data

    return props

def _geojson(srcs, root):
    """
    Serialize srcs data as returned from method '_extract_sources' to
    geojson files (one per tectonic region type).
    """
    for typology, trt_geo_srcs in srcs.items():
        for trt, geo_srcs in trt_geo_srcs.items():
            features = []

            for geo, srcs in geo_srcs.items():
                geometry = _get_geojson_geometry(geo)
                src_ids = [src.id for src in srcs]
                src_names = [src.name for src in srcs]
                mfds, tot_occur_rate = _get_mfds(srcs)
                src_properties = _get_source_properties(srcs)

                if geometry:
                    feature = geojson.Feature(
                        geometry=geometry,
                        properties={
                            'typology': typology,
                            'trt': trt,
                            'tot_occur_rate': tot_occur_rate,
                            'src_ids': src_ids,
                            'src_names': src_names,
                            'mfds': mfds,
                            'src_properties': src_properties
                        }
                    )
                    features.append(feature)

            features = geojson.FeatureCollection(features)
            dump = geojson.dumps(features, sort_keys=True)
            f = open(
                '%s_%s_%s.geojson' % \
                (root, typology, trt.replace(' ', '_')), 'w'
            )
            f.write(dump)
            f.close()

if __name__ == '__main__':

    parser = set_up_arg_parser()
    args = parser.parse_args()

    root, _ = os.path.splitext(args.nrml_source_model_file)

    srcs = _extract_sources(args.nrml_source_model_file)
    _geojson(srcs, root)
