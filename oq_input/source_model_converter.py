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
Convert NRML source model file to ESRI shapefile (and vice versa).
"""
import os
import subprocess
import sys
import numpy
import argparse
import shapefile
from copy import copy
from shapely import wkt
from argparse import RawTextHelpFormatter
from collections import OrderedDict

from openquake.hazardlib import scalerel
from openquake.hazardlib.source.rupture import Rupture as HazardlibRupture
from openquake.hazardlib import geo, mfd, pmf, source
from openquake.hazardlib.geo.surface import (ComplexFaultSurface,
                                             SimpleFaultSurface,
                                             PlanarSurface,
                                             MultiSurface)
from openquake.hazardlib.tom import PoissonTOM

#from openquake.commonlib.source import parse_source_model
from openquake.commonlib import nrml
#from openquake.commonlib.nrml import nodefactory
from openquake.commonlib.node import read_nodes, LiteralNode, striptag
from openquake.commonlib.sourceconverter import SourceConverter
#from openquake.commonlib.sourcewriter import write_source_model

from input_utils import fix_source_node, AKI_RICH_ERR_MSG, WRONG_ORDER_ERR_MSG

# maximum field size allowed by shapefile
FIELD_SIZE = 255
# maximum number of occurrence rates that can be stored for incremental MFD
MAX_RATES = 50
# maximum number of nodal planes
MAX_NODAL_PLANES = 20
# maximum number of hypocentral depths
MAX_HYPO_DEPTHS = 20
# maximum number of planes (for sources described by multi-surface)
MAX_PLANES = 10

# each triplet contains nrmllib parameter name, shapefile field name and
# data type
BASE_PARAMS = [
    ('id', 'id', 'c'), ('name', 'name', 'c'),
    ('tectonicRegion', 'trt', 'c'),
    ('magScaleRel', 'msr', 'c'),
    ('ruptAspectRatio', 'rar', 'f'),
    ('rake', 'rake', 'f')
]
GEOMETRY_PARAMS = [
    ('upperSeismoDepth', 'usd', 'f'),
    ('lowerSeismoDepth', 'lsd', 'f'),
    ('dip', 'dip', 'f')
]
MFD_PARAMS = [
    ('minMag', 'min_mag', 'f'), ('maxMag', 'max_mag', 'f'),
    ('aValue', 'a_val', 'f'), ('bValue', 'b_val', 'f'),
    ('binWidth', 'bin_width', 'f')
]


# shapefile specific fields
RATE_PARAMS = [('rate%s' % (i+1), 'f') for i in range(MAX_RATES)]
STRIKE_PARAMS = [('strike%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
DIP_PARAMS = [('dip%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
RAKE_PARAMS = [('rake%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
NPW_PARAMS = [('npweight%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
HDEPTH_PARAMS = [('hd%s' % (i+1), 'f') for i in range(MAX_HYPO_DEPTHS)]
HDW_PARAMS = [('hdweight%s' % (i+1), 'f') for i in range(MAX_HYPO_DEPTHS)]
PLANES_STRIKES_PARAM = [('pstrike%s' % (i+1), 'f') for i in range(MAX_PLANES)]
PLANES_DIPS_PARAM = [('pdip%s' % (i+1), 'f') for i in range(MAX_PLANES)]

def get_taglist(node):
    """
    Return a list of tags (with NRML namespace removed) representing the
    order of the nodes within a node
    """
    return [striptag(subnode.tag) for subnode in node]
    #return [re.sub(r'\{[^}]*\}', "", copy(subnode.tag))
    #        for subnode in node.nodes]

def register_fields(w):
    """
    Register shapefile fields.
    """
    PARAMS_LIST = [BASE_PARAMS, GEOMETRY_PARAMS, MFD_PARAMS]
    for PARAMS in PARAMS_LIST:
        for _, param, dtype in PARAMS:
            w.field(param, fieldType=dtype, size=FIELD_SIZE)

    PARAMS_LIST = [
        RATE_PARAMS, STRIKE_PARAMS, DIP_PARAMS, RAKE_PARAMS, NPW_PARAMS,
        HDEPTH_PARAMS, HDW_PARAMS, PLANES_STRIKES_PARAM, PLANES_DIPS_PARAM
    ]
    for PARAMS in PARAMS_LIST:
        for param, dtype in PARAMS:
            w.field(param, fieldType=dtype, size=FIELD_SIZE)

    # source typology
    w.field('sourcetype', 'C')

def expand_src_param(values, shp_params):
    """
    Expand hazardlib source attribute (defined through list of values)
    into dictionary of shapefile parameters.
    """
    if values is None:
        return dict([(key, None) for key, _ in shp_params])
    else:
        num_values = len(values)
        return OrderedDict(
            [(key, float(values[i]) if i < num_values else None)
                for i, (key, _) in enumerate(shp_params)]
        )

def extract_source_params(src):
    """
    Extract params from source object.
    """
    tags = get_taglist(src)
    data = []
    for key, param, vtype in BASE_PARAMS:
        if key in src.attrib:
            if vtype == "c":
                data.append((param, src.attrib[key]))
            elif vtype == "f":
                data.append((param, float(src.attrib[key])))
            else:
                data.append((param, None))
        elif key in tags:
            if vtype == "c":
                data.append((param, src.nodes[tags.index(key)].text))
            elif vtype == "f":
                data.append((param, float(src.nodes[tags.index(key)].text)))
            else:
                data.append((param, None))
        else:
            data.append((param, None))
    return OrderedDict(data)

def parse_area_geometry(node):
    """

    """
    geometry = {}
    for subnode in node:
        if "upperSeismoDepth" in subnode.tag:
            geometry["upperSeismoDepth"] = subnode.text
        elif "lowerSeismoDepth" in subnode.tag:
            geometry["lowerSeismoDepth"] = subnode.text
        elif "Polygon" in subnode.tag:
            crds = subnode.nodes[0].nodes[0].nodes[0].text
            geometry["polygon"] = numpy.array([[crds[i], crds[i + 1]]
                                              for i in range(0, len(crds), 2)])
        else:
            pass
    return geometry

def parse_point_geometry(node):
    """

    """
    geometry = {}
    for subnode in node:
        if "upperSeismoDepth" in subnode.tag:
            geometry["upperSeismoDepth"] = subnode.text
        elif "lowerSeismoDepth" in subnode.tag:
            geometry["lowerSeismoDepth"] = subnode.text
        elif "Point" in subnode.tag:
            lon, lat = subnode.nodes[0].text
            geometry["point"] = numpy.array([lon, lat])
        else:
            pass
    return geometry
    

def parse_simple_fault_geometry(node):
    """
    Parses a simple fault geometry node returning both the attributes and
    parameters in a dictionary
    """
    assert "simpleFaultGeometry" in node.tag
    # Get general attributes
    geometry = {}
    for subnode in node:
        if "dip" in subnode.tag:
            geometry["dip"] = subnode.text
        elif "upperSeismoDepth" in subnode.tag:
            geometry["upperSeismoDepth"] = subnode.text
        elif "lowerSeismoDepth" in subnode.tag:
            geometry["lowerSeismoDepth"] = subnode.text
        elif "LineString" in subnode.tag:
            crds = subnode.nodes[0].text
            geometry["trace"] = numpy.array([[crds[i], crds[i + 1]]
                                            for i in range(0, len(crds), 2)])
        else:
            pass
    return geometry

def parse_complex_fault_geometry(node):
    """
    Parses a complex fault geometry node returning both the attributes and
    parameters in a dictionary
    """
    assert "complexFaultGeometry" in node.tag
    # Get general attributes
    geometry = {"intermediateEdges": []}
    for subnode in node:
        crds = subnode.nodes[0].nodes[0].text
        if "faultTopEdge" in subnode.tag:
            geometry["faultTopEdge"] = numpy.array(
                [[crds[i], crds[i + 1], crds[i + 2]]
                for i in range(0, len(crds), 3)])
            geometry["upperSeismoDepth"] = numpy.min(
                geometry["faultTopEdge"][:, 2])
        elif "faultBottomEdge" in subnode.tag:

            geometry["faultBottomEdge"] = numpy.array(
                [[crds[i], crds[i + 1], crds[i + 2]]
                for i in range(0, len(crds), 3)])
            geometry["lowerSeismoDepth"] = numpy.max(
                geometry["faultBottomEdge"][:, 2])
        elif "intermediateEdge" in subnode.tag:
            geometry["intermediateEdges"].append(
                numpy.array([[crds[i], crds[i + 1], crds[i + 2]]
                          for i in range(0, len(crds), 3)]))
        else:
            pass
    geometry["dip"] = None
    return geometry

def parse_planar_fault_geometry(node):
    """
    Parses a planar fault geometry node returning both the attributes and
    parameters in a dictionary
    """
    assert "planarSurface" in node.tag
    geometry = {"strike": node.attrib["strike"],
                "dip": node.attrib["dip"]}
    upper_depth = numpy.inf
    lower_depth = 0.0
    tags = get_taglist(node)
    corner_points = []
    for locn in ["topLeft", "topRight", "bottomRight", "bottomLeft"]:
        plane = node.nodes[tags.index(locn)]
        upper_depth = plane["depth"] if plane["depth"] < upper_depth else\
            upper_depth
        lower_depth = plane["depth"] if plane["depth"] > lower_depth else\
            lower_depth
        corner_points.append([plane["lon"], plane["lat"], plane["depth"]])
    geometry["upperSeismoDepth"] = upper_depth
    geometry["lowerSeismoDepth"] = lower_depth
    geometry["corners"] = numpy.array(corner_points)
    return geometry

def get_attrs_dict(attrs, params):
    """

    """
    data = []
    for key, param, _ in params:
        if key in attrs and attrs[key] is not None:
            data.append((param, attrs[key]))
        else:
            data.append((param, None))
    return OrderedDict(data)
    #return OrderedDict(
    #    [(param, attrs[key] if key in attrs and attrs[key] is not None else None)
    #     for key, param, _ in params])

def extract_geometry_params(src):
    """

    """
    tags = get_taglist(src)
    if "areaSource" in src.tag:
        geom = src.nodes[tags.index("areaGeometry")]
        area_attrs = parse_area_geometry(geom)
        return get_attrs_dict(area_attrs, GEOMETRY_PARAMS)
        
        #OrderedDict([
        #    (param, area_attrs[key] if key in area_attrs else None)
        #    for key, param, _ in GEOMETRY_PARAMS])
    elif "pointSource" in src.tag:
        geom = src.nodes[tags.index("pointGeometry")]
        point_attrs = parse_point_geometry(geom)
        return get_attrs_dict(point_attrs, GEOMETRY_PARAMS)
        #OrderedDict([
        #    (param, point_attrs[key] if key in point_attrs else None)
        #    for key, param, _ in GEOMETRY_PARAMS])
    elif "simpleFaultSource" in src.tag:
        simple_attrs = parse_simple_fault_geometry(
            src.nodes[tags.index("simpleFaultGeometry")])
        return get_attrs_dict(simple_attrs, GEOMETRY_PARAMS)
        #return OrderedDict([(param, simple_attrs[key])
        #                     for key, param, _ in GEOMETRY_PARAMS])
    elif "complexFaultSource" in src.tag:
        complex_attrs = parse_complex_fault_geometry(
            src.nodes[tags.index("complexFaultGeometry")])
        return get_attrs_dict(complex_attrs, GEOMETRY_PARAMS)
        #return OrderedDict([(param, complex_attrs[key])
        #                     for key, param, _ in GEOMETRY_PARAMS])
    elif "characteristicFaultSource" in src.tag:
        surface_node = src.nodes[tags.index("surface")]
        upper_depth = numpy.inf
        lower_depth = 0.0
        dip = 0.0
        counter = 0.0
        for surface in surface_node:
            if "simpleFaultGeometry" in surface.tag:
                surf_attrs = parse_simple_fault_geometry(surface)
            elif "complexFaultGeometry" in surface.tag:
                surf_attrs = parse_complex_fault_geometry(surface)
            elif "planarSurface" in surface.tag:
                surf_attrs = parse_planar_fault_geometry(surface)
            else:
                continue
            upper_depth = surf_attrs["upperSeismoDepth"] \
                if surf_attrs["upperSeismoDepth"] < upper_depth else\
                upper_depth
            lower_depth = surf_attrs["lowerSeismoDepth"] \
                if surf_attrs["lowerSeismoDepth"] > lower_depth else\
                lower_depth
            if surf_attrs["dip"]:
                dip += surf_attrs["dip"]
                counter += 1.0
            else:
                dip = None

        if dip and counter:
           dip /= counter
            
        return OrderedDict([("usd", upper_depth),
                            ("lsd", lower_depth),
                            ("dip", dip)])
    else:
        return OrderedDict()
            

def extract_mfd_params(src):
    """
    Extracts the MFD parameters from an object
    """
    tags = get_taglist(src)
    if "incrementalMFD" in tags:
        mfd_node = src.nodes[tags.index("incrementalMFD")]
    elif "truncGutenbergRichterMFD":
        mfd_node = src.nodes[tags.index("truncGutenbergRichterMFD")]
    else:
        raise ValueError("Source %s contains no supported MFD type!" % src.tag)
    data = []
    rates = []
    for key, param, vtype in MFD_PARAMS:
        if key in mfd_node.attrib and mfd_node.attrib[key] is not None:
            data.append((param, mfd_node.attrib[key]))
        else:
            data.append((param, None))
    if "incrementalMFD" in mfd_node.tag:
        # Extract Rates
        assert "occurRates" in mfd_node.nodes[0].tag
        rates = mfd_node.nodes[0].text

        n_r = len(rates)
        if n_r > MAX_RATES:
            raise ValueError("Number of rates in source %s too large "
                             "to be placed into shapefile" % src.tag)
        rate_dict = OrderedDict([(key, rates[i] if i < n_r else None)
                                  for i, (key, _) in enumerate(RATE_PARAMS)
                                  ])
    else:
        rate_dict = OrderedDict([(key, None)
                                 for i, (key, _) in enumerate(RATE_PARAMS)])

            
    return OrderedDict(data), rate_dict


def extract_source_nodal_planes(src):
    """

    """
    if not "pointSource" in src.tag and not "areaSource" in src.tag:
        strikes = dict([(key, None) for key, _ in STRIKE_PARAMS])
        dips = dict([(key, None) for key, _ in DIP_PARAMS])
        rakes = dict([(key, None) for key, _ in RAKE_PARAMS])
        np_weights = dict([(key, None) for key, _ in NPW_PARAMS])
        return strikes, dips, rakes, np_weights
    tags = get_taglist(src)
    npd_nodeset = src.nodes[tags.index("nodalPlaneDist")]
    if len(npd_nodeset) > MAX_NODAL_PLANES:
        raise ValueError("Number of nodal planes %s exceeds stated maximum "
                         "of %s" % (str(len(npd_nodeset)),
                                    str(MAX_NODAL_PLANES)))
    if len(npd_nodeset):
        strikes = []
        dips = []
        rakes = []
        np_weights = []
        for npd_node in npd_nodeset:
            strikes.append(float(npd_node.attrib["strike"]))
            dips.append(float(npd_node.attrib["dip"]))
            rakes.append(float(npd_node.attrib["rake"]))
            np_weights.append(float(npd_node.attrib["probability"]))
        strikes = expand_src_param(strikes, STRIKE_PARAMS)
        dips = expand_src_param(dips, DIP_PARAMS)
        rakes = expand_src_param(rakes, RAKE_PARAMS)
        np_weights = expand_src_param(np_weights, NPW_PARAMS)
    else:
        strikes = dict([(key, None) for key, _ in STRIKE_PARAMS])
        dips = dict([(key, None) for key, _ in DIP_PARAMS])
        rakes = dict([(key, None) for key, _ in RAKE_PARAMS])
        np_weights = dict([(key, None) for key, _ in NPW_PARAMS])
    return strikes, dips, rakes, np_weights
   
def extract_source_hypocentral_depths(src):
    """
    Extract source hypocentral depths.
    """
    if not "pointSource" in src.tag and not "areaSource" in src.tag:
        hds = dict([(key, None) for key, _ in HDEPTH_PARAMS])
        hdsw = dict([(key, None) for key, _ in HDW_PARAMS])
        return hds, hdsw
    
    tags = get_taglist(src)
    hdd_nodeset = src.nodes[tags.index("hypoDepthDist")]
    if len(hdd_nodeset) > MAX_HYPO_DEPTHS:
        raise ValueError("Number of hypocentral depths %s exceeds stated maximum "
                         "of %s" % (str(len(hdd_nodeset)),
                                    str(MAX_HYPO_DEPTHS)))
    if len(hdd_nodeset):
        hds = []
        hdws = []
        for hdd_node in hdd_nodeset:
            hds.append(float(hdd_node.attrib["depth"]))
            hdws.append(float(hdd_node.attrib["probability"]))

        hds = expand_src_param(hds, HDEPTH_PARAMS)
        hdsw = expand_src_param(hdws, HDW_PARAMS)
    else:
        hds = dict([(key, None) for key, _ in HDEPTH_PARAMS])
        hdsw = dict([(key, None) for key, _ in HDW_PARAMS])

    return hds, hdsw


def extract_source_planes_strikes_dips(src):
    """
    Extract strike and dip angles for source defined by multiple planes.
    """
    if not "characteristicFaultSource" in src.tag: 
        strikes = dict([(key, None) for key, _ in PLANES_STRIKES_PARAM])
        dips = dict([(key, None) for key, _ in PLANES_DIPS_PARAM])
        return strikes, dips
    tags = get_taglist(src)
    surface_set = src.nodes[tags.index("surface")]
    strikes = []
    dips = []
    num_planes = 0
    for surface in surface_set:
        if "planarSurface" in surface.tag:
            strikes.append(float(surface.attrib["strike"]))
            dips.append(float(surface.attrib["dip"]))
            num_planes += 1
    if num_planes > MAX_PLANES:
        raise ValueError("Number of planes in sourcs %s exceededs maximum "
                         "of %s" % (str(num_planes), str(MAX_PLANES)))
    if num_planes:
        strikes = expand_src_param(strikes, PLANES_STRIKES_PARAM)
        dips = expand_src_param(dips, PLANES_DIPS_PARAM)
    else:
        strikes = dict([(key, None) for key, _ in PLANES_STRIKES_PARAM])
        dips = dict([(key, None) for key, _ in PLANES_DIPS_PARAM])

    return strikes, dips

def set_params(w, src):
    """
    Set source parameters.
    """
    params = extract_source_params(src)
    # this is done because for characteristic sources geometry is in
    # 'surface' attribute
    params.update(extract_geometry_params(src))
        
    #params.update(extract_source_params(
    #    src.geometry if getattr(src, 'geometry', None) else src.surface,
    #    GEOMETRY_PARAMS)
    #)
    mfd_pars, rate_pars = extract_mfd_params(src)
    params.update(mfd_pars)
    params.update(rate_pars)

    strikes, dips, rakes, np_weights = extract_source_nodal_planes(src)
    params.update(strikes)
    params.update(dips)
    params.update(rakes)
    params.update(np_weights)

    hds, hdsw = extract_source_hypocentral_depths(src)
    params.update(hds)
    params.update(hdsw)

    pstrikes, pdips = extract_source_planes_strikes_dips(src)
    params.update(pstrikes)
    params.update(pdips)
    params['sourcetype'] = striptag(src.tag)
    w.record(**params)

def set_area_geometry(w, src):
    """
    Set area polygon as shapefile geometry
    """
    assert "areaSource" in src.tag
    geometry_node = src.nodes[get_taglist(src).index("areaGeometry")]
    area_attrs = parse_area_geometry(geometry_node)
    w.poly(parts=[area_attrs["polygon"].tolist()])


def set_point_geometry(w, src):
    """
    Set point location as shapefile geometry.
    """
    assert "pointSource" in src.tag
    geometry_node = src.nodes[get_taglist(src).index("pointGeometry")]
    point_attrs = parse_point_geometry(geometry_node)
    w.point(point_attrs["point"][0], point_attrs["point"][1])


def set_simple_fault_geometry(w, src):
    """
    Set simple fault trace coordinates as shapefile geometry.

    :parameter w:
        Writer
    :parameter src:
        
    """
    assert "simpleFaultSource" in src.tag
    geometry_node = src.nodes[get_taglist(src).index("simpleFaultGeometry")]
    fault_attrs = parse_simple_fault_geometry(geometry_node)
    w.line(parts=[fault_attrs["trace"].tolist()])

def build_polygon_from_fault_attrs(w, fault_attrs):
    """

    """
    fault_trace = geo.line.Line([geo.point.Point(row[0], row[1])
                                 for row in fault_attrs["trace"]])
    lon, lat = SimpleFaultSurface.get_surface_vertexes(
        fault_trace,
        fault_attrs["upperSeismoDepth"],
        fault_attrs["lowerSeismoDepth"],
        fault_attrs["dip"])
    # Reorder the vertexes
    lons = numpy.concatenate([lon[::2], numpy.flipud(lon[1::2])])
    lats = numpy.concatenate([lat[::2], numpy.flipud(lat[1::2])])
    depths = numpy.concatenate([
        numpy.ones_like(lon[::2]) * fault_attrs["upperSeismoDepth"],
        numpy.ones_like(lon[::2]) * fault_attrs["lowerSeismoDepth"]])
    # Create the 3D polygon
    w.poly(parts=[[[tlon, tlat, tdep] for tlon, tlat, tdep
                  in zip(list(lons), list(lats), list(depths))]])

def set_simple_fault_geometry_3D(w, src):
    """
    Builds a 3D polygon from a node instance
    """
    assert "simpleFaultSource" in src.tag
    geometry_node = src.nodes[get_taglist(src).index("simpleFaultGeometry")]
    fault_attrs = parse_simple_fault_geometry(geometry_node)
    build_polygon_from_fault_attrs(w, fault_attrs)

def build_lineset_from_complex_fault_attrs(w, fault_attrs):
    """

    """
    if len(fault_attrs["intermediateEdges"]):
        intermediate_edges = [edge.tolist()
                              for edge in fault_attrs["intermediateEdges"]]
    else:
        intermediate_edges = []
    edges = [fault_attrs["faultTopEdge"].tolist()] +\
             intermediate_edges +\
             [fault_attrs["faultBottomEdge"].tolist()]
    w.line(parts=edges)

def set_complex_fault_geometry(w, src):
    """

    """
    assert "complexFaultSource" in src.tag
    geometry_node = src.nodes[get_taglist(src).index("complexFaultGeometry")]
    fault_attrs = parse_complex_fault_geometry(geometry_node)
    build_lineset_from_complex_fault_attrs(w, fault_attrs)

def set_characteristic_geometry(w_simple, w_simple_3D, w_complex, w_planar,
        src):
    """

    """
    assert "characteristicFaultSource" in src.tag
    # Build from surface node
    pparts = []
    tags = get_taglist(src)
    surface = src.nodes[tags.index("surface")]
    for node in surface:
        if "complexFaultGeometry" in node.tag:
            cf_parts = []
            for subnode in node.nodes:
                crds = numpy.array(
                    subnode.nodes[0].nodes[0].text).reshape(-1, 3)
                cf_parts.append(crds.tolist())
            w_complex.line(parts=cf_parts)
        elif "simpleFaultGeometry" in node.tag:
            fault_attrs = parse_simple_fault_geometry(node)
            w_simple.line(parts=[fault_attrs["trace"].tolist()])
            build_polygon_from_fault_attrs(w_simple_3D, fault_attrs)
        elif "planarSurface" in node.tag:
            tag_list = get_taglist(node)
            #tag_list = [subnode.tag[37:] for subnode in node.nodes]
            indices = [tag_list.index("topLeft"),
                       tag_list.index("topRight"),
                       tag_list.index("bottomRight"),
                       tag_list.index("bottomLeft")]
           
            lons = [float(node.nodes[iloc].attrib["lon"]) for iloc in indices]
            lats = [float(node.nodes[iloc].attrib["lat"]) for iloc in indices]
            depths = [float(node.nodes[iloc].attrib["depth"])
                      for iloc in indices]
            pparts.append(
                [[lon, lat, depth]
                for lon, lat, depth in zip(lons, lats, depths)])
        else:
            # Nope - nothing left
            pass
    if len(pparts):
        w_planar.poly(parts=pparts)


def record_to_dict(record, fields):
    """

    """
    data = []
    for name, attr in zip(fields, record):
        value = attr.strip()
        if value:
            data.append((name, value))
    return OrderedDict(data)


def area_geometry_from_shp(shape, record):
    """
    """
    assert record["sourcetype"] == "areaSource"
    geom = []
    for row in shape.points:
        geom.extend([row[0], row[1]])
    poslist_node = LiteralNode("posList", text=geom)
    linear_ring_node = LiteralNode("LinearRing", nodes=[poslist_node])
    exterior_node = LiteralNode("exterior", nodes=[linear_ring_node])
    polygon_node = LiteralNode("Polygon", nodes=[exterior_node])

    upper_depth_node = LiteralNode("upperSeismoDepth",
                                   text=float(record["usd"]))
    lower_depth_node = LiteralNode("lowerSeismoDepth",
                                   text=float(record["lsd"]))
    return LiteralNode(
        "areaGeometry",
        nodes=[polygon_node, upper_depth_node, lower_depth_node])

def point_geometry_from_shp(shape, record):
    """
    Retrieves a point geometry
    """
    assert record["sourcetype"] == "pointSource"
    xy = shape.points[0][0], shape.points[0][1]
    pos_node = LiteralNode("pos", text=xy)
    point_node = LiteralNode("Point", nodes=[pos_node])
    upper_depth_node = LiteralNode("upperSeismoDepth",
                                   text=float(record["usd"]))
    lower_depth_node = LiteralNode("lowerSeismoDepth",
                                   text=float(record["lsd"]))
    return LiteralNode(
        "pointGeometry",
        nodes=[point_node, upper_depth_node, lower_depth_node])

def simple_fault_geometry_from_shp(shape, record):
    """

    """
    assert record["sourcetype"] == "simpleFaultSource"
    geom = []
    for row in shape.points:
        geom.extend([row[0], row[1]])
    poslist_node = LiteralNode("posList", text=geom)
    trace_node = LiteralNode("LineString", nodes=[poslist_node])
    dip_node = LiteralNode("dip", text=float(record["dip"]))
    upper_depth_node = LiteralNode("upperSeismoDepth",
                                   text=float(record["usd"]))
    lower_depth_node = LiteralNode("lowerSeismoDepth",
                                   text=float(record["lsd"]))
    return LiteralNode("simpleFaultGeometry",
                       nodes=[trace_node, dip_node, upper_depth_node,
                              lower_depth_node])

def complex_fault_geometry_from_shp(shape, record):
    """

    """
    assert record["sourcetype"] == "complexFaultSource"
    breakers = shape.parts
    breakers.append(len(shape.z))
    indices = [range(breakers[i], breakers[i + 1])
               for i in range(0, len(breakers) - 1)]
    edges = []
    for iloc, idx in enumerate(indices):
        geom = []
        for j in idx:
            geom.extend([shape.points[j][0], shape.points[j][1],
                         shape.z[j]])
        poslist_node = LiteralNode("posList", text=geom)
        linestring_node = LiteralNode("LineString", nodes=[poslist_node])
        if iloc == 0:
            # Fault top edge
            edges.append(LiteralNode("faultTopEdge", nodes=[linestring_node]))
        elif iloc == (len(indices) - 1):
            # Fault bottom edges
            edges.append(LiteralNode("faultBottomEdge",
                                     nodes=[linestring_node]))
        else:
            edges.append(LiteralNode("intermediateEdge",
                                     nodes=[linestring_node]))
    return LiteralNode("complexFaultGeometry", nodes=edges)

def build_incremental_mfd_from_shp(record):
    """

    """
    rates = []
    for i in range(1, MAX_RATES + 1):
        key = "rate{:s}".format(str(i))
        if key in record:
            rates.append(record[key])
    occur_rates = LiteralNode("occurRates", text=" ".join(rates))
    return LiteralNode("incrementalMFD",
        {"binWidth": float(record["bin_width"]),
         "minMag": float(record["min_mag"])},
        nodes=[occur_rates])

def build_trunc_gr_from_shp(record):
    """

    """
    attribs = {"aValue": float(record["a_val"]),
               "bValue": float(record["b_val"]),
               "minMag": float(record["min_mag"]),
               "maxMag": float(record["max_mag"])}
    return LiteralNode("truncGutenbergRichterMFD", attribs)

def build_mfd_from_shp(record):
    """

    """
    if ("a_val" in record) and ("b_val" in record):
        # Is truncated GR
        return build_trunc_gr_from_shp(record)
    elif ("min_mag" in record) and ("bin_width" in record):
        return build_incremental_mfd_from_shp(record)
    else:
        raise ValueError("MFD type unsupported or incomplete for source %s"
                         % record["id"])

def build_npd_from_shp(record):
    """

    """
    npds = []
    for iloc in range(1, MAX_NODAL_PLANES + 1):
        strike_key = "strike{:s}".format(str(iloc))
        dip_key = "dip{:s}".format(str(iloc))
        rake_key = "rake{:s}".format(str(iloc))
        weight_key = "npweight{:s}".format(str(iloc))
        if (strike_key in record) and (dip_key in record) and\
            (rake_key in record) and (weight_key in record):
            attribs = {"strike": float(record[strike_key]),
                       "dip": float(record[dip_key]),
                       "rake": float(record[rake_key]),
                       "probability": float(record[weight_key])}
            npds.append(LiteralNode("nodalPlane", attribs))
    return LiteralNode("nodalPlaneDist", nodes=npds)

def build_hdd_from_shp(record):
    """

    """
    hdds = []
    for iloc in range(0, MAX_HYPO_DEPTHS):
        hd_key = "hd{:s}".format(str(iloc))
        weight_key = "hdweight{:s}".format(str(iloc))
        if (hd_key in record) and (weight_key in record):
            hd_text = (float(record[weight_key]), float(record[hd_key]))
            hdds.append(LiteralNode("hypoDepth",
                                    {"depth": record[hd_key],
                                     "probability": record[weight_key]},
                                     text=hd_text))
    return LiteralNode("hypoDepthDist", nodes=hdds)

def build_point_source_from_shp(shape, record):
    """
    Builds the full point source from shapefile
    """
    attribs = {"id": record["id"], "name": record["name"],
               "tectonicRegion": record["trt"]}
    nodes = [point_geometry_from_shp(shape, record)]
    # Scaling relation
    nodes.append(LiteralNode("magScaleRel", text=record["msr"]))
    # Aspect ratio
    nodes.append(LiteralNode("ruptAspectRatio", text=float(record["rar"])))
    # MFD
    nodes.append(build_mfd_from_shp(record))
    # Nodal Plane Distribution
    nodes.append(build_npd_from_shp(record))
    # Hypocentre Depth Distribtion
    nodes.append(build_hdd_from_shp(record))
    return LiteralNode("pointSource", attribs, nodes=nodes)


def build_area_source_from_shp(shape, record):
    """
    Builds the full area source from shapefile
    """
    attribs = {"id": record["id"], "name": record["name"],
               "tectonicRegion": record["trt"]}
    nodes = [area_geometry_from_shp(shape, record)]
    # Scaling relation
    nodes.append(LiteralNode("magScaleRel", text=record["msr"]))
    # Aspect ratio
    nodes.append(LiteralNode("ruptAspectRatio", text=float(record["rar"])))
    # MFD
    nodes.append(build_mfd_from_shp(record))
    # Nodal Plane Distribution
    nodes.append(build_npd_from_shp(record))
    # Hypocentre Depth Distribtion
    nodes.append(build_hdd_from_shp(record))
    return LiteralNode("areaSource", attribs, nodes=nodes)

def build_simple_fault_source_from_shp(shape, record):
    """
    
    """
    attribs = {"id": record["id"], "name": record["name"],
               "tectonicRegion": record["trt"]}
    nodes = [simple_fault_geometry_from_shp(shape, record)]
    # Scaling relation
    nodes.append(LiteralNode("magScaleRel", text=record["msr"]))
    # Aspect ratio
    nodes.append(LiteralNode("ruptAspectRatio", text=float(record["rar"])))
    # MFD
    nodes.append(build_mfd_from_shp(record))
    # Rake
    nodes.append(LiteralNode("rake", text=float(record["rake"])))
    return LiteralNode("simpleFaultSource", attribs, nodes=nodes)


def build_complex_fault_source_from_shp(shape, record):
    """
    
    """
    attribs = {"id": record["id"], "name": record["name"],
               "tectonicRegion": record["trt"]}
    nodes = [complex_fault_geometry_from_shp(shape, record)]
    # Scaling relation
    nodes.append(LiteralNode("magScaleRel", text=record["msr"]))
    # Aspect ratio
    nodes.append(LiteralNode("ruptAspectRatio", text=float(record["rar"])))
    # MFD
    nodes.append(build_mfd_from_shp(record))
    # Rake
    nodes.append(LiteralNode("rake", text=float(record["rake"])))
    return LiteralNode("complexFaultSource", attribs, nodes=nodes)



class SourceModel(object):
    """

    """
    def __init__(self, sources, name=None):
        """

        """
        self.sources = sources
        self.name = name
        self.has_area_source = False
        self.has_point_source = False
        self.has_simple_fault_geometry = False
        self.has_complex_fault_geometry = False
        self.has_planar_geometry = False
        self.has_mfd_gr = False
        self.has_mfd_incremental = False
        self.num_r = 0
        self.num_np = 0
        self.num_hd = 0
        self.num_p = 0
        self.appraise_source_model()

    def appraise_source_model(self):
        """
        Identify parameters defined in NRML source model file, so that
        shapefile contains only source model specific fields.
        """
        for src in self.sources:
            # source params
            src_taglist = get_taglist(src)
            if "areaSource" in src.tag:
                # this is true also for area sources
                self.has_area_source = True
                npd_node = src.nodes[src_taglist.index("nodalPlaneDist")]
                npd_size = len(npd_node)
                hdd_node = src.nodes[src_taglist.index("hypoDepthDist")]
                hdd_size = len(hdd_node)
                self.num_np = npd_size if npd_size > self.num_np else\
                    self.num_np
                self.num_hd = hdd_size if hdd_size > self.num_hd else\
                    self.num_hd
            elif "pointSource" in src.tag:
                self.has_point_source = True
                npd_node = src.nodes[src_taglist.index("nodalPlaneDist")]
                npd_size = len(npd_node)
                hdd_node = src.nodes[src_taglist.index("hypoDepthDist")]
                hdd_size = len(hdd_node)
                self.num_np = npd_size if npd_size > self.num_np else\
                    self.num_np
                self.num_hd = hdd_size if hdd_size > self.num_hd else\
                    self.num_hd
            elif "simpleFaultSource" in src.tag:
                self.has_simple_fault_geometry = True
            elif "complexFaultSource" in src.tag:
                self.has_complex_fault_geometry = True
            elif "characteristicFaultSource" in src.tag:
                # Get the surface node
                surface_node = src.nodes[src_taglist.index("surface")]
                p_size = 0
                for surface in surface_node.nodes:
                    if "simpleFaultGeometry" in surface.tag:
                        self.has_simple_fault_geometry = True
                    elif "complexFaultGeometry" in surface.tag:
                        self.has_complex_fault_geometry = True
                    elif "planarSurface" in surface.tag:
                        self.has_planar_geometry = True
                        p_size += 1
                self.num_p = p_size if p_size > self.num_p else self.num_p
            else:
                pass

            # MFD params
            if "truncGutenbergRichterMFD" in src_taglist:
                self.has_mfd_gr = True

            elif "incrementalMFD" in src_taglist:
                self.has_mfd_incremental = True
                # Get rate size
                mfd_node = src.nodes[src_taglist.index("incrementalMFD")]
                r_size = len(mfd_node.nodes[0].text)
                self.num_r = r_size if r_size > self.num_r else self.num_r
                
            else:
                pass

    def __len__(self):
        """
        Return number of sources
        """
        return len(self.sources)

    def __iter__(self):
        """
        Iterate over sources
        """
        for source in self.sources:
            yield(source)


class SourceModelParser(object):
    """
    Base class executes simple export to NRML
    """
    def __init__(self):
        """
        
        """
        self.source_file = None
        self.destination = None 
    
    def read(self, nrml_file, validate=False,
            simple_fault_spacing=1.0, complex_mesh_spacing=5.0,
            mfd_spacing=0.1):
        """
        Build the source model from nrml format
        """
        self.source_file = nrml_file
        if validate:
            converter = SourceConverter(1.0, simple_fault_spacing,
                                        complex_mesh_spacing,
                                        mfd_spacing,
                                        10.0)
        src_nodes = read_nodes(nrml_file, lambda elem: "Source" in elem.tag,
                               nrml.nodefactory["sourceModel"])
        sources = []
        for no, src_node in enumerate(src_nodes, 1):
            if validate:
                print "Validating Source %s" % src_node.attrib["id"]
                _ = converter.convert_node(src_node)
            sources.append(src_node)
        return SourceModel(sources)

    def write(self, destination, source_model, name=None):
        """
        Exports to NRML
        """
        if os.path.exists(destination):
            os.remove(destination)
        self.destination = destination
        #assert isinstance(source_model, SourceModel) and len(source_model)
        if name:
            source_model.name = name
        output_source_model = LiteralNode("sourceModel", {"name": name},
                                          nodes=source_model.sources)
        print "Exporting Source Model to %s" % self.destination
        with open(self.destination, "w") as f:
            nrml.write([output_source_model], f, "%s")


class ShapefileParser(SourceModelParser):
    """

    """ 
    def filter_params(self, src_mod):
        """
        Remove params uneeded by source_model
        """
        # point and area related params
        STRIKE_PARAMS[src_mod.num_np:] = []
        DIP_PARAMS[src_mod.num_np:] = []
        RAKE_PARAMS[src_mod.num_np:] = []
        NPW_PARAMS[src_mod.num_np:] = []
        HDEPTH_PARAMS[src_mod.num_hd:] = []
        HDW_PARAMS[src_mod.num_hd:] = []
        # planar rupture related params
        PLANES_STRIKES_PARAM[src_mod.num_p:] = []
        PLANES_DIPS_PARAM[src_mod.num_p:] = []
        # rate params
        RATE_PARAMS[src_mod.num_r:] = []

        if src_mod.has_simple_fault_geometry is False:
            GEOMETRY_PARAMS.remove(('dip', 'dip', 'f'))

        
        if src_mod.has_simple_fault_geometry is False and\
            src_mod.has_complex_fault_geometry is False and\
            src_mod.has_planar_geometry is False:
            BASE_PARAMS.remove(('rake', 'rake', 'f'))

        if src_mod.has_simple_fault_geometry is False and\
            src_mod.has_complex_fault_geometry is False and\
                src_mod.has_area_source is False and\
                src_mod.has_point_source is False:
            GEOMETRY_PARAMS[:] = []

        if src_mod.has_mfd_gr is False:
            MFD_PARAMS.remove(('max_mag', 'max_mag', 'f'))
            MFD_PARAMS.remove(('a_val', 'a_val', 'f'))
            MFD_PARAMS.remove(('b_val', 'b_val', 'f'))

        if src_mod.has_mfd_incremental is False:
            MFD_PARAMS.remove(('binWidth', 'bin_width', 'f'))

    def read(self, input_shapefile, validate=False,
            simple_fault_spacing=1.0, complex_mesh_spacing=5.0,
            mfd_spacing=0.1):
        """
        Build the source model from nrml format
        """
        reader = shapefile.Reader(input_shapefile)
        fields = [field[0] for field in reader.fields[1:]]
        shapes = reader.shapes()
        records = reader.records()
        sources = []
        if validate:
            converter = SourceConverter(1.0, simple_fault_spacing,
                                        complex_mesh_spacing,
                                        mfd_spacing,
                                        10.0)

        for iloc in range(0, reader.numRecords):
            # Build record dictionary
            record = record_to_dict(records[iloc], fields)
            shape = shapes[iloc]
            if "pointSource" in record["sourcetype"]:
                src = build_point_source_from_shp(shape, record)
            elif "areaSource" in record["sourcetype"]:
                src = build_area_source_from_shp(shape, record)
            elif "simpleFaultSource" in record["sourcetype"]:
                src = build_simple_fault_source_from_shp(shape, record)
            elif "complexFaultSource" in record["sourcetype"]:
                src = build_complex_fault_source_from_shp(shape, record)
            elif "characteristicFaultSource" in record["sourcetype"]:
                print "Characteristic Fault Source Not Yet Supported - Sorry!"
            else:
                pass
            if validate:
                print "Validating Source %s" % src.attrib["id"]
                _ = converter.convert_node(src)
            sources.append(src)
        return SourceModel(sources)

    def write(self, destination, source_model, name=None):
        """
        Save sources - to multiple
        shapefiles corresponding to different source typolgies/geometries
        ('_point', '_area', '_simple', '_complex', '_planar')
        """
        if os.path.exists(destination + ".shp"):
            os.system("rm %s.*" % destination)
        self.destination = destination
        #field_flags = appraise_nrml_source_model(source_model)
        #filter_params(*field_flags)
        self.filter_params(source_model)

        w_area = shapefile.Writer(shapefile.POLYGON)
        w_point = shapefile.Writer(shapefile.POINT)
        w_simple = shapefile.Writer(shapefile.POLYLINE)
        w_simple3d = shapefile.Writer(shapefile.POLYGONZ)
        w_complex = shapefile.Writer(shapefile.POLYLINEZ)
        w_planar = shapefile.Writer(shapefile.POLYGONZ)

        register_fields(w_area)
        register_fields(w_point)
        register_fields(w_simple)
        register_fields(w_simple3d)
        register_fields(w_complex)
        register_fields(w_planar)

        #srcm = SourceModelParser(source_model).parse()
        for src in source_model.sources:

            # Order is important here
            if "areaSource" in src.tag:
                set_params(w_area, src)
                set_area_geometry(w_area, src)
            elif "pointSource" in src.tag:
                set_params(w_point, src)
                set_point_geometry(w_point, src)
            elif "complexFaultSource" in src.tag:
                set_params(w_complex, src)
                set_complex_fault_geometry(w_complex, src)
            elif "simpleFaultSource" in src.tag:
                set_params(w_simple, src)
                set_simple_fault_geometry(w_simple, src)
                # Create the 3D polygon
                set_params(w_simple3d, src)
                set_simple_fault_geometry_3D(w_simple3d, src)
            elif "characteristicFaultSource" in src.tag:
                src_taglist = get_taglist(src)
                surface_node = src.nodes[src_taglist.index("surface")]
                for subnode in surface_node:
                    if "simpleFaultGeometry" in subnode.tag:
                        set_params(w_simple, src)
                        set_params(w_simple3d, src)
                    elif "complexFaultGeometry" in subnode.tag:
                        set_params(w_complex, src)
                    elif "planarSurface" in subnode.tag:
                        set_params(w_planar, src)
                    else:
                        raise ValueError(
                            'Geometry class %s not recognized'
                            % subnode.tag)
                    set_characteristic_geometry(w_simple, w_simple3d,
                                                w_complex, w_planar, src)

            else:
                raise ValueError('Source type %s not recognized'
                                 % src.tag)

        root = self.destination

        if len(w_area.shapes()) > 0:
            w_area.save('%s_area' % root)
        if len(w_point.shapes()) > 0:
            w_point.save('%s_point' % root)
        if len(w_complex.shapes()) > 0:
            w_complex.save('%s_complex' % root)
        if len(w_simple.shapes()) > 0:
            w_simple.save('%s_simple' % root)
            w_simple3d.save('%s_simple3d' % root)
        if len(w_planar.shapes()) > 0:
            w_planar.save('%s_planar' % root)

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """

    description = ('Convert NRML source model file to ESRI Shapefile(s) and '
                   'vice versa.\n\nTo convert from NRML to shapefile type: '
                   '\npython source_model_converter.py '
                   '--input-nrml-file PATH_TO_SOURCE_MODEL_NRML_FILE. '
                   '--output-file PATH_TO_OUTPUT_FILE'
                   '\n\nFor each type of source geometry defined in the NRML '
                   'file (point, area, simple fault, complex fault, planar) '
                   'a separate shapefile is created. Each shapefile is '
                   'differentiated by a specific ending'
                   '(\'_point\', \'_area\', \'_simple\', \'_complex\', '
                   '\'_planar\')'
                   '\n\nTo convert from shapefile(s) to NRML type: '
                   '\npython source_model_converter.py '
                   '--input-shp-files PATH_TO_SOURCE_MODEL_SHP_FILE1 '
                   'PATH_TO_SOURCE_MODEL_SHP_FILE2 ...'
                   '--output-file PATH_TO_OUTPUT_FILE'
                   '\n\nSources defined in different shapefile are saved'
                   ' into a single NRML file.')
    parser = argparse.ArgumentParser(description=description,
                                     add_help=False,
                                     formatter_class=RawTextHelpFormatter)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--output-file', help='path to output file (root' +
                       ' name only)',
                       default=None,
                       required=True)
    group = flags.add_mutually_exclusive_group()
    group.add_argument('--input-nrml-file',
                       help='path to source model NRML file',
                       default=None)
    parser.add_argument('--validate',
                        dest='validate',
                        action="store_true",
                        help='Apply validation to input model (can be slow)',
                        default=False,
                        required=False)

    group.add_argument('--input-shp-files',
                       help='path(s) to source model ESRI shapefile(s)' +
                            '(file root only - no extension)',
                       nargs='+',
                       default=None)

    return parser

if __name__ == "__main__":
    parser = set_up_arg_parser()
    args = parser.parse_args()
    if args.input_nrml_file:
        input_parser = SourceModelParser()
        source_model = input_parser.read(args.input_nrml_file,
                                         args.validate)
        # Shapefile parser
        shapefile_parser = ShapefileParser()
        shapefile_parser.write(args.output_file, source_model)
    elif args.input_shp_files:
        input_parser = ShapefileParser()
        for iloc, filename in enumerate(args.input_shp_files):
            if iloc == 0:
                source_model = input_parser.read(filename,
                                                 args.validate)
            else:
                tmp = input_parser.read(filename, args.validate)
                source_model.sources.extend(tmp.sources)
        nrml_writer = SourceModelParser()
        nrml_writer.write(args.output_file, source_model)
    else:
        parser.print_usage()
