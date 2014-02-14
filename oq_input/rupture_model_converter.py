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
Convert rupture model from NRML to ESRI shapefile (and viceversa)
"""
import argparse
import shapefile
from collections import OrderedDict
from argparse import RawTextHelpFormatter
from lxml import etree

from openquake.nrmllib.hazard.parsers import RuptureModelParser
from openquake.nrmllib.models import SimpleFaultRuptureModel, \
    ComplexFaultRuptureModel
from openquake.nrmllib.hazard.writers import SourceModelXMLWriter
from openquake.nrmllib import NRMLFile, SERIALIZE_NS_MAP

from source_model_converter import set_simple_fault_geometry, \
    set_complex_fault_geometry, create_simple_fault_geometry, \
    create_complex_fault_geometry


# maximum field size allowed by shapefile
FIELD_SIZE = 255

# each triplet contains nrmllib parameter name, shapefile field name and
# data type
BASE_PARAMS = [
    ('magnitude', 'mag', 'f'), ('rake', 'rake', 'f')
]
HYPO_PARAMS = [
    ('lon', 'lon', 'f'), ('lat', 'lat', 'f'), ('depth', 'depth', 'f')
]
GEOMETRY_PARAMS = [
    ('upper_seismo_depth', 'usd', 'f'), ('lower_seismo_depth', 'lsd', 'f'),
    ('dip', 'dip', 'f')
]


class RuptureModelXMLWriter(object):
    def __init__(self, dest):
        self.dest = dest

    def serialize(self, rm):
        """
        Serialize rupture model to NRML file
        """
        srcm = SourceModelXMLWriter(None)
        with NRMLFile(self.dest, 'w') as fh:
            root = etree.Element(
                'nrml', nsmap=SERIALIZE_NS_MAP
            )

            if isinstance(rm, ComplexFaultRuptureModel):
                rm_elem = etree.SubElement(root, 'complexFaultRupture')
                mag = etree.SubElement(rm_elem, 'magnitude')
                mag.text = str(rm.magnitude)
                rake = etree.SubElement(rm_elem, 'rake')
                rake.text = str(rm.rake)
                hypocenter = etree.SubElement(rm_elem, 'hypocenter',
                    attrib={'lon': rm.hypocenter[0], 'lat': rm.hypocenter[1],
                    'depth': rm.hypocenter[2]})
                srcm._append_complex_fault_geom(rm_elem, rm.geometry)
            elif isinstance(rm, SimpleFaultRuptureModel):
                rm_elem = etree.SubElement(root, 'simpleFaultRupture')
                mag = etree.SubElement(rm_elem, 'magnitude')
                mag.text = str(rm.magnitude)
                rake = etree.SubElement(rm_elem, 'rake')
                rake.text = str(rm.rake)
                hypocenter = etree.SubElement(rm_elem, 'hypocenter',
                    attrib={'lon': rm.hypocenter[0], 'lat': rm.hypocenter[1],
                    'depth': rm.hypocenter[2]})
                srcm._append_simple_fault_geom(rm_elem, rm.geometry)
            else:
                raise ValueError('Rupture model %s not recognized' % rm.__class__)

            fh.write(etree.tostring(root, pretty_print=True,
                xml_declaration=True, encoding='UTF-8'))


def register_fields(w):
    """
    Register shapefile fields.
    """
    PARAMS_LIST = [BASE_PARAMS, HYPO_PARAMS, GEOMETRY_PARAMS]
    for PARAMS in PARAMS_LIST:
        for _, param, dtype in PARAMS:
            w.field(param, fieldType=dtype, size=FIELD_SIZE)

    # source typology
    w.field('rup_type', 'C')

def nrml2shp(rupture_model, output_file):
    """
    Save NRML rupture model to ESRI shapefile.
    """
    w_simple = shapefile.Writer(shapefile.POLYLINE)
    w_complex = shapefile.Writer(shapefile.POLYLINEZ)

    register_fields(w_simple)
    register_fields(w_complex)

    rm = RuptureModelParser(rupture_model).parse()

    params = OrderedDict(
        [(p_shp, getattr(rm, p_nrml)) for p_nrml, p_shp, _ in BASE_PARAMS]
    )

    lon, lat, depth = rm.hypocenter
    params.update(
        OrderedDict([('lon', lon), ('lat', lat), ('depth', depth)])
    )

    params.update(
        [(p_shp, getattr(rm.geometry, p_nrml, None))
         for p_nrml, p_shp, _ in GEOMETRY_PARAMS]
    )

    params['rup_type'] = rm.__class__.__name__

    w_simple.record(**params)
    w_complex.record(**params)

    # order is important here
    if isinstance(rm, ComplexFaultRuptureModel):
        set_complex_fault_geometry(w_complex, rm.geometry)
    elif isinstance(rm, SimpleFaultRuptureModel):
        set_simple_fault_geometry(w_simple, rm.geometry)
    else:
        raise ValueError('Rupture model %s not recognized' % rm.__class__)

    if len(w_simple.shapes()) > 0:
        w_simple.save(output_file)
    if len(w_complex.shapes()) > 0:
        w_complex.save(output_file)

def shp2nrml(rupture_model, output_file):
    """
    Convert rupture model ESRI shapefiles to NRML.
    """
    sf = shapefile.Reader(rupture_model)

    assert len(sf.shapes()) == 1

    shape = sf.shapes()[0]
    record = sf.records()[0]

    params = OrderedDict()
    PARAMS_LIST = [BASE_PARAMS, HYPO_PARAMS, GEOMETRY_PARAMS]
    idx_p = 0
    for PARAMS in PARAMS_LIST:
        idx_v = 0
        for nrmllib_p, _, _ in PARAMS:
            params[nrmllib_p] = record[idx_p + idx_v]
            idx_v += 1
        idx_p += len(PARAMS)

    magnitude = params['magnitude']
    rake = params['rake']
    hypocenter = [params['lon'], params['lat'], params['depth']]

    if record[-1] == 'SimpleFaultRuptureModel':
        geometry = create_simple_fault_geometry(shape, params)
        rup = SimpleFaultRuptureModel(
            None, magnitude, rake, hypocenter, geometry
        )
    elif record[-1] == 'ComplexFaultRuptureModel':
        geometry = create_complex_fault_geometry(shape, params)
        rup = ComplexFaultRuptureModel(
            None, magnitude, rake, hypocenter, geometry
        )
    else:
        raise ValueError('Rupture model %s not recognized' % record[-1])

    w = RuptureModelXMLWriter('%s.xml' % output_file)
    w.serialize(rup)

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML rupture model file to ESRI Shapefile and '
            'vice versa.\n\nTo convert from NRML to shapefile type: '
            '\npython rupture_model_converter.py '
            '--input-nrml-file PATH_TO_RUPTURE_MODEL_NRML_FILE. '
            '--output-file PATH_TO_OUTPUT_FILE'
            '\n\nTo convert from shapefile to NRML type: '
            '\npython rupture_model_converter.py '
            '--input-shp-files PATH_TO_RUPTUE_MODEL_SHP_FILE '
            '--output-file PATH_TO_OUTPUT_FILE',
            add_help=False, formatter_class=RawTextHelpFormatter)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--output-file', help='path to output file (root name only)',
    default=None,
    required=True)
    group = flags.add_mutually_exclusive_group()
    group.add_argument('--input-nrml-file',
        help='path to rupture model NRML file',
        default=None)
    group.add_argument('--input-shp-file',
        help='path to rupture model ESRI shapefile (file root only - no extension)',
        default=None)
    

    return parser

if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_nrml_file:
        nrml2shp(args.input_nrml_file, args.output_file)
    elif args.input_shp_file:
        shp2nrml(args.input_shp_file, args.output_file)
    else:
        parser.print_usage()