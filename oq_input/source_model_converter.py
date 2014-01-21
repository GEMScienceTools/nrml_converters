"""
Convert NRML source model file to ESRI shapefile.
"""
import os
import argparse
import shapefile
from shapely import wkt

from openquake.nrmllib.hazard.parsers import SourceModelParser
from openquake.nrmllib.models import (PointSource, PointGeometry, AreaSource,
    AreaGeometry, SimpleFaultSource, SimpleFaultGeometry, ComplexFaultSource,
    ComplexFaultGeometry, IncrementalMFD, TGRMFD, NodalPlane, HypocentralDepth,
    CharacteristicSource, PlanarSurface, Point)

# maximum field size allowed by shapefile
FIELD_SIZE = 255
# maximum number of occurrence rates that can be stored for incremental MFD
MAX_RATES = 50
# maximum number of nodal planes
MAX_NODAL_PLANES = 20
# maximum number of hypocentral depths
MAX_HYPO_DEPTHS = 20

# each triplet contains hazardlib parameter name, shapefile field name and
# data type
BASE_PARAMS = [
    ('id', 'id', 'c'), ('name', 'name', 'c'), ('trt', 'trt', 'c'),
    ('mag_scale_rel', 'msr', 'c'), ('rupt_aspect_ratio', 'rar', 'f')
]
GEOMETRY_PARAMS = [
    ('upper_seismo_depth', 'usd', 'f'), ('lower_seismo_depth', 'lsd', 'f'),
    ('dip', 'dip', 'f')
]
MFD_PARAMS = [
    ('min_mag', 'min_mag', 'f'), ('max_mag', 'max_mag', 'f'),
    ('a_val', 'a_val', 'f'), ('b_val', 'b_val', 'f'),
    ('bin_width', 'bin_width', 'f')
]

# shapefile specific fields
RATE_PARAMS = [('rate%s' % (i+1), 'f') for i in range(MAX_RATES)]
STRIKE_PARAMS = [('strike%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
DIP_PARAMS = [('dip%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
RAKE_PARAMS = [('rake%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
NPW_PARAMS = [('np_weight%s' % (i+1), 'f') for i in range(MAX_NODAL_PLANES)]
HDEPTH_PARAMS = [('hd%s' % (i+1), 'f') for i in range(MAX_HYPO_DEPTHS)]
HDW_PARAMS = [('hd_weight%s' % (i+1), 'f') for i in range(MAX_HYPO_DEPTHS)]

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
        HDEPTH_PARAMS, HDW_PARAMS
    ]
    for PARAMS in PARAMS_LIST:
        for param, dtype in PARAMS:
            w.field(param, fieldType=dtype, size=FIELD_SIZE)

def expand_src_param(values, shp_params):
    """
    Expand hazardlib source attribute (defined through list of values)
    into dictionary of shapefile parameters.
    """
    if values is None:
        return dict([(key, None) for key, _ in shp_params])
    else:
        num_values = len(values)
        return dict(
            [(key, float(values[i]) if i < num_values else None)
            for i, (key, _) in enumerate(shp_params)]
        )

def check_size(values, name, MAX):
    """
    Raise error if size of a list is larger than allowed.
    """
    num_values = len(values)
    if values is not None and num_values > MAX:
        raise ValueError('Number of values in NRML file for %s'
            'is too large for being saved in shapefile.' % name)

def extract_source_params(obj, PARAMS):
    """
    Extract params from source object.
    """
    return dict(
        [(param, getattr(obj, key, None)) for key, param, _ in PARAMS]
    )

def extract_source_rates(src):
    """
    Extract source occurrence rates.
    """
    rates = getattr(src.mfd, 'occur_rates', None)
    check_size(rates, 'occur_rates', MAX_RATES)

    return expand_src_param(rates, RATE_PARAMS)

def extract_source_nodal_planes(src):
    """
    Extract source nodal planes.
    """
    nodal_planes = getattr(src, 'nodal_plane_dist', None)
    check_size(nodal_planes, 'nodal_plane_dist', MAX_NODAL_PLANES)

    strikes = [np.strike for np in nodal_planes]
    dips = [np.dip for np in nodal_planes]
    rakes = [np.rake for np in nodal_planes]
    np_weights = [np.probability for np in nodal_planes]

    strikes = expand_src_param(strikes, STRIKE_PARAMS)
    dips = expand_src_param(dips, DIP_PARAMS)
    rakes = expand_src_param(rakes, RAKE_PARAMS)
    np_weights = expand_src_param(np_weights, NPW_PARAMS)

    return strikes, dips, rakes, np_weights

def extract_hypocentral_depths(src):
    """
    Extrac source hypocentral depths.
    """
    hypo_depths = getattr(src, 'hypo_depth_dist', None)
    check_size(hypo_depths, 'hypo_depths', MAX_HYPO_DEPTHS)

    hds = [hd.depth for hd in hypo_depths]
    hdws = [hd.probability for hd in hypo_depths]

    hds = expand_src_param(hds, HDEPTH_PARAMS)
    hdsw = expand_src_param(hdws, HDW_PARAMS)

    return hds, hdsw

def set_params(w, src):
    """
    Set source parameters.
    """
    params = extract_source_params(src, BASE_PARAMS)

    params.update(extract_source_params(src.geometry, GEOMETRY_PARAMS))

    params.update(extract_source_params(src.mfd, MFD_PARAMS))

    params.update(extract_source_rates(src))

    strikes, dips, rakes, np_weights = extract_source_nodal_planes(src)
    params.update(strikes)
    params.update(dips)
    params.update(rakes)
    params.update(np_weights)

    hds, hdsw = extract_hypocentral_depths(src)
    params.update(hds)
    params.update(hdsw)

    w.record(**params)

def set_area_geometry(w, src):
    """
    Set area polygon as shapefile geometry
    """
    coords = wkt.loads(src.geometry.wkt)
    lons, lats = coords.exterior.xy

    w.poly(parts=[[[lon, lat] for lon, lat in zip(lons, lats)]])

def save_area_srcs_to_shp(source_model):
    """
    Save area sources to ESRI shapefile.
    """
    w = shapefile.Writer(shapefile.POLYGON)
    register_fields(w)

    srcm = SourceModelParser(source_model).parse()
    for src in srcm:
        if isinstance(src, AreaSource):
            set_params(w, src)
            set_area_geometry(w, src)

    if w.shapes > 0:
        root = os.path.splitext(source_model)[0]
        w.save('%s_area' % root)

def nrml2shp(source_model):
    """
    Convert NRML source model file to ESRI shapefile
    """
    save_area_srcs_to_shp(source_model)

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML source model file to ESRI Shapefile. '
            'To run just type: python source_model_converter.py '
            '--input-nrml-file=PATH_TO_SOURCE_MODEL_NRML_FILE ',
            add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-nrml-file',
        help='path to source model NRML file (Required)',
        default=None,
        required=True)

    return parser

if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_nrml_file:
        nrml2shp(args.input_nrml_file)
    else:
        parser.print_usage()
