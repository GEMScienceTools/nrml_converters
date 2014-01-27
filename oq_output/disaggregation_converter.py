#!/usr/bin/env/python
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
'''
Convert NRML disaggregation file to .csv file and plot disaggregation matrices
using GMT
'''

import os
import argparse
import numpy
from lxml import etree
from collections import OrderedDict
from subprocess import call

NRML='{http://openquake.org/xmlns/nrml/0.4}'

def parse_nrml_disaggregation_file(nrml_disaggregation):
    """
    Parse NRML disaggregation file.
    """
    metadata = OrderedDict()
    matrices = {}

    parse_args = dict(source=nrml_disaggregation)
    for _, element in etree.iterparse(**parse_args):
        if element.tag == '%sdisaggMatrices' % NRML:
            a = element.attrib
            metadata['smlt_path'] = a.get('sourceModelTreePath')
            metadata['gsimlt_path'] = a.get('gsimTreePath')
            metadata['imt'] = a['IMT']
            metadata['investigation_time'] = a['investigationTime']
            metadata['sa_period'] = a.get('saPeriod')
            metadata['sa_damping'] = a.get('saDamping')
            metadata['lon'] = a.get('lon')
            metadata['lat'] = a.get('lat')
            metadata['Mag'] = \
                numpy.array(a.get('magBinEdges').split(','), dtype=float)
            metadata['Dist'] = \
                numpy.array(a.get('distBinEdges').split(','), dtype=float)
            metadata['Lon'] = \
                numpy.array(a.get('lonBinEdges').split(','), dtype=float)
            metadata['Lat'] = \
                numpy.array(a.get('latBinEdges').split(','), dtype=float)
            metadata['Eps'] = \
                numpy.array(a.get('epsBinEdges').split(','), dtype=float)
            metadata['TRT'] = \
                numpy.array(
                    a.get('tectonicRegionTypes').split(','), dtype=object
                )
        elif element.tag == '%sdisaggMatrix' % NRML:
            a = element.attrib
            disag_type = a.get('type')
            dims = tuple(map(int, a.get('dims').split(',')))
            poe = float(a.get('poE'))
            iml = float(a.get('iml'))

            matrix = numpy.zeros(dims)
            for e in element:
                a = e.attrib
                idx = tuple(map(int, a.get('index').split(',')))
                value = float(a.get('value'))
                matrix[idx] = value

            matrices[disag_type] = (poe, iml, matrix)

    return metadata, matrices

def save_disagg_to_csv(nrml_disaggregation, output_dir, plot):
    """
    Save disaggregation matrices to multiple .csv files.
    """
    metadata, matrices = parse_nrml_disaggregation_file(nrml_disaggregation)

    skip_keys = ('Mag', 'Dist', 'Lon', 'Lat', 'Eps', 'TRT')

    base_header = ','.join(
        '%s=%s' % (key, value) for key, value in metadata.items()
        if value is not None and key not in skip_keys
    )

    for disag_type, (poe, iml, matrix) in matrices.items():
        header = '# %s,poe=%s,iml=%s\n' % (base_header, poe, iml)

        if disag_type == 'Mag,Lon,Lat':
            matrix = numpy.swapaxes(matrix, 0, 1)
            matrix = numpy.swapaxes(matrix, 1, 2)
            disag_type = 'Lon,Lat,Mag'

        variables = tuple(disag_type.split(','))

        axis = [metadata[v] for v in variables]

        header += ','.join(v for v in variables)
        header += ',poe'

        # compute axis mid points
        axis = [(ax[: -1] + ax[1:]) / 2.
                if ax.dtype==float else ax for ax in axis]

        values = None
        if len(axis) == 1:
            values = numpy.array([axis[0], matrix.flatten()]).T
        else:
            grids = numpy.meshgrid(*axis, indexing='ij')
            values = [g.flatten() for g in grids]
            values.append(matrix.flatten())
            values = numpy.array(values).T

        output_file = '%s/%s.csv' % (output_dir, disag_type.replace(',', '_'))
        numpy.savetxt(output_file, values, fmt='%s', delimiter=',',
                      header=header, comments='')

        if plot:
            if disag_type == 'Mag':
                plot_1d_hist(output_file, 'Magnitude', '')
            elif disag_type == 'Dist':
                plot_1d_hist(output_file, 'JB distance', '')
            elif disag_type == 'TRT':
                ntrt = metadata['TRT'].size
                bin_edges = numpy.linspace(0, ntrt, ntrt)
                annotation_file = open("annotation.dat",'w')
                for i in range(ntrt):
                    annotation_file.write("%s %s %s %s %s %s %s\n" % 
                        (bin_edges[i],
                        numpy.max(matrix) + 0.05 * numpy.max(matrix),
                        12, 0.0, 0, 'MC', metadata['TRT'][i]))
                annotation_file.close()
                plot_1d_hist(output_file, 'Tectonic Region',
                             '', annotation_file.name)
            elif disag_type == 'Mag,Dist':
                plot_2d_hist(output_file, 'Magnitude', 'JB distance', '')
            elif disag_type == 'Lon,Lat':
                plot_2d_hist(output_file, 'Longitude', 'Latitude', '')
            elif disag_type == 'Mag,Dist,Eps':
                plot_3d_hist(output_file, 'Magnitude', 'Distance', 'Epsilon', '')
            elif disag_type == 'Lon,Lat,Eps':
                plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Epsilon', '')
            elif disag_type == 'Lon,Lat,Mag':
                plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Magnitude', '')
            elif disag_type == 'Lon,Lat,TRT':
                plot_3d_hist(output_file, 'Longitude', 'Latitude', '', '')

def plot_1d_hist(hist_file, xlabel, title, annotation_file=None):
    """
    Plot 1D histogram
    """
    _, tail = os.path.split(hist_file)

    name = os.path.splitext(hist_file)[0]
    plot_file = open('%s.ps' % name,'w')

    if tail == 'TRT.csv':
        data = numpy.loadtxt(hist_file, delimiter=',', skiprows=2, dtype=object)
        data = data.reshape(-1, 2)
        x = numpy.linspace(0, 1, data.shape[0])
        y = data[:, 1].astype(numpy.float)
        bin_width = 1
    else:
        data = numpy.loadtxt(hist_file, delimiter=',', skiprows=2)
        data = data.reshape(-1, 2)
        x = data[:, 0]
        y = data[:, 1]
        bin_width = numpy.diff(x)[0]

    region = '-R%s/%s/0/%s' % \
        (numpy.min(x) - bin_width, numpy.max(x) + bin_width, numpy.max(y))
    projection = '-JX15/15'
    annotation = '-B:%s:%s/:Probability:%s:.%s:WS' % \
        (xlabel, 2 * bin_width, numpy.round(numpy.max(y)/4, 4), title)

    if tail == 'TRT.csv':
        numpy.savetxt('trt_hist.dat', numpy.array([x, y], dtype=float).T)
        annotation = '-B:%s:/:Probability:%s:.%s:WS' % \
            (xlabel, numpy.round(numpy.max(y)/4, 4), title)
        call(['psxy', 'trt_hist.dat', region, projection, annotation, '-Xc',
            '-Yc', '-Sb%su' % bin_width, '-Ggray', '-W1p', '-K'],
            stdout=plot_file)
        call(['rm', 'trt_hist.dat'])
    else:
        call(['psxy', hist_file, region, projection, annotation,'-H2', '-Xc',
            '-Yc', '-Sb%su' % bin_width, '-Ggray', '-W1p', '-K'],
            stdout=plot_file)

    if annotation_file:
        call(['pstext', annotation_file, region, projection, '-N', '-O', '-K'],
              stdout=plot_file)
        call(['rm', annotation_file])

def plot_2d_hist(hist_file, xlabel, ylabel, title):
    """
    Plot 2D histogram
    """
    name = os.path.splitext(hist_file)[0]
    plot_file = open('%s.ps' % name,'w')

    x, y, z = numpy.loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True)
    bin_width1 = numpy.diff(numpy.unique(x))[0]
    bin_width2 = numpy.diff(numpy.unique(y))[0]

    region = '-R%s/%s/%s/%s/%s/%s' % \
        (x[0] - bin_width1, x[-1] + bin_width1,
         y[0] - bin_width2, y[-1] + bin_width2,
         0.0, numpy.max(z))
    projection = '-JX18/15'
    annotation = '-B:%s:%s/:%s:%s/%s:Probability::.%s:wEsNZ' % \
        (xlabel, 2 * bin_width1, ylabel, 2 * bin_width2,
        numpy.round(numpy.max(y) / 4, 4), title)

    call(['psxyz',hist_file,region,projection,
          '-JZ8c', annotation, '-E45/35', '-m', '-So0.5', '-Wthinnest',
          '-Ggray'], stdout=plot_file)

def plot_3d_hist(hist_file, xlabel, ylabel, zlabel, title):
    """
    Plot 3d histogram
    """
    _, tail = os.path.split(hist_file)

    name = os.path.splitext(hist_file)[0]
    plot_file = open('%s.ps' % name,'w')

    if tail == 'Lon_Lat_TRT.csv':
        x, y, p = numpy.loadtxt(
            hist_file, delimiter=',', skiprows=2, unpack=True, usecols=(0, 1, 3)
        )
        z = numpy.loadtxt(
            hist_file, delimiter=',', skiprows=2, unpack=True, usecols=(2,),
            dtype=str
        )
        x_axis = numpy.unique(x)
        y_axis = numpy.unique(y)
        z_axis = numpy.arange(len(numpy.unique(z)))
    else:
        x, y, z, p = numpy.loadtxt(
            hist_file, delimiter=',', skiprows=2, unpack=True
        )
        x_axis = numpy.unique(x)
        y_axis = numpy.unique(y)
        z_axis = numpy.unique(z)

    p = p.reshape((len(x_axis), len(y_axis), len(z_axis)))

    bin_width1 = numpy.diff(x_axis)[0] if len(x_axis) > 1 else 0
    bin_width2 = numpy.diff(y_axis)[0] if len(y_axis) > 1 else 0
    bin_width3 = numpy.diff(z_axis)[0] if len(z_axis) > 1 else 1

    z_bin_edges = z_axis - bin_width3 / 2.
    z_bin_edges = numpy.append(z_bin_edges, [z_axis[-1] + bin_width3 / 2.])

    max_high = numpy.max(numpy.sum(p, axis=2))

    region = '-R%s/%s/%s/%s/%s/%s' % \
        (x_axis[0] - bin_width1, x_axis[-1] + bin_width1,
         y_axis[0] - bin_width2, y_axis[-1] + bin_width2,
         0.0, max_high)
    projection = '-JX12/12'
    annotation = '-B:%s:%s/:%s:%s/%s:Probability::.%s:wEsNZ' % \
        (xlabel, 2 * bin_width1, ylabel, 2 * bin_width2,
        numpy.round(max_high / 4, 4), title)

    cpt = open('colors.cpt', 'w')
    call(['makecpt', '-Cjet',
          '-T%s/%s/%s' % (z_bin_edges[0], z_bin_edges[-1], bin_width3), '-N'],
          stdout=cpt)
    cpt.close()

    # modify cpt file to create custom annotation
    if tail == 'Lon_Lat_TRT.csv':
        cpt_table = numpy.loadtxt('colors.cpt')
        if len(cpt_table.shape) == 1:
            cpt_table = cpt_table.reshape(-1, cpt_table.size)
        annotations = numpy.array([';%s' % trt for trt in numpy.unique(z)])
        annotations = annotations.reshape(-1, 1)
        annotations = annotations.astype(object)
        cpt_table = numpy.concatenate((cpt_table, annotations), axis=1)
        cpt = open('colors.cpt', 'w')
        for (v1, r1, g1, b1, v2, r2, g2, b2, label) in cpt_table:
            cpt.write('%f %i %i %i %f %i %i %i %s' %
                      (v1, r1, g1, b1, v2, r2, g2, b2, label))
        cpt.close()

    first = False
    for i, v1 in enumerate(x_axis):
        for j, v2 in enumerate(y_axis):
            for k, v3 in enumerate(z_axis):

                if k == 0:
                    base_hight = 0.
                else:
                    base_hight = numpy.sum(p[i, j, :k])

                f = open('values.dat', 'w')
                f.write('%s %s %s %s' % (v1, v2, base_hight + p[i, j, k], v3))
                f.close()

                if p[i, j, k] != 0 and first is False:
                    first = True
                    call(['psxyz', 'values.dat', region, projection,
                          '-JZ8c', annotation,'-E45/35','-Wthinnest',
                          '-K', '-SO0.5b%s' % base_hight, '-C%s' % cpt.name],
                          stdout=plot_file)
                elif p[i, j, k] != 0:
                    call(['psxyz', 'values.dat', region, projection,
                          '-JZ8c','-E45/35','-Wthinnest',
                          '-O', '-K', '-SO0.5b%s' % base_hight,
                          '-C%s' % cpt.name], stdout=plot_file)
                else:
                    continue

    call(['psscale', '-B:%s:' % zlabel, '-D18/7.5/10/0.5', '-Li0.3',
              '-C%s' % cpt.name, '-O'], stdout=plot_file)

    call(['rm', 'values.dat', 'colors.cpt'])

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert Disaggregation file from Nrml to .csv file. '
            'Optionally plot results using GMT.'
            'To run just type: python disaggregation_converter.py ' 
            '--input-file=/PATH/TO/INPUT_FILE '
            '--output-file=/PATH/TO/OUTPUT_FILE', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--plot', help='plot disaggregation matrices using GMT',
                        action='store_true')
    flags.add_argument('--input-file',
                        help='path to NRML disaggregation file (Required)',
                        default=None,
                        required=True)
    flags.add_argument('--output-dir',
                        help='path to output directory (Required, raise an '
                             'error if it already exists)',
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

        save_disagg_to_csv(args.input_file, args.output_dir, args.plot)
    else:
        parser.print_usage()
