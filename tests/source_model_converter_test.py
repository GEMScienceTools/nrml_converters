import os
import unittest
import numpy
from lxml import etree
from shapely import wkt
from shapely.geometry.polygon import Polygon
from glob import glob
from subprocess import call

from oq_input.source_model_converter import nrml2shp, shp2nrml

from openquake.nrmllib.hazard.parsers import SourceModelParser
from openquake.nrmllib.models import PointSource, AreaSource, SimpleFaultSource, \
    ComplexFaultSource, CharacteristicSource, IncrementalMFD, TGRMFD, \
    ComplexFaultGeometry, SimpleFaultGeometry, PlanarSurface

DATA_PATH = '%s/data/' % os.path.dirname(__file__)

class TestSourceModelConverter(unittest.TestCase):

    def tearDown(self):
        files = glob('%s*.shp' % DATA_PATH)
        files.extend(glob('%s*.shx' % DATA_PATH))
        files.extend(glob('%s*.dbf' % DATA_PATH))
        files.extend(glob('%s*_complex.xml' % DATA_PATH))
        files.extend(glob('%s*_simple.xml' % DATA_PATH))
        files.extend(glob('%s*_area.xml' % DATA_PATH))
        files.extend(glob('%s*_point.xml' % DATA_PATH))
        files.extend(glob('%s*_planar.xml' % DATA_PATH))

        for f in files:
            call(['rm', f])

    def assert_coordinates_equal(self, wkt1, wkt2):
        """
        Check that geometries (in WKT) are equal
        """
        geo1 = wkt.loads(wkt1)
        geo2 = wkt.loads(wkt2)

        if isinstance(geo1, Polygon):
            geo1 = [p for p in geo1.exterior.coords]
            geo2 = [p for p in geo2.exterior.coords]
        else:
            geo1 = [p for p in geo1.coords]
            geo2 = [p for p in geo2.coords]

        geo1 = numpy.array(geo1)
        geo2 = numpy.array(geo2)

        numpy.testing.assert_allclose(geo1, geo2)

    def assert_source_model_equals(self, source_model1, source_model2):
        """
        Check that source models saved in two NRML source models files
        are equal.
        """
        srcm1 = SourceModelParser(source_model1).parse()
        srcm2 = SourceModelParser(source_model2).parse()

        srcs1 = [src for src in srcm1]
        srcs2 = [src for src in srcm2]

        self.assertEqual(len(srcs1), len(srcs2))

        for src1, src2 in zip(srcs1, srcs2):
            assert isinstance(
                src1,
                (PointSource, AreaSource, SimpleFaultSource, ComplexFaultSource,
                 CharacteristicSource)
            )
            assert isinstance(
                src2,
                (PointSource, AreaSource, SimpleFaultSource, ComplexFaultSource,
                 CharacteristicSource)
            )

            self.assertEqual(src1.__class__, src2.__class__)
            self.assertEqual(src1.id, src2.id)
            self.assertEqual(src1.name, src2.name)
            self.assertEqual(src1.trt, src2.trt)

            if not isinstance(src1, CharacteristicSource):
                self.assertEqual(src1.mag_scale_rel, src2.mag_scale_rel)
                self.assertEqual(src1.rupt_aspect_ratio, src2.rupt_aspect_ratio)

            self.assertEqual(src1.mfd.__class__, src2.mfd.__class__)
            if isinstance(src1.mfd, IncrementalMFD):
                self.assertEqual(src1.mfd.min_mag, src2.mfd.min_mag)
                self.assertEqual(src1.mfd.bin_width, src2.mfd.bin_width)
                self.assertEqual(src1.mfd.occur_rates, src2.mfd.occur_rates)
            elif isinstance(src1.mfd, TGRMFD):
                self.assertEqual(src1.mfd.a_val, src2.mfd.a_val)
                self.assertEqual(src1.mfd.b_val, src2.mfd.b_val)
                self.assertEqual(src1.mfd.min_mag, src2.mfd.min_mag)
                self.assertEqual(src1.mfd.max_mag, src2.mfd.max_mag)
            else:
                raise ValueError('MFD class %s not recognized' % mfd.__class__)

            # this covers both the case of point and area sources
            if isinstance(src1, PointSource):
                self.assert_coordinates_equal(
                    src1.geometry.wkt, src2.geometry.wkt
                )
                self.assertEqual(
                    src1.geometry.upper_seismo_depth,
                    src2.geometry.upper_seismo_depth
                )
                self.assertEqual(
                    src1.geometry.lower_seismo_depth,
                    src2.geometry.lower_seismo_depth
                )
                for np1, np2 in \
                    zip(src1.nodal_plane_dist, src2.nodal_plane_dist):
                    self.assertEqual(np1.probability, np2.probability)
                    self.assertEqual(np1.strike, np2.strike)
                    self.assertEqual(np1.dip, np2.dip)
                    self.assertEqual(np1.rake, np2.rake)
                for hd1, hd2 in zip(src1.hypo_depth_dist, src2.hypo_depth_dist):
                    self.assertEqual(hd1.probability, hd2.probability)
                    self.assertEqual(hd1.depth, hd2.depth)
            # first we check complex fault
            elif isinstance(src1, ComplexFaultSource) or \
                (isinstance(src1, CharacteristicSource) and
                 isinstance(src1.surface, ComplexFaultGeometry)):
                geo1 = src1.geometry if getattr(src1, 'geometry', None) \
                       else src1.surface
                geo2 = src2.geometry if getattr(src2, 'geometry', None) \
                       else src2.surface
                self.assert_coordinates_equal(
                    geo1.top_edge_wkt,
                    geo2.top_edge_wkt
                )
                self.assert_coordinates_equal(
                    geo1.bottom_edge_wkt,
                    geo2.bottom_edge_wkt
                )
                self.assertEqual(
                    len(geo1.int_edges),
                    len(geo2.int_edges)
                )
                for edge1, edge2 in \
                    zip(geo1.int_edges, geo2.int_edges):
                    self.assert_coordinates_equal(edge1, edge2)
                self.assertEqual(src1.rake, src2.rake)
            # then simple fault
            elif isinstance(src1, SimpleFaultSource) or\
                (isinstance(src1, CharacteristicSource) and
                 isinstance(src1.surface, SimpleFaultGeometry)):
                geo1 = src1.geometry if getattr(src1, 'geometry', None) \
                       else src1.surface
                geo2 = src2.geometry if getattr(src2, 'geometry', None) \
                       else src2.surface
                self.assert_coordinates_equal(
                    geo1.wkt,
                    geo2.wkt
                )
                self.assertEqual(geo1.dip, geo2.dip)
                self.assertEqual(
                    geo1.upper_seismo_depth, 
                    geo2.upper_seismo_depth
                )
                self.assertEqual(
                    geo1.lower_seismo_depth,
                    geo2.lower_seismo_depth
                )
                self.assertEqual(src1.rake, src2.rake)
            # then CharacteristicSource with list of PlanarSurface
            elif isinstance(src1, CharacteristicSource) and \
                isinstance(src1.surface, list):
                for surf1, surf2 in zip(src1.surface, src2.surface):
                    assert isinstance(surf1, PlanarSurface)
                    assert isinstance(surf2, PlanarSurface)
                    self.assertEqual(surf1.strike, surf2.strike)
                    self.assertEqual(surf1.dip, surf2.dip)
                    # top_left
                    self.assertEqual(surf1.top_left.longitude, surf2.top_left.longitude)
                    self.assertEqual(surf1.top_left.latitude, surf2.top_left.latitude)
                    self.assertEqual(surf1.top_left.depth, surf2.top_left.depth)
                    # top_right
                    self.assertEqual(surf1.top_right.longitude, surf2.top_right.longitude)
                    self.assertEqual(surf1.top_right.latitude, surf2.top_right.latitude)
                    self.assertEqual(surf1.top_right.depth, surf2.top_right.depth)
                    # bottom_left
                    self.assertEqual(surf1.bottom_left.longitude, surf2.bottom_left.longitude)
                    self.assertEqual(surf1.bottom_left.latitude, surf2.bottom_left.latitude)
                    self.assertEqual(surf1.bottom_left.depth, surf2.bottom_left.depth)
                    # bottom_right
                    self.assertEqual(surf1.bottom_right.longitude, surf2.bottom_right.longitude)
                    self.assertEqual(surf1.bottom_right.latitude, surf2.bottom_right.latitude)
                    self.assertEqual(surf1.bottom_right.depth, surf2.bottom_right.depth)
            else:
                raise ValueError(
                    'Source class %s not recognized' % src1.__class__
                )

    def test_complex_fault(self):
        # check that by converting original xml file to shapefile
        # and then converting back to xml, we get the same file
        nrml2shp('%ssource_model_cf.xml' % DATA_PATH)
        shp2nrml('%ssource_model_cf_complex' % DATA_PATH)

        self.assert_source_model_equals(
            '%ssource_model_cf.xml' % DATA_PATH,
            '%ssource_model_cf_complex.xml' % DATA_PATH
        )

    def test_simple_fault(self):
        # check that by converting original xml file to shapefile
        # and then converting back to xml, we get the same file
        nrml2shp('%ssource_model_sf.xml' % DATA_PATH)
        shp2nrml('%ssource_model_sf_simple' % DATA_PATH)

        self.assert_source_model_equals(
            '%ssource_model_sf.xml' % DATA_PATH,
            '%ssource_model_sf_simple.xml' % DATA_PATH
        )

    def test_area_source(self):
        # check that by converting original xml file to shapefile
        # and then converting back to xml, we get the same file
        nrml2shp('%ssource_model_as.xml' % DATA_PATH)
        shp2nrml('%ssource_model_as_area' % DATA_PATH)

        self.assert_source_model_equals(
            '%ssource_model_as.xml' % DATA_PATH,
            '%ssource_model_as_area.xml' % DATA_PATH
        )

    def test_point_source(self):
        # check that by converting original xml file to shapefile
        # and then converting back to xml, we get the same file
        nrml2shp('%ssource_model_ps.xml' % DATA_PATH)
        shp2nrml('%ssource_model_ps_point' % DATA_PATH)

        self.assert_source_model_equals(
            '%ssource_model_ps.xml' % DATA_PATH,
            '%ssource_model_ps_point.xml' % DATA_PATH
        )

    def test_planar_source(self):
        # check that by converting original xml file to shapefile
        # and then converting back to xml, we get the same file
        nrml2shp('%ssource_model_pl.xml' % DATA_PATH)
        shp2nrml('%ssource_model_pl_planar' % DATA_PATH)

        self.assert_source_model_equals(
            '%ssource_model_pl.xml' % DATA_PATH,
            '%ssource_model_pl_planar.xml' % DATA_PATH
        )