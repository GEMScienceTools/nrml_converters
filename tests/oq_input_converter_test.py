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
import os
import unittest
from tests import gem_run_script, gem_rmtree, gem_unlink

BASEPATH = os.path.join(os.path.dirname(__file__), os.path.pardir, "oq_input")


class SiteModelConverterTestCase(unittest.TestCase):
    """
    Round trip test for the Site Model converter
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "site_model_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "sample_site_model.xml")

    def test_xml_to_csv(self):
        """

        """
        gem_unlink("dummy_site.csv")

        input_file = os.path.join(os.path.dirname(__file__),
                                  "..", "sample_data",
                                  "sample_site_model.xml")
        gem_run_script(self.prog, ["--input-xml-file",
                                   input_file,
                                   "--output-csv-file",
                                   "dummy_site.csv"])
        # cleanup
        gem_unlink("dummy_site.csv", True)

    def test_csv_to_xml(self):
        """

        """
        gem_unlink("dummy_site.xml")
        input_file = os.path.join(os.path.dirname(__file__),
                                  "..", "sample_data",
                                  "sample_site_model.csv")
        gem_run_script(self.prog, ["--input-csv-file",
                                   input_file,
                                   "--output-xml-file",
                                   "dummy_site.xml"])
        # cleanup
        gem_unlink("dummy_site.xml", True)

    def test_roundtrip(self):
        """

        """
        gem_unlink("dummy_site.csv")
        gem_unlink("dummy_site2.xml")
        input_file = os.path.join(os.path.dirname(__file__),
                                  "..", "sample_data",
                                  "sample_site_model.xml")
        gem_run_script(self.prog, ["--input-xml-file",
                                   input_file,
                                   "--output-csv-file",
                                   "dummy_site.csv"])

        gem_run_script(self.prog, ["--input-csv-file",
                                   "dummy_site.csv",
                                   "--output-xml-file",
                                   "dummy_site2.xml"])

        gem_unlink("dummy_site.csv", True)
        gem_unlink("dummy_site2.xml", True)


class SourceModelShapefileConverterTestCase(unittest.TestCase):
    """
    Simple conversion test for the Source Model to shapefile converter
    - more tests will follow
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "source_model_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "sample_source_model.xml")

    def test_xml_to_shp_no_validation(self):
        """
        Tests the conversion to shapefile - without validation
        """
        #output_dir = os.path.join(os.path.dirname(__file__), "outputs")
        gem_rmtree("source_shp")
        os.mkdir("source_shp")

        gem_run_script(self.prog, [
            "--input-nrml-file",
            self.input_file,
            "--output-file",
            os.path.join("source_shp", "src_model_files")])
        # Cleanup
        gem_rmtree("source_shp")

    def test_xml_to_shp_with_validation(self):
        """
        Tests the conversion to shapefile - with validation
        """
        #output_dir = os.path.join(os.path.dirname(__file__), "outputs")
        gem_rmtree("source_shp")
        os.mkdir("source_shp")
        gem_run_script(self.prog, [
            "--input-nrml-file",
            self.input_file,
            "--output-file",
            os.path.join("source_shp", "src_model_files"),
            "--validate"])
        # Cleanup
        gem_rmtree("source_shp", True)

    def test_shp_to_xml_no_validation(self):
        """
        """
        gem_unlink("dummy1.xml")
        gem_rmtree("source_shp")
        os.mkdir("source_shp")

        # Build shapefiles
        gem_run_script(self.prog, [
            "--input-nrml-file",
            self.input_file,
            "--output-file",
            os.path.join("source_shp", "src_model_files")])

        # Run test for all shapefiles
        point_file = os.path.join("source_shp", "src_model_files" + "_point")
        area_file = os.path.join("source_shp", "src_model_files" + "_area")
        simple_file = os.path.join("source_shp", "src_model_files" + "_simple")
        complex_file = os.path.join("source_shp",
                                    "src_model_files" + "_complex")
        gem_run_script(self.prog, [
            "--input-shp-files",
            point_file,
            area_file,
            simple_file,
            complex_file,
            "--output-file",
            "dummy1.xml"])

        # cleanup
        gem_rmtree("source_shp", True)
        gem_unlink("dummy1.xml", True)

    def test_shp_to_xml_validation(self):
        """
        """
        # Build shapefiles
        gem_unlink("dummy1.xml")
        gem_rmtree("source_shp")
        os.mkdir("source_shp")

        gem_run_script(self.prog, [
            "--input-nrml-file",
            self.input_file,
            "--output-file",
            os.path.join("source_shp", "src_model_files"),
            "--validate"])

        # Run test for all shapefiles
        point_file = os.path.join("source_shp", "src_model_files" + "_point")
        area_file = os.path.join("source_shp", "src_model_files" + "_area")
        simple_file = os.path.join("source_shp", "src_model_files" + "_simple")
        complex_file = os.path.join("source_shp",
                                    "src_model_files" + "_complex")
        gem_run_script(self.prog, ["--input-shp-files",
                                   point_file,
                                   area_file,
                                   simple_file,
                                   complex_file,
                                   "--output-file",
                                   "dummy1.xml",
                                   "--validate"])
        # cleanup
        gem_unlink("dummy1.xml", True)
        gem_rmtree("source_shp", True)
