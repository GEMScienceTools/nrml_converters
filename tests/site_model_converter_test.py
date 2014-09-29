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
import numpy as np
from shapely import wkt
from oq_input.site_model_converter import csv_to_xml, xml_to_csv 
from openquake.nrmllib.hazard.parsers import SiteModelParser


class TestSiteModelCsv2Nrml(unittest.TestCase):
    """
    Tests the correct parsing of a csv site model file to nrml
    """
    def setUp(self):
        """

        """
        self.input_csv = os.path.join(".", "tests/data/sample_site_model.csv")
        self.output_xml = "temp_site_xml.xml"
        self.reference_xml = os.path.join(".",
                                          "tests/data/sample_site_model.xml")


    def _compare_site_models(self, model1, model2):
        """
        Utility function to test two site model classes
        """
        self.assertEqual(len(model1), len(model2))
        for iloc in range(0, len(model1)):
            site1 = model1[iloc]
            site2 = model2[iloc]
            locn1 = wkt.loads(site1.wkt)
            locn2 = wkt.loads(site2.wkt)
            self.assertAlmostEqual(locn1.x, locn2.x)
            self.assertAlmostEqual(locn1.y, locn2.y)
            self.assertAlmostEqual(site1.vs30, site2.vs30)
            self.assertAlmostEqual(site1.z1pt0, site2.z1pt0)
            self.assertAlmostEqual(site1.z2pt5, site2.z2pt5)
            self.assertEqual(site1.vs30_type, site2.vs30_type)

    def test_csv_to_xml(self):
        """
        Tests the conversion of csv to xml
        """
        # Convert csv file
        csv_to_xml(self.input_csv, self.output_xml)
        # Load in both reference and exported xml files
        parser1 = SiteModelParser(self.output_xml)
        model1 = list(parser1.parse())
        parser2 = SiteModelParser(self.reference_xml)
        model2 = list(parser2.parse())
        # Test comparison
        self._compare_site_models(model1, model2)
        # Cleanup
        os.remove(self.output_xml)
