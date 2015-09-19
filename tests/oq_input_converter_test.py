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
import subprocess
import unittest
import numpy as np

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
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "sample_site_model.xml")
        exit_code = subprocess.call(["python",
                                     self.prog,
                                     "--input-xml-file",
                                     self.input_file,
                                     "--output-csv-file",
                                     "dummy_site.csv"])
        self.assertEqual(exit_code, 0)
        # cleanup
        subprocess.call(["rm", "dummy_site.csv"])

    def test_csv_to_xml(self):
        """

        """
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "sample_site_model.csv")
        exit_code = subprocess.call(["python",
                                     self.prog,
                                     "--input-csv-file",
                                     self.input_file,
                                     "--output-xml-file",
                                     "dummy_site.xml"])
        self.assertEqual(exit_code, 0)
        # cleanup
        subprocess.call(["rm", "dummy_site.xml"])

    def test_roundtrip(self):
        """

        """
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "sample_site_model.xml")
        exit_code = subprocess.call(["python",
                                     self.prog,
                                     "--input-xml-file",
                                     self.input_file,
                                     "--output-csv-file",
                                     "dummy_site.csv"])
        self.assertEqual(exit_code, 0)
        
        exit_code = subprocess.call(["python",
                                     self.prog,
                                     "--input-csv-file",
                                     "dummy_site.csv",
                                     "--output-xml-file",
                                     "dummy_site2.xml"])
        self.assertEqual(exit_code, 0)
        subprocess.call(["rm", "dummy_site.csv"])
        subprocess.call(["rm", "dummy_site2.xml"])


class SourceModelShapefileConverterTestCase

        

