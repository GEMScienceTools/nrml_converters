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
Test suite for the NRML Converters for OpenQuake output

Each test tests only for successful execution, not for correctness or speed
"""

import os
import unittest
from tests import gem_run_script, gem_rmtree, gem_unlink

BASEPATH = os.path.join(os.path.dirname(__file__), os.path.pardir, "oq_output")


class HazardCurveConverterTestCase(unittest.TestCase):
    """
    Test case for the Hazard Curve Converter
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "hazard_curve_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "hazard_curves_short.xml")

    def test_wo_plotting(self):
        """
        Tests the execution without plotting
        """
        gem_unlink("dummy_hazard_curve.csv")

        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-file",
                                   "dummy_hazard_curve"])

        # Cleanup
        gem_unlink("dummy_hazard_curve.csv", True)

    def test_plotting(self):
        """
        Tests the execution with plotting
        """
        gem_rmtree("dummy_hazard_curve")
        gem_unlink("dummy_hazard_curve.csv")

        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-file",
                                   "dummy_hazard_curve",
                                   "--plot-curves",
                                   "True"])

        # Cleanup
        gem_rmtree("dummy_hazard_curve", True)
        gem_unlink("dummy_hazard_curve.csv", True)


class HazardMapConverterTestCase(unittest.TestCase):
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "hazard_map_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "hazard_map.xml")

    def test_hazard_map_converter(self):
        """
        Tests the hazard map converter
        """
        gem_unlink("dummy_hazard_map.csv")

        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-file",
                                   "dummy_hazard_map"])

        # Cleanup
        gem_unlink("dummy_hazard_map.csv", True)


class UHSConverterTestCase(unittest.TestCase):
    """
    Test case for the UHS converter
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "uhs_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "uniform_hazard_spectra_short.xml")

    def test_uhs_converter_wo_plotting(self):
        """
        Tests the execution without plotting
        """
        gem_unlink("dummy_uhs.csv")

        gem_run_script(self.prog, ["--input-file", self.input_file,
                                   "--output-file", "dummy_uhs"])

        # Cleanup
        gem_unlink("dummy_uhs.csv", True)

    def test_uhs_converter_plotting(self):
        """
        Tests the execution with plotting
        """
        gem_rmtree("dummy_uhs")
        gem_unlink("dummy_uhs.csv"
)
        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-file",
                                   "dummy_uhs",
                                   "--plot-spectra",
                                   "True"])

        # Cleanup
        gem_rmtree("dummy_uhs", True)
        gem_unlink("dummy_uhs.csv", True)


class ScenarioGMFConverterTestCase(unittest.TestCase):
    """
    Test case for the Scenario GMF converter
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "scenario_gmf_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "gmf_scenario.xml")

    def test_gmf_scenario_converter(self):
        """
        Tests scenario gmf execution
        """
        gem_rmtree("dummy_scenario_gmf")

        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-dir",
                                   "dummy_scenario_gmf"])

        # Cleanup
        gem_rmtree("dummy_scenario_gmf", True)


class EventSetGMFConverterTestCase(unittest.TestCase):
    """
    Test case for the Event Set GMF converter
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "gmfset_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "gmfs_event_based.xml")

    def test_gmf_set_converter(self):
        """
        Tests the execution of the gmf set converter
        """
        gem_rmtree("dummy_gmf_set")

        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-dir",
                                   "dummy_gmf_set"])

        # Cleanup
        gem_rmtree("dummy_gmf_set", True)


class EventSetConverterTestCase(unittest.TestCase):
    """
    Test case for the Event Set converter
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "eventset_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "event_set.xml")

    def test_gmf_set_converter(self):
        """
        Tests the execution of the gmf set converter
        """
        gem_rmtree("dummy_event_set")

        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-dir",
                                   "dummy_event_set"])

        # Cleanup
        gem_rmtree("dummy_event_set", True)


class DisaggregationConverterTestCase(unittest.TestCase):
    """
    Test case for the disaggregation converter
    """
    def setUp(self):
        self.prog = os.path.join(BASEPATH, "disaggregation_converter.py")
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       "..", "sample_data",
                                       "disaggregation.xml")

    def test_disaggregation_converter(self):
        """
        Tests the execution of the disaggregation converter
        """
        gem_rmtree("dummy_disag")

        gem_run_script(self.prog, ["--input-file",
                                   self.input_file,
                                   "--output-dir",
                                   "dummy_disag"])

        # Cleanup
        gem_rmtree("dummy_disag", True)
