import unittest
import tempfile
import filecmp
import os

import cobra

from gem_utils import names

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files')


class TestNames(unittest.TestCase):
    def test_find_names_w_compartment_suffix(self):
        """Test find_names_w_compartment_suffix function."""
        # Read in the E coli core model
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, 'e_coli_core_compartment_names.xml'))

        # Find metabolites with a compartment suffix
        names_w_compartment_suffix = names.find_names_w_compartment_suffix(model)

        # Assert that the correct metabolites are found
        self.assertEqual(len(names_w_compartment_suffix), 1)
        self.assertEqual(names_w_compartment_suffix['gln__L_e'], 'L-Glutamine[e]')

    def test_trim_name(self):
        """Test trim_name function."""
        # Trim the compartment suffix
        trimmed_name = names.trim_name('L-Glutamine[e]')

        # Assert that the compartment suffix is trimmed
        self.assertEqual(trimmed_name, 'L-Glutamine')

    def test_fix_names_w_compartment_suffix(self):
        """Test fix_names_w_compartment_suffix function."""
        # Read in the E coli core model
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, 'e_coli_core_compartment_names.xml'))

        # Fix the metabolite names
        names.fix_names_w_compartment_suffix(model)

        # Write the model to a temporary file
        tmp_out = tempfile.mkstemp(suffix='.xml')[1]
        cobra.io.write_sbml_model(model, tmp_out)

        # Compare the model with the expected model
        assert filecmp.cmp(tmp_out, os.path.join(TESTFILE_DIR, 'e_coli_core_missing_formulas_fixed.xml'))
