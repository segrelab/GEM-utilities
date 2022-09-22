import unittest
import tempfile
import filecmp

import cobra

from gem_utils import names


class TestNames(unittest.TestCase):
    def test_find_names_w_compartment_suffix(self):
        """Test find_names_w_compartment_suffix function."""
        # Read in the E coli core model
        model = cobra.io.read_sbml_model('test_files/e_coli_core_compartment_names.xml')

        # Find metabolites with a compartment suffix
        names_w_compartment_suffix = names.find_names_w_compartment_suffix(model)

        # Assert that the correct metabolites are found
        self.assertEqual(len(names_w_compartment_suffix), 1)
        self.assertIn('L-Glutamine[e]', names_w_compartment_suffix)

    def test_trim_name(self):
        """Test trim_name function."""
        # Trim the compartment suffix
        trimmed_name = names.trim_name('L-Glutamine[e]')

        # Assert that the compartment suffix is trimmed
        self.assertEqual(trimmed_name, 'L-Glutamine')


