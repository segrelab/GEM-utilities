import unittest
import os

import cobra

from gem_utilities import missing_formulas

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files')


class TestFindMatch(unittest.TestCase):

    def test_find_matching_metabolite(self):
        """Test find_matching_metabolite function."""
        # Read in the E coli core model where the external acetate (M_ac_e) is missing a chemical formula
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, 'e_coli_core_missing_formulas.xml'))

        # Get the external acetate metabolite
        ac_e = model.metabolites.get_by_id('ac_e')

        # Find a metabolite with the same name in a different compartment
        matching_met = missing_formulas.find_matching_metabolite(model, ac_e)

        # Check that the matching metabolite is the internal acetate metabolite
        self.assertEqual(matching_met.id, 'ac_c')


if __name__ == '__main__':
    unittest.main()
