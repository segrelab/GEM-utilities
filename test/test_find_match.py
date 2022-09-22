import unittest
import tempfile
import filecmp

import cobra

from gem_utils import missing_formulas


class TestFindMatch(unittest.TestCase):

    def test_find_matching_metabolite(self):
        """Test find_matching_metabolite function."""
        # Read in the E coli core model where the external acetate (M_ac_e) is missing a chemical formula
        model = cobra.io.read_sbml_model('test_files/e_coli_core_missing_formulas.xml')

        # Get the external acetate metabolite
        ac_e = model.metabolites.get_by_id('ac_e')

        # Find a metabolite with the same name in a different compartment
        matching_met = missing_formulas.find_matching_metabolite(model, ac_e)

        # Check that the matching metabolite is the internal acetate metabolite
        self.assertEqual(matching_met.id, 'ac_c')


if __name__ == '__main__':
    unittest.main()
