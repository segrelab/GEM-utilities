import unittest
import tempfile
import filecmp

import cobra

from gem_utils import copy_formulas


class TestCopyFormulas(unittest.TestCase):

    def test_copy_formulas(self):
        """Test copy_formulas function."""
        # Read in the E coli core model where the external acetate (M_ac_e) is missing a chemical formula
        model = cobra.io.read_sbml_model('test_files/e_coli_core_missing_formulas.xml')

        # Get the chemical formula for external acetate
        ac_e_formula = model.metabolites.get_by_id('ac_e').formula

        # Assert that the chemical formula is missing (by checking that it is False)
        self.assertFalse(ac_e_formula)

        # Copy the chemical formula from the internal acetate (M_ac_c) to the external acetate (M_ac_e)
        copy_formulas.copy_formulas(model, 'ac_c', 'ac_e')

        # Now get the chemical formula for external acetate
        ac_e_formula = model.metabolites.get_by_id('ac_e').formula

        # Assert that the chemical formula is now present (by checking that it is True)
        self.assertTrue(ac_e_formula)

        # Assert that the chemical formula is correct
        self.assertEqual(ac_e_formula, 'C2H3O2')

        # Write the model to a temporary file
        tmp_out = tempfile.mkstemp(suffix='.xml')[1]
        cobra.io.write_sbml_model(model, tmp_out)

        # Compare the model with the expected model
        assert filecmp.cmp(tmp_out, 'test_files/e_coli_core_missing_formulas_fixed.xml')


if __name__ == '__main__':
    unittest.main()
