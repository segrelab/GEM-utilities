import unittest
import tempfile
import filecmp
import os

import cobra

from gem_utils import missing_formulas

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files')


class TestCopyFormulas(unittest.TestCase):

    def test_copy_formulas(self):
        """Test copy_formulas function."""
        # Read in the E coli core model where the external acetate (M_ac_e) is missing a chemical formula
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, 'e_coli_core_missing_formulas.xml'))

        # Get the chemical formula for external acetate
        ac_e_formula = model.metabolites.get_by_id('ac_e').formula

        # Assert that the chemical formula is missing (by checking that it is False)
        self.assertFalse(ac_e_formula)

        # Copy the chemical formula from the internal acetate (M_ac_c) to the external acetate (M_ac_e)
        missing_formulas.copy_formulas(model, 'ac_c', 'ac_e')

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
        assert filecmp.cmp(tmp_out, os.path.join(TESTFILE_DIR, 'e_coli_core_missing_formulas_fixed.xml'))

    def test_pro_model(self):
        """Test fix_names_w_compartment_suffix function on a PRO model."""
        # Read in the PRO model
        model = cobra.io.read_sbml_model('../../GEM-repos/Prochlorococcus_Model/Model_files/iSO595v7.xml')

        # Find all the metabolites with missing formulas
        mets_wo_formula = [met for met in model.metabolites if not met.formula]
        # Find a match and copy the formula
        for met in mets_wo_formula:
            match = missing_formulas.find_matching_metabolite(model, met)
            if match:
                missing_formulas.copy_formulas(model, match.id, met.id)

        # Write the model to a temporary file
        tmp_out = tempfile.mkstemp(suffix='.xml')[1]
        cobra.io.write_sbml_model(model, tmp_out)

        # Compare the model with the expected model
        assert filecmp.cmp(tmp_out, os.path.join(TESTFILE_DIR, 'pro_model_fixed.xml'))


if __name__ == '__main__':
    unittest.main()
