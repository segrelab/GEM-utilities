import unittest
import os

from cobra.io import read_sbml_model

from gem_utilities import maintenance

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'test_files')


class TestMaintenance(unittest.TestCase):

    def test_is_maintenance_reaction(self):
        """Test is_maintenance_reaction function."""
        # Read in the E coli core model
        model = read_sbml_model(os.path.join(TESTFILE_DIR,
                                             'ecc_missing_formulas.xml'))

        # Get the maintenance reaction
        maintenance_rxn = model.reactions.get_by_id('ATPM')

        # Assert that the maintenance reaction is a maintenance reaction
        self.assertTrue(maintenance.is_maintenance_reaction(model,
                                                            maintenance_rxn,
                                                            'BiGG'))

        # Get the biomass reaction
        biomass_rxn = model.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM')

        # Assert that the biomass reaction is not a maintenance reaction
        self.assertFalse(maintenance.is_maintenance_reaction(model,
                                                             biomass_rxn,
                                                             'BiGG'))

        # Check that you get an error if you use an invalid notation
        with self.assertRaises(ValueError):
            maintenance.is_maintenance_reaction(model, maintenance_rxn,
                                                'invalid notation')
            
        # Check that you get no maintenance reaction if you use the wrong
        # notation
        self.assertFalse(maintenance.is_maintenance_reaction(model,
                                                             maintenance_rxn,
                                                             'ModelSEED'))

        # TODO: Check that you get no maintenance reaction if there is none


if __name__ == '__main__':
    unittest.main()
