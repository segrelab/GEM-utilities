import unittest
import os
import json
import cobra

from gem_utilities import escher

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'test_files')


class TestEscher(unittest.TestCase):

    def test_create_comparison_data_file(self):
        """Test create_comparison_data_file function."""
        # Test 1: Give two minimal models and check that the reaction
        # data dictionary is correct
        # TODO: Make more minimal models, right now, these are pretty big
        # Load in the models
        model1 = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR,
                                                       'mit1002_core_model.xml'))
        model2 = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR,
                                                       '102897.3_CoreModel.xml'))
        # Get the reaction data dictionary
        rxn_data = escher.create_comparison_data_file(model1, model2)
        # Load in the saved reaction data dictionary
        with open(os.path.join(TESTFILE_DIR, 'reaction_data.json')) as f:
            saved_rxn_data = json.load(f)
        # Check that the reaction data dictionary is correct
        self.assertEqual(rxn_data, saved_rxn_data)

        # Test 2: Give something that isn't a COBRA model and check that
        # you get an error
        with self.assertRaises(TypeError):
            escher.create_comparison_data_file('not a COBRA model')

        # TODO: Test 3: Give models with different compartment
        # nomenclatures and check that you get a warning