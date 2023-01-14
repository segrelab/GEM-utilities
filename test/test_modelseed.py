import unittest
import json

from gem_utils import modelseed


class TestModelSEED(unittest.TestCase):
    def test_generate_rxn(self):
        """Test generate_rxn function."""
        # Open the JSON file from the ModelSEED database and load it as a dictionary
        reactions_db = json.load(open('../../ModelSEEDDatabase/Biochemistry/reactions.json'))

        # Get the reaction data for the reaction with the EC number
        rxn = modelseed.generate_rxn('3.6.1.1', reactions_db)

        # Check that the reaction is as expected
        self.assertEqual(rxn.id, '3.6.1.1')
        self.assertEqual(rxn.name, 'diphosphate phosphohydrolase')
        self.assertEqual(rxn.lower_bound, -1000.0)
        self.assertEqual(rxn.upper_bound, 1000.0)
        # self.assertEqual(len(rxn.metabolites), 4)
        # self.assertEqual(rxn.metabolites['cpd00002_c0'], -1.0)
        # self.assertEqual(rxn.metabolites['cpd00001_c0'], -1.0)
        # self.assertEqual(rxn.metabolites['cpd00008_c0'], 1.0)
        # self.assertEqual(rxn.metabolites['cpd00009_c0'], 1.0)
