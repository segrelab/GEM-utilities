import unittest
import tempfile
import filecmp
import os

from cobra.io import read_sbml_model, load_model

from gem_utilities import names

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'test_files')


class TestNames(unittest.TestCase):
    def test_find_names_w_compartment_suffix(self):
        """Test find_names_w_compartment_suffix function."""
        # Read in the E coli core model
        model = read_sbml_model(os.path.join(TESTFILE_DIR,
                                'ecc_compartment_names.xml'))

        # Find metabolites with a compartment suffix
        names_w_compartment_suffix = names.find_names_w_compartment_suffix(
            model)

        # Assert that the correct metabolites are found
        self.assertEqual(len(names_w_compartment_suffix), 1)
        self.assertEqual(names_w_compartment_suffix['gln__L_e'],
                         'L-Glutamine[e]')

    def test_trim_name(self):
        """Test trim_name function."""
        # Trim the compartment suffix
        trimmed_name = names.trim_name('L-Glutamine[e]')

        # Assert that the compartment suffix is trimmed
        self.assertEqual(trimmed_name, 'L-Glutamine')

    def test_fix_names_w_compartment_suffix(self):
        """Test fix_names_w_compartment_suffix function."""
        # Read in the E coli core model
        model = read_sbml_model(os.path.join(TESTFILE_DIR,
                                             'ecc_compartment_names.xml'))

        # Assert that the metabolite name for L-glutamine includes the
        # compartment suffix
        self.assertEqual(model.metabolites.get_by_id('gln__L_e').name,
                         'L-Glutamine[e]')

        # Fix the metabolite names
        names.fix_names_w_compartment_suffix(model)

        # Assert that the metabolite name for L-glutamine no longer includes
        # the compartment suffix
        self.assertEqual(model.metabolites.get_by_id('gln__L_e').name,
                         'L-Glutamine')


class TestBuildReactionString(unittest.TestCase):
    def test_build_reaction_string(self):
        """Test build_reaction_string function."""
        # Read in the E coli core model
        model = load_model("textbook")

        # Find a reaction
        reaction = model.reactions.get_by_id('EX_glc__D_e')

        # Build the reaction string
        reaction_string = names.build_reaction_string(reaction)

        # Assert that the reaction string is correct
        self.assertEqual(reaction_string, '1.0 D-Glucose <--> ')

        # Use a reaction that is irreversible
        reaction = model.reactions.get_by_id('ATPM')
        reaction_string = names.build_reaction_string(reaction)
        self.assertEqual(reaction_string, '1.0 ATP + 1.0 H2O --> 1.0 ADP + 1.0 H+ + 1.0 Phosphate')