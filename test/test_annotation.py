import os
import unittest

import cobra
from cobra import Metabolite, Model, Reaction

from gem_utilities.annotation import add_kos_to_model, get_ko_for_kegg_reaction

TESTFILE_DIR = os.path.join(os.path.dirname(__file__), "test_files")


class TestAddKOs(unittest.TestCase):
    def test_get_ko_for_kegg_reaction(self):
        """Test the get_ko_for_kegg_reaction function using a couple of known
        KEGG reactions."""
        # Test the function on ATP:D-glucose 6-phosphotransferase
        # The entry on the website is: https://www.kegg.jp/entry/R00299
        kos = get_ko_for_kegg_reaction("R00299")
        # Sort the output to avoid order issues
        kos.sort()

        # The reaction R00299 is known to be associated with the following KOs
        # K00844, K00845, K12407, K25026
        self.assertEqual(kos, ["K00844", "K00845", "K12407", "K25026"])

    def test_add_kos_to_model(self):
        """Test the add_kos_to_model function using a simple model, with
        just one reaction."""
        # Load the model file
        model = cobra.io.read_sbml_model(
            os.path.join(TESTFILE_DIR, "ecc_all_formulas.xml")
        )

        # Load the expected model file
        expected_model = cobra.io.read_sbml_model(
            os.path.join(TESTFILE_DIR, "ecc_all_formulas_with_kos.xml")
        )

        # Run the function
        new_model = add_kos_to_model(model)

        # Check that every reaction that has a KEGG annotation also has a KO annotation
        for reaction in new_model.reactions:
            if "kegg.reaction" in reaction.annotation:
                self.assertTrue("kegg.orthology" in reaction.annotation, msg=f'Reaction {reaction.id} does not have a KO annotation')

        # Check that every ko annotation is correct
        for reaction in new_model.reactions:
            if "ko" in reaction.annotation:
                self.assertEqual(
                    reaction.annotation["kegg.orthology"],
                    expected_model.reactions.get_by_id(reaction.id).annotation[
                        "kegg.orthology"
                    ],
                    msg=f'KO annotation for reaction {reaction.id} is incorrect')
