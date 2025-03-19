import os
import unittest

import cobra
from cobra import Metabolite, Model, Reaction

from gem_utilities.annotation import add_kos_to_model, get_ko_for_kegg_reaction


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
