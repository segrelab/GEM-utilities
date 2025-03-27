import os
import unittest

import cobra

from gem_utilities.maps import Mapper, download_kegg_pathways

TESTFILE_DIR = os.path.join(os.path.dirname(__file__), "test_files")


class TestMapper(unittest.TestCase):
    def test_mapper(self):
        """Test the mapper function using a simple model
        Before you run this- you need a copy of the KEGG pathway files locally,
        which can be downloaded using the download_kegg_pathways function.
        e.g.
        download_kegg_pathways("/path/to/kegg_data")
        """
        # Load the model file (must laready have KOs added)
        model = cobra.io.read_sbml_model(
            os.path.join(TESTFILE_DIR, "ecc_all_formulas_with_kos.xml")
        )

        # Ensure KO annotations are in the model
        for reaction in model.reactions:
            if "kegg.orthology" not in reaction.annotation:
                reaction.annotation["kegg.orthology"] = []  # Initialize if needed

        # Create the mapper
        mapper = Mapper(
            model, kegg_dir="/Users/helenscott/Documents/PhD/Segre-lab/kegg_data/kgml"
        )

        # Draw maps for a single model
        mapper.map_model_kos(
            output_dir="output/pathway_maps",
            pathway_numbers=["00300"],
            color_hexcode="#2ca02c",  # Use a green color
        )

    def test_compare_models(self):
        # For comparing multiple models
        model = cobra.io.read_sbml_model("my_model1.xml")
        model2 = cobra.io.read_sbml_model("my_model2.xml")
        models = {"Model1": model, "Model2": model2}

        # Create the mapper
        mapper = Mapper(model, kegg_dir="/path/to/kegg_data")

        mapper.map_models_kos(
            models=models,
            output_dir="output/comparison_maps",
            draw_individual_files=True,
            draw_grid=True,
        )
