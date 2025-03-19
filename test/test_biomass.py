import os
import unittest
from io import StringIO

import cobra
import pandas as pd
from cobra import Metabolite, Model, Reaction

from gem_utilities.biomass import check_biomass_producibility, unlump_biomass


class TestUnlumpBiomass(unittest.TestCase):
    def setUp(self):
        # Create a simple test model
        self.model = Model("test_model")

        # Create metabolites
        dna = Metabolite("cpd11461_c0", name="DNA", compartment="c0")
        damp = Metabolite("dAMP_c0", name="dAMP", compartment="c0")
        dcmp = Metabolite("dCMP_c0", name="dCMP", compartment="c0")
        dgmp = Metabolite("dGMP_c0", name="dGMP", compartment="c0")
        dtmp = Metabolite("dTMP_c0", name="dTMP", compartment="c0")

        # Create a reaction that synthesizes DNA from its components
        dna_synth = Reaction("DNA_synth")
        dna_synth.add_metabolites(
            {damp: -1.0, dcmp: -1.0, dgmp: -1.0, dtmp: -1.0, dna: 1.0}
        )

        # Add the reaction to the model
        self.model.add_reactions([dna_synth])

        # Biomass compounds list
        self.biomass_compounds = ["cpd11461_c0"]

    def test_unlump_biomass(self):
        # Expected output
        expected_unlumped = ["dAMP_c0", "dCMP_c0", "dGMP_c0", "dTMP_c0"]

        # Run the unlump_biomass function
        result = unlump_biomass(
            self.model, self.biomass_compounds, lumped_metabolites=["cpd11461_c0"]
        )

        # Check if the result matches the expected output
        self.assertEqual(sorted(result), sorted(expected_unlumped))


class TestBiomassComponentProducibility(unittest.TestCase):
    def setUp(self):
        # Make a simple growth phenotypes dataframe
        phenotypes_csv = """
minimal_media,c_source,met_id,growth
minimal,control,,No
minimal,glucose,glc__D,Yes
minimal,acetate,ac,Yes
        """

        # Use StringIO to simulate a file
        self.phenotypes_df = pd.read_csv(StringIO(phenotypes_csv.strip()))

        # Define the minimal media
        self.media_definitions = {
            "minimal": {
                "EX_co2_e": 1000.0,
                "EX_h_e": 1000.0,
                "EX_h2o_e": 1000.0,
                "EX_nh4_e": 1000.0,
                "EX_o2_e": 1000.0,
                "EX_pi_e": 1000.0,
            }
        }

    def test_biomass_producibility(self):
        """Test the test_biomass_producibility function using the E. coli core
        model and a simple media definition"""
        model = cobra.io.load_model("textbook")
        check_biomass_producibility(
            model,
            self.phenotypes_df,
            self.media_definitions,
            biomass_rxn="Biomass_Ecoli_core",
            external_compartment="e",
            lumped_biomass_components=None,
        )

        # Compare the resulting CSV file with the expected output
        expected_csv = pd.read_csv(
            "test/test_files/e_coli_core_biomass_producibility.csv"
        )
        result_csv = pd.read_csv("e_coli_core_biomass_producibility.csv")
        pd.testing.assert_frame_equal(expected_csv, result_csv)

    def tearDown(self):
        # Clean up the test files
        if os.path.exists("e_coli_core_biomass_producibility_heatmap.png"):
            os.remove("e_coli_core_biomass_producibility_heatmap.png")
        if os.path.exists("e_coli_core_biomass_producibility.csv"):
            os.remove("e_coli_core_biomass_producibility.csv")


if __name__ == "__main__":
    unittest.main()
