import os
import unittest
from io import StringIO

import cobra
import pandas as pd
from cobra import Metabolite, Model, Reaction

from gem_utilities.biomass import (
    calculate_biomass_weight,
    check_biomass_producibility,
    unlump_biomass,
)

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_files")


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

        # Create a biomass reaction that consumes the DNA
        self.biomass = Reaction("Biomass")
        self.biomass.add_metabolites({dna: -1.0})

        # Add the reaction to the model
        self.model.add_reactions([dna_synth, self.biomass])

    def test_unlump_biomass(self):
        # Expected output
        expected_unlumped = {"dAMP_c0": -1, "dCMP_c0": -1, "dGMP_c0": -1, "dTMP_c0": -1}

        # Run the unlump_biomass function
        result = unlump_biomass(
            self.biomass.metabolites, self.model, lumped_metabolites=["cpd11461_c0"]
        )

        # Check if the result matches the expected output
        # Use the metabolite IDs as keys to avoid issues with the metabolite
        # object location changing
        self.assertEqual({m.id: s for m, s in result.items()}, expected_unlumped)


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
            os.path.join(
                TESTFILE_DIR, "e_coli_core_biomass_producibility_all_sinks.csv"
            )
        )
        result_csv = pd.read_csv("e_coli_core_biomass_producibility_all_sinks.csv")
        pd.testing.assert_frame_equal(expected_csv, result_csv)

    def tearDown(self):
        # Clean up the test files
        if os.path.exists("e_coli_core_biomass_producibility_heatmap_all_sinks.png"):
            os.remove("e_coli_core_biomass_producibility_heatmap_all_sinks.png")
        if os.path.exists("e_coli_core_biomass_producibility_all_sinks.csv"):
            os.remove("e_coli_core_biomass_producibility_all_sinks.csv")


class TestBiomassWeight(unittest.TestCase):
    def setUp(self):
        # Create metabolites each with a formula weight of ~1.0
        # Because you cannot set the formula weight directly (it is calculated
        # from the formula), we create the metabolites with a formula of "H"
        # to give them a weight of 1.00794.
        met_a = Metabolite("A", name="Metabolite A", compartment="c0", formula="H")
        met_b = Metabolite("B", name="Metabolite B", compartment="c0", formula="H")
        met_c = Metabolite("C", name="Metabolite C", compartment="c0", formula="H")

        # Make a metabolite with no formula weight
        met_d = Metabolite("D", name="Metabolite D", compartment="c0")

        # Make ATP, ADP, and Pi metabolites with real formulas/weights
        atp = Metabolite("ATP", name="ATP", compartment="c0", formula="C10H12N5O13P3")
        adp = Metabolite("ADP", name="ADP", compartment="c0", formula="C10H12N5O10P2")
        pi = Metabolite("Pi", name="Phosphate", compartment="c0", formula="HO4P")
        h2o = Metabolite("H2O", name="Water", compartment="c0", formula="H2O")
        h = Metabolite("H", name="H+ Ion", compartment="c0", formula="H")

        # Create a biomass reaction that consumes the metabolites
        biomass_no_GAM = Reaction("Biomass_no_GAM")
        biomass_no_GAM.add_metabolites(
            {
                met_a: -1.0,
                met_b: -1.0,
                met_c: -1.0,
            }
        )

        # Create a biomass reaction that includes a metabolite with no formula weight
        biomass_no_GAM_missing_weight = Reaction("Biomass_no_GAM_missing_weight")
        biomass_no_GAM_missing_weight.add_metabolites(
            {
                met_a: -1.0,
                met_b: -1.0,
                met_c: -1.0,
                met_d: -1.0,  # This metabolite has no formula weight
            }
        )

        # Create a biomass reaction that includes the growth associated maintenance
        biomass_with_GAM = Reaction("Biomass_with_GAM")
        biomass_with_GAM.add_metabolites(
            {
                met_a: -1.0,
                met_b: -1.0,
                met_c: -1.0,
                atp: -1.0,
                h2o: -1.0,
                adp: 1.0,
                pi: 1.0,
                h: 1.0,
            }
        )

        # Create models for the two biomass reactions
        self.model_no_GAM = Model("test_model_no_GAM")
        self.model_no_GAM.add_reactions([biomass_no_GAM])
        self.model_with_GAM = Model("test_model_with_GAM")
        self.model_with_GAM.add_reactions([biomass_with_GAM])
        self.model_missing_weight = Model("test_model_missing_weight")
        self.model_missing_weight.add_reactions([biomass_no_GAM_missing_weight])

    def test_calculate_biomass_weight_toy_models(self):
        """Test the calculate_biomass_weight function using toy models"""
        exp_weight = 3.0

        # Calculate the biomass weight without GAM
        weight_no_GAM = calculate_biomass_weight(
            self.model_no_GAM, "Biomass_no_GAM", lumped_biomass_components=None
        )
        self.assertAlmostEqual(weight_no_GAM, exp_weight, places=1)

        # Calculate the biomass weight with GAM
        weight_with_GAM = calculate_biomass_weight(
            self.model_with_GAM, "Biomass_with_GAM", lumped_biomass_components=None
        )
        self.assertAlmostEqual(weight_with_GAM, exp_weight, places=1)

        # Check that the function raises an error for a biomass reaction with a
        # metabolite that has no formula weight
        with self.assertRaises(ValueError):
            calculate_biomass_weight(
                self.model_missing_weight,
                "Biomass_no_GAM_missing_weight",
                lumped_biomass_components=None,
            )

    def test_calculate_biomass_weight_e_coli(self):
        """Test the calculate_biomass_weight function using the E. coli core model"""
        # Set what the expected outputs are
        exp_weight = 999.0
        exp_table = pd.read_csv(
            os.path.join(TESTFILE_DIR, "iML1515_biomass_weight_work_table.csv")
        )

        # Load the full E. coli from COBRApy
        model = cobra.io.load_model("iML1515")

        # Calculate the biomass weight and save the work table
        weight = calculate_biomass_weight(
            model,
            "BIOMASS_Ec_iML1515_WT_75p37M",
            lumped_biomass_components=None,
            save_work_table=True,
            out_dir=".",
        )

        # Compare the returned weight with the expected value
        self.assertAlmostEqual(weight, exp_weight, places=0)
        # Check that the work table matches the expected output
        work_table = pd.read_csv("iML1515_biomass_weight_work_table.csv")
        pd.testing.assert_frame_equal(work_table, exp_table)

    def tearDown(self):
        # Clean up the test files
        if os.path.exists("iML1515_biomass_weight_work_table.csv"):
            os.remove("iML1515_biomass_weight_work_table.csv")


if __name__ == "__main__":
    unittest.main()
