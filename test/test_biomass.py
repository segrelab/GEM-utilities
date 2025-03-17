import unittest

from cobra import Metabolite, Model, Reaction

from gem_utilities.biomass import unlump_biomass


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


if __name__ == "__main__":
    unittest.main()
