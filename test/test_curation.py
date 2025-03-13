import filecmp
import json
import os
import tempfile
import unittest

import cobra

from gem_utilities import curation

TESTFILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_files")


class TestCuration(unittest.TestCase):

    def test_parse_aliases(self):
        aliases = [
            "Name: ATP; Adenosine 5'-triphosphate; adenosine-5'-triphosphate; adenosine-triphosphate; adenylpyrophosphate",
            "AraCyc: ATP",
            "BiGG: atp; fake-id",
            "BrachyCyc: ATP",
            "KEGG: C00002",
            "MetaCyc: ATP",
        ]
        expected_output = {
            "aracyc": ["ATP"],
            "bigg": ["atp", "fake-id"],
            "brachycyc": ["ATP"],
            "kegg": ["C00002"],
            "metacyc": ["ATP"],
        }
        result = curation.parse_aliases(aliases)
        self.assertEqual(result, expected_output)

    def test_create_cobra_metabolite(self):
        """Test create_cobra_metabolite function, by creating a metabolite
        from the ModelSEED database."""

        # Load the (subset of the) modelseed databases
        met_db = json.load(
            open(os.path.join(TESTFILE_DIR, "mini_modelseed_db_compounds.json"))
        )

        # Convert the ModelSEED databases to dictionaries for easy searching
        template_met_db = {met["id"]: met for met in met_db if not met["is_obsolete"]}

        # Create a new metabolite
        new_met = curation.create_cobra_metabolite(template_met_db, "cpd00002", "c0")
        # Check that the metabolite was created
        self.assertIsInstance(new_met, cobra.Metabolite)
        # Check that the metabolite has the correct information
        self.assertEqual(new_met.id, "cpd00002_c0")
        self.assertEqual(new_met.name, "ATP")
        self.assertEqual(new_met.formula, "C10H13N5O13P3")
        self.assertEqual(new_met.charge, -3)
        self.assertEqual(
            new_met.annotation,
            {
                "aracyc": ["ATP"],
                "bigg": ["atp"],
                "brachycyc": ["ATP"],
                "kegg": ["C00002"],
                "metacyc": ["ATP"],
            },
        )
        self.assertEqual(new_met.notes, {})

        # Add it to an empty model
        model = cobra.Model()
        model.add_metabolites([new_met])

        # Save the model to a temporary file
        with tempfile.NamedTemporaryFile(suffix=".xml") as temp_file:
            cobra.io.write_sbml_model(model, temp_file.name)
            # Check that the model file is valid SBML
            results = cobra.io.validate_sbml_model(temp_file.name)
            errors = results[1]['SBML_ERROR']
            self.assertTrue(0, len(errors), msg=f"SBML validation errors: {errors}")

    def test_create_cobra_reaction(self):
        """Test create_cobra_reaction function, by adding a reaction
        (and the corresponding metabolites) to the E.coli core model."""

        # Load the E. coli core model
        model = cobra.io.load_model("textbook")

        # Load the (subset of the) modelseed databases
        rxn_db = json.load(
            open(os.path.join(TESTFILE_DIR, "mini_modelseed_db_reactions.json"))
        )
        met_db = json.load(
            open(os.path.join(TESTFILE_DIR, "mini_modelseed_db_compounds.json"))
        )

        # Convert the ModelSEED databases to dictionaries for easy searching
        template_rxn_db = {rxn["id"]: rxn for rxn in rxn_db if not rxn["is_obsolete"]}
        template_met_db = {met["id"]: met for met in met_db if not met["is_obsolete"]}

        # Create a new reaction
        new_rxn, new_mets = curation.create_cobra_reaction(
            model, template_rxn_db, template_met_db, "EX_glc__D_e"
        )
        # Add the new reaction to the model
        model.add_reactions([new_rxn])
        model.add_metabolites(new_mets)
        # Check that the reaction was added
        self.assertIn("EX_glc__D_e", [rxn.id for rxn in model.reactions])
