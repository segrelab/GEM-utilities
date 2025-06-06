import json
import os
import tempfile
import unittest
from io import StringIO

import cobra
import pandas as pd

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

        # Add it to my test model
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, "test_model.xml"))
        model.add_metabolites([new_met])

        # Save the model to a temporary file
        with tempfile.NamedTemporaryFile(suffix=".xml") as temp_file:
            cobra.io.write_sbml_model(model, temp_file.name)
            # Check that the model file is valid SBML
            results = cobra.io.validate_sbml_model(temp_file.name)
            errors = results[1]["SBML_ERROR"]
            self.assertEqual(0, len(errors), msg=f"SBML validation errors: {errors}")

    def test_add_ms_reaction_from_id(self):
        """Test the add reaction function, by adding a reaction
        (and the corresponding metabolites) to my minimal model."""

        # Load the test model
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, "test_model.xml"))

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

        # Create a new reaction and add it to a new model object
        new_model = curation.add_ms_reaction_from_id(
            model, template_rxn_db, template_met_db, "rxn00001", "c0"
        )

        # Check that the reaction was added
        # I.e. that it wasn't in the original model, but is in the new model
        self.assertTrue(
            ("rxn00001_c0" in [rxn.id for rxn in new_model.reactions])
            and ("rxn00001_c0" not in [rxn.id for rxn in model.reactions])
        )

        # Save the model to a temporary file
        with tempfile.NamedTemporaryFile(suffix=".xml") as temp_file:
            cobra.io.write_sbml_model(new_model, temp_file.name)
            # Check that the model file is valid SBML
            results = cobra.io.validate_sbml_model(temp_file.name)
            errors = results[1]["SBML_ERROR"]
            self.assertEqual(0, len(errors), msg=f"SBML validation errors: {errors}")

    def test_add_transport_rxn(self):
        """Use the create_cobra_reaction function to add a transport reaction
        (rxn05145, a phosphate-transporting ATPase)
        from the ModelSEED database to my minimal model. Check that it handles
        the external metabolite and exchange reaction correctly."""
        # Load the test model
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, "test_model.xml"))

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

        # Create a new reaction and add it to the model
        # TODO: Check that transport reactions are normally called _c0 reactions
        new_model = curation.add_ms_reaction_from_id(
            model, template_rxn_db, template_met_db, "rxn05145", "c", "c", "e"
        )

        # Check that the reaction was added
        # I.e. that it wasn't in the original model, but is in the new model
        self.assertTrue(
            ("rxn05145_c" in [rxn.id for rxn in new_model.reactions])
            and ("rxn05145_c" not in [rxn.id for rxn in model.reactions])
        )

        # Check that the external metabolites were added
        self.assertIn("cpd00009_e", [met.id for met in new_model.metabolites])

        # Check that the corresponding exchange reactions were added
        self.assertTrue(
            ("EX_cpd00009_e" in [rxn.id for rxn in new_model.reactions])
            and ("EX_cpd00009_e" not in [rxn.id for rxn in model.reactions])
        )

        # Save the model to a temporary file
        with tempfile.NamedTemporaryFile(suffix=".xml") as temp_file:
            cobra.io.write_sbml_model(new_model, temp_file.name)
            # Check that the model file is valid SBML
            results = cobra.io.validate_sbml_model(temp_file.name)
            errors = results[1]["SBML_ERROR"]
            self.assertEqual(0, len(errors), msg=f"SBML validation errors: {errors}")

    def test_process_intervention_df(self):
        """Test function to process interventions from a DataFrame."""
        intervention_csv = """
        reaction, change, gene
        rxn00001,+,unknown
        rxn05145,+,unknown
        """

        # Use StringIO to simulate a file
        df = pd.read_csv(StringIO(intervention_csv.strip()))

        # Load the test model
        model = cobra.io.read_sbml_model(os.path.join(TESTFILE_DIR, "test_model.xml"))

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

        # Process the interventions
        new_model = curation.process_intervention_df(
            model, df, template_rxn_db, template_met_db, "c", "c", "e"
        )

        # Check that the reactions from the DF, and the associate exchange
        # reaction were added
        self.assertTrue(
            ("rxn00001_c" in [rxn.id for rxn in new_model.reactions])
            and ("rxn00001_c" not in [rxn.id for rxn in model.reactions])
        )
        self.assertTrue(
            ("rxn05145_c" in [rxn.id for rxn in new_model.reactions])
            and ("rxn05145_c" not in [rxn.id for rxn in model.reactions])
        )
        self.assertTrue(
            ("EX_cpd00009_e" in [rxn.id for rxn in new_model.reactions])
            and ("EX_cpd00009_e" not in [rxn.id for rxn in model.reactions])
        )

        # Save the model to a temporary file
        with tempfile.NamedTemporaryFile(suffix=".xml") as temp_file:
            cobra.io.write_sbml_model(new_model, temp_file.name)
            # Check that the model file is valid SBML
            results = cobra.io.validate_sbml_model(temp_file.name)
            errors = results[1]["SBML_ERROR"]
            self.assertEqual(0, len(errors), msg=f"SBML validation errors: {errors}")