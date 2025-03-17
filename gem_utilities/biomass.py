import os
from typing import List

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from gem_utilities.media import clean_media


def unlump_biomass(
    model: cobra.Model,
    biomass_compounds: List[str],
    lumped_metabolites=["cpd11461_c0", "cpd11463_c0", "cpd11462_c0"],
) -> List[str]:
    """
    "Un-lump" the biomass so that any lumped biomass component (e.g. DNA) is
    separated into its constituent metabolites (e.g. dAMP, dCMP, dGMP, dTMP)

    Args:
    model (cobra.Model): The model to test
    biomass_compounds (list): List of biomass compounds to test
    lumped_metabolites (list): List of metabolites that are lumped in the
        biomass reaction (defaults to
        ['cpd11461_c0', 'cpd11463_c0', 'cpd11462_c0'] which is DNA, protein,
        and RNA respectively)

    Returns:
    list: List of individual metabolites that make up the biomass
    """
    # Make a copy of the biomass compounds
    unlumped_compounds = biomass_compounds.copy()

    # Loop through the lumped metabolites
    for lumped_metabolite in lumped_metabolites:
        # Remove the metabolite from the list of biomass compounds
        unlumped_compounds.remove(lumped_metabolite)
        # Find the reaction in the model that makes the lumped metabolite
        synth_rxn = [
            r
            for r in model.reactions
            if model.metabolites.get_by_id(lumped_metabolite) in r.products
        ]
        # Throw an error if there is not exactly one reaction that makes the lumped metabolite
        if len(synth_rxn) != 1:
            raise ValueError(
                "There should be exactly one reaction that makes the lumped metabolite"
            )
        # Get the reaction
        synth_rxn = synth_rxn[0]
        # Get the metabolites that are consumed in the reaction
        substrates = [m.id for m in synth_rxn.reactants]
        # Add the metabolites to the list of biomass compounds
        unlumped_compounds += substrates

    # Return the list of individual metabolites that make up the biomass
    return unlumped_compounds


def test_biomass_producibility(
    model, growth_phenotypes, media_definitions, biomass_rxn="bio1_biomass", out_dir="."
) -> None:
    """
    Add exchange reactions for all metabolites in the model and test the producibility of the biomass components on the given media

    Args:
    model (cobra.Model): The model to test
    growth_phenotypes (pandas.DataFrame): A dataframe with columns "minimal_media", "c_source", "met_id", and "growth"
    media_definitions (dict): A dictionary of media definitions
    biomass_rxn (str): The ID of the biomass reaction in the model. Default is "bio1_biomass" (which is used in KBase models).
    out_dir (str): The directory to save the results to. Default is the current directory.

    Returns:
    None. The results are saved to a CSV file and a heatmap plot is saved to a file
    """
    # Check that the growth phenotypes dataframe has the expected columns
    expected_columns = ["minimal_media", "c_source", "met_id", "growth"]
    if not all(col in growth_phenotypes.columns for col in expected_columns):
        raise ValueError(
            "growth_phenotypes dataframe must have columns "
            + ", ".join(expected_columns)
        )

    # Add sink reactions for all metabolites, but set the lower bound to 0
    # because by default the sink reactions are reversible, and so can be
    # used to import metabolites that are not in the media
    for metabolite in model.metabolites:
        # Check if there is already a sink reaction for this metabolite
        if "SK_" + metabolite.id not in [r.id for r in model.reactions]:
            model.add_boundary(metabolite, type="sink", lb=0)

    # Get the biomass composition from the model
    biomass_rxn = model.reactions.get_by_id(biomass_rxn)
    biomass_compounds = [
        met.id for met in biomass_rxn.metabolites if biomass_rxn.metabolites[met] < 0
    ]

    # "Un-lump" the biomass so that any lumped biomass component (e.g. DNA) is
    # separated into its constituent metabolites (e.g. dAMP, dCMP, dGMP, dTMP)
    unlumped_compounds = unlump_biomass(model, biomass_compounds)

    # Make a dictionary to store the producibility results
    biomass_producibility = {}

    # Check the producibility of the biomass componets on the different carbon sources
    for index, row in growth_phenotypes.iterrows():
        # Make an ID for the results that is combination of the minimal media name and the carbon source
        c_source = row["minimal_media"] + "_" + row["c_source"]
        # Make a dictionary to store the results for just this carbon source
        biomass_producibility[c_source] = {}
        # Set the model media to match the experimental media
        medium = media_definitions[row["minimal_media"]].copy()
        medium["EX_" + row["met_id"] + "_e0"] = (
            1000.0  # FIXME: I should set this to a consistent, lower value
        )
        # Test it
        biomass_producibility[c_source] = try_biomass_in_one_medium(
            medium, unlumped_compounds, model
        )

    # Make a dataframe of the producibility results and save it to a CSV file
    df = pd.DataFrame.from_dict(biomass_producibility)
    # Save the dataframe to a CSV file and make the file name specific the the model.id
    df.to_csv(os.path.join(out_dir, model.id + "_biomass_producibility.csv"))

    # Plot the producibility results
    plot_biomass_prodcubility(model, df, out_dir=out_dir)


def try_biomass_in_one_medium(medium_dict, biomass_compounds, model):
    """Runs FBA for each biomass_compound as objective under `medium_dict`.
    Returns a dict {compound_id: True/False}."""
    results = {}
    # Set the medium
    model.medium = clean_media(model, medium_dict)

    # Loop over each biomass compound
    for cpd_id in biomass_compounds:
        # Get the human-friendly name of the metabolite
        cpd_name = model.metabolites.get_by_id(cpd_id).name

        # Set the objective to the sink/demand for cpd_id
        # e.g., "SK_cpd_id"
        model.objective = {model.reactions.get_by_id("SK_" + cpd_id): 1}

        # Optimize
        sol = model.optimize()

        # Check feasibility
        if (sol.status == "optimal") and (sol.objective_value > 1e-6):
            results[cpd_name] = True
        else:
            results[cpd_name] = False
    return results


def plot_biomass_prodcubility(model: cobra.Model, df: pd.DataFrame, out_dir="."):
    """
    Plot the producibility results as a heatmap

    Args:
    model (cobra.Model): The model that was tested
    df (pandas.DataFrame): The producibility results as a dataframe

    Returns:
    None. The plot is saved to a file.
    """
    # Convert the producibility results to a boolean
    df_binary = df.replace({True: 1, False: 0})

    # Define my own color palette
    my_palette = sns.color_palette(["red", "green"])

    # Plot the producibility results as a heatmap
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(
        df_binary,
        annot=False,
        cmap=my_palette,
        cbar_kws={"ticks": [0.25, 0.75], "label": "Producibile"},
    )

    # Now fix the colorbar tick labels
    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(["No", "Yes"])

    # Add white lines to separate the different cells
    for i in range(df_binary.shape[0]):
        plt.axhline(i, color="white", linewidth=0.5)
    for i in range(df_binary.shape[1]):
        plt.axvline(i, color="white", linewidth=1)

    # Make sure all y ticks/component names are shown
    plt.yticks(
        [x + 0.5 for x in range(df_binary.shape[0])], df_binary.index, fontsize=5
    )

    plt.xlabel("Biomass composition")
    plt.ylabel("Biomass compound")

    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(out_dir, model.id + "_biomass_producibility_heatmap.png"))
