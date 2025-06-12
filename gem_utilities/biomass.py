import os
from typing import List

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from gem_utilities.media import clean_media


def unlump_biomass(
    biomass_metabolites: dict,
    model: cobra.Model,
    lumped_metabolites=["cpd11461_c0", "cpd11463_c0", "cpd11462_c0"],
) -> List[str]:
    """
    "Un-lump" the biomass so that any lumped biomass component (e.g. DNA) is
    separated into its constituent metabolites (e.g. dAMP, dCMP, dGMP, dTMP)

    Args:
    biomass_metabolites (dict): Metabolites attribute of the model's biomass
        reaction, where the keys are metabolite IDs and the values are their
        stoichiometric coefficients
    model (cobra.Model): The model to use for the biomass reaction
    lumped_metabolites (list): List of metabolites that are lumped in the
        biomass reaction (defaults to
        ['cpd11461_c0', 'cpd11463_c0', 'cpd11462_c0'] which is DNA, protein,
        and RNA respectively)

    Returns:
    list: List of individual metabolites that make up the biomass
    """
    # Make a new metabolite/coefficient dictionary for the biomass reaction
    # with "unlumped" components
    unlumped_metabolites = {}

    for metabolite, coeff in biomass_metabolites.items():
        # If the metabolite is a lumped component, skip it
        if metabolite.id in lumped_metabolites:
            # Find the reactions that produces the lumped metabolite
            # Find the reaction in the model that makes the lumped metabolite
            synth_rxn = [r for r in model.reactions if metabolite in r.products]
            # Throw an error if there is not exactly one reaction that makes the lumped metabolite
            if len(synth_rxn) != 1:
                raise ValueError(
                    "There should be exactly one reaction that makes the lumped metabolite"
                )
            # Get the reaction
            synth_rxn = synth_rxn[0]
            # Loop trhough all of the metabolites in the reaction
            for subcomponent, sub_coeff in synth_rxn.metabolites.items():
                # If the metabolite is not a lumped component, add it to the list of biomass compounds
                if subcomponent.id not in lumped_metabolites:
                    # Add the metabolite and coefficient to the dictionary
                    update_dict(
                        unlumped_metabolites,
                        {subcomponent: sub_coeff * abs(coeff)},
                    )
        # Otherwise, add the metabolite and coefficient to the dictionary
        else:
            update_dict(unlumped_metabolites, {metabolite: coeff})

    # Return the new dictionary of unlumped metabolites
    return unlumped_metabolites


def update_dict(d, new_items):
    """Update a dictionary with new items."""
    for key, value in new_items.items():
        if key in d:
            d[key] += value
        else:
            d[key] = value


def check_biomass_producibility(
    model: cobra.Model,
    growth_phenotypes: pd.DataFrame,
    media_definitions: dict,
    media_negative_controls: bool = True,
    sinks_for_all: bool = True,
    biomass_rxn: str = "bio1_biomass",
    external_compartment: str = "e0",
    lumped_biomass_components: List[str] = [
        "cpd11461_c0",
        "cpd11463_c0",
        "cpd11462_c0",
    ],
    out_dir: str = ".",
) -> None:
    """
    _summary_

    Parameters
    ----------
    model : cobra.Model
        The model to test
    growth_phenotypes : pd.DataFrame
        A dataframe with columns "minimal_media", "c_source", "met_id", and
        "growth"
    media_definitions : dict
        A dictionary of media definitions
    media_negative_controls : bool, optional
        Do you want to show negative controls for each media on the heatmap?,
        by default True
    sinks_for_all : bool, optional
        Do you want to add sink reactions for all metabolites in the model or
        just for the other biomass compontents?, by default True
    biomass_rxn : str, optional
        the reaction ID for the biomass reaction in the model, by default
        "bio1_biomass"
    external_compartment : str, optional
        the comparement ID for the extracellular compartment in the model, by
        default "e0"
    lumped_biomass_components : List[str], optional
        List of the biomass components which are pseudo-metabolites to break
        into their constituent parts (e.g. DNA, RNA, and protein), by default
        [ "cpd11461_c0", "cpd11463_c0", "cpd11462_c0"]
    out_dir : str, optional
        the directory in which to save the results, by default "."

    Raises
    ------
    ValueError
        The growth_phenotypes dataframe does not have the required columns:
        minimal_media, c_source, met_id, and growth
    """
    # Check that the growth phenotypes dataframe has the expected columns
    expected_columns = ["minimal_media", "c_source", "met_id", "growth"]
    if not all(col in growth_phenotypes.columns for col in expected_columns):
        raise ValueError(
            "growth_phenotypes dataframe must have columns "
            + ", ".join(expected_columns)
        )

    # Get the biomass composition from the model
    biomass_rxn = model.reactions.get_by_id(biomass_rxn)

    # "Un-lump" the biomass so that any lumped biomass component (e.g. DNA) is
    # separated into its constituent metabolites (e.g. dAMP, dCMP, dGMP, dTMP)
    if lumped_biomass_components:
        unlumped_compounds = [
            met.id
            for met in unlump_biomass(model, biomass_rxn.metabolites)
            if biomass_rxn.metabolites[met] < 0
        ]
    else:
        unlumped_compounds = [
            met.id
            for met in biomass_rxn.metabolites
            if biomass_rxn.metabolites[met] < 0
        ]

    # Add sinks, either for all metabolites, or juts for the biomass components
    if sinks_for_all:
        # Add sink reactions for all metabolites, but set the lower bound to 0
        # because by default the sink reactions are reversible, and so can be
        # used to import metabolites that are not in the media
        for metabolite in model.metabolites:
            # Check if there is already a sink reaction for this metabolite
            if "SK_" + metabolite.id not in [r.id for r in model.reactions]:
                model.add_boundary(metabolite, type="sink", lb=0)
        # Set the string about sinks to be added to the results files
        sinks = "all_sinks"
    else:
        # Add sink reactions for just the biomass components
        for met_id in unlumped_compounds:
            metabolite = model.metabolites.get_by_id(met_id)
            # Check if there is already a sink reaction for this metabolite
            if "SK_" + metabolite.id not in [r.id for r in model.reactions]:
                model.add_boundary(metabolite, type="sink", lb=0)
        # Set the string about sinks to be added to the results files
        sinks = "biomass_sinks_only"

    # Make a dictionary to store the producibility results
    biomass_producibility = {}

    # Add negative controls if requested
    if media_negative_controls:
        controls = {"empty": {}}
        # Loop through the media definitions and add a negative control for each
        unique_media = growth_phenotypes["minimal_media"].unique()
        for media in unique_media:
            controls[media] = media_definitions[media].copy()
        # Test the negative controls
        for cntrl_name, cntrl_medium in controls.items():
            biomass_producibility[cntrl_name] = try_biomass_in_one_medium(
                cntrl_medium, unlumped_compounds, model
            )

    # Check the producibility of the biomass componets on the different carbon sources
    for index, row in growth_phenotypes.iterrows():
        # Make an ID for the results that is combination of the minimal media name and the carbon source
        c_source = row["minimal_media"] + "_" + row["c_source"]
        # Make a dictionary to store the results for just this carbon source
        biomass_producibility[c_source] = {}
        # Set the model media to match the experimental media
        medium = media_definitions[row["minimal_media"]].copy()
        if not pd.isna(row["met_id"]):
            medium["EX_" + row["met_id"] + "_" + external_compartment] = (
                1000.0  # FIXME: I should set this to a consistent, lower value
            )
        # Test it
        biomass_producibility[c_source] = try_biomass_in_one_medium(
            medium, unlumped_compounds, model
        )

    # Make a dataframe of the producibility results and save it to a CSV file
    df = pd.DataFrame.from_dict(biomass_producibility)
    # Save the dataframe to a CSV file and make the file name specific the the model.id
    df.to_csv(
        os.path.join(out_dir, model.id + "_biomass_producibility_" + sinks + ".csv")
    )

    # Plot the producibility results
    plot_biomass_prodcubility(model, df, sinks, out_dir=out_dir)


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


def plot_biomass_prodcubility(model: cobra.Model, df: pd.DataFrame, sinks, out_dir="."):
    """
    Plot the producibility results as a heatmap

    Args:
    model (cobra.Model): The model that was tested
    df (pandas.DataFrame): The producibility results as a dataframe
    sinks (str): The name of the sink reactions that were added

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
    plt.savefig(
        os.path.join(
            out_dir, model.id + "_biomass_producibility_heatmap_" + sinks + ".png"
        )
    )


def calculate_biomass_weight(
    model: cobra.Model,
    biomass_rxn: str = "bio1_biomass",
    lumped_biomass_components: List[str] = [
        "cpd11461_c0",
        "cpd11463_c0",
        "cpd11462_c0",
    ],
    save_work_table: bool = False,
    out_dir: str = None,
) -> float:
    """
    Calculate the weight of the biomass reaction in the model.

    Parameters
    ----------
    model : cobra.Model
        The model to use for the biomass reaction.
    biomass_rxn : str, optional
        The ID of the biomass reaction, by default "bio1_biomass"
    lumped_biomass_components : List[str], optional
        List of the biomass components which are pseudo-metabolites to break
        into their constituent parts (e.g. DNA, RNA, and protein), by default
        [ "cpd11461_c0", "cpd11463_c0", "cpd11462_c0"]
    save_work_table : bool, optional
        If True, save the work table to a CSV file, by default False
    out_dir : str, optional
        The directory in which to save the work table, by default None

    Returns
    -------
    float
        Weight of the biomass reaction in grams per mole (which is the unit of
        molecular mass).
    """
    # If save_work_table is True, make sure that out_dir is set
    if save_work_table and out_dir is None:
        raise ValueError(
            "If save_work_table is True, out_dir must be set to a valid directory."
        )

    # Get the biomass reaction
    biomass_rxn = model.reactions.get_by_id(biomass_rxn)

    # "Un-lump" the biomass so that any lumped biomass component (e.g. DNA) is
    # separated into its constituent metabolites (e.g. dAMP, dCMP, dGMP, dTMP)
    if lumped_biomass_components:
        # Check that the lumped biomass components are in the model
        for lumped_met in lumped_biomass_components:
            if lumped_met not in [m.id for m in model.metabolites]:
                raise ValueError(
                    f"Lumped biomass component {lumped_met} is not in the model."
                )
        # Unlump the biomass reaction metabolites to get the new stoichiometry
        unlumped_stoichiometry = unlump_biomass(
            biomass_rxn.metabolites,
            model,
            lumped_metabolites=lumped_biomass_components,
        )
    else:
        unlumped_stoichiometry = biomass_rxn.metabolites

    # Make sure that the stoichiometry is a dictionary
    if not isinstance(unlumped_stoichiometry, dict):
        raise ValueError(
            "The stoichiometry of the biomass reaction is not a dictionary."
        )

    # Make sure that all of the metabolites in the stoichiometry have a formula weight that is not 0
    for metabolite in unlumped_stoichiometry:
        if not hasattr(metabolite, "formula_weight") or metabolite.formula_weight == 0:
            raise ValueError(
                f"The metabolite {metabolite.id} does not have a formula weight."
            )

    # Calculate the weight of the biomass reaction
    weight = 0.0
    if save_work_table:
        # Create a list to store the work table
        work_table = []
    # Loop through the metabolites in the biomass reaction
    for metabolite, coeff in unlumped_stoichiometry.items():
        # Multiply the formula weight of the metabolite by its coefficient
        # Use the opposite sign of the coefficient because the biomass weight
        # should include the consumed metabolites (negative coefficient) and
        # not the produced ones (positive coefficient)
        weight += metabolite.formula_weight * (-1 * coeff)
        # Save the information to the work table if requested
        if save_work_table:
            work_table.append(
                {
                    "metabolite": metabolite.id,
                    "coefficient": coeff,
                    "formula": metabolite.formula,
                    "formula_weight": metabolite.formula_weight,
                    "weight_contribution": metabolite.formula_weight * (-1 * coeff),
                }
            )

    # If requested, save the work table to a CSV file
    if save_work_table:
        # Convert the work table to a DataFrame
        work_table_df = pd.DataFrame(work_table)
        # Add a row for the total weight of the biomass reaction
        work_table_df = work_table_df.append(
            {
                "metabolite": "Total",
                "coefficient": "",
                "formula": "",
                "formula_weight": "",
                "weight_contribution": weight,
            },
            ignore_index=True,
        )
        # Save the DataFrame to a CSV file
        work_table_df.to_csv(
            os.path.join(out_dir, model.id + "_biomass_weight_work_table.csv"),
            index=False,
        )
    # Return the weight of the biomass reaction
    if weight < 0:
        raise ValueError("The biomass weight cannot be negative. Check the model.")
    if weight == 0:
        raise ValueError("The biomass weight cannot be zero. Check the model.")
    return weight
