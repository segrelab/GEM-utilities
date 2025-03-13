import warnings

import cobra
import pandas as pd


def process_intervention_df(
    model: cobra.Model,
    df: pd.DataFrame,
    modelseed_rxn_db: dict,
    modelseed_cpd_db: dict,
    rxn_compartment: str = "c0",
    internal_compartment: str = "c0",
    external_compartment: str = "e0",
) -> cobra.Model:
    """
    Process interventions on a metabolic model based on a DataFrame of changes.

    Parameters:
        model (cobra.Model): The original metabolic model.
        df (pandas.DataFrame): DataFrame containing intervention details with
            columns 'reaction', 'change', and 'gene'.
        modelseed_rxn_db (dict): Dictionary of template reactions for adding
            new reactions.
        model_seed_cpd_db (dict): Dictionary of template compounds for adding
            new metabolites.
        rxn_compartment (str): Compartment code for the reactions. Default is
            "c0". Note- there is currently no way to specify different
            compartments for different reactions.
        internal_compartment (str): Compartment code for internal metabolites.
            Default is "c0".
        external_compartment (str): Compartment code for external metabolites.
            Default is "e0".

    Returns:
        cobra.Model: The modified metabolic model after applying interventions.
    """
    model_new = model.copy()
    report = []  # to store a report of interventions
    for idx, row in df.iterrows():
        rxn_id = row["reaction"].strip()  # remove any extra whitespace
        change = row[" change"].strip() if " change" in row else row["change"].strip()
        gene = row[" gene"].strip() if " gene" in row else row["gene"].strip()
        # Give a warning if the gene is not empty or "unknown"
        if gene and gene.lower() != "unknown":
            warnings.warn(
                f"A gene name, {gene}, was given for reaction {rxn_id}, but will not be added to the model."
            )
        if change == "-":
            # Removal intervention: check if the reaction exists
            if rxn_id in model_new.reactions:
                rxn_obj = model_new.reactions.get_by_id(rxn_id)
                model_new.remove_reactions([rxn_obj])
                report.append(f"Removed reaction {rxn_id}")
            else:
                report.append(f"Reaction {rxn_id} not found; nothing removed")
        elif change == "+":
            # Addition intervention: check if the reaction is already present
            if rxn_id in model_new.reactions:
                report.append(f"Reaction {rxn_id} already exists; not added")
            else:
                report.append(
                    f"Adding reaction {rxn_id} to compartment {rxn_compartment}"
                )
                # TODO: Add a column to the dataframe for reaciton compartment
                # and pass that to this funciton as the rxn_compartment
                new_model = add_ms_reaction_from_id(
                    model_new,
                    modelseed_rxn_db,
                    modelseed_cpd_db,
                    rxn_id,
                    rxn_compartment=rxn_compartment,
                    internal_compartment=internal_compartment,
                    external_compartment=external_compartment,
                )
        else:
            report.append(f"Unrecognized change '{change}' for reaction {rxn_id}")
        # Update which model to use for the next iteration
        model_new = new_model if "new_model" in locals() else model_new
    # --- Print the report ---
    print("Intervention Report:")
    for line in report:
        print(line)
    return model_new


def add_ms_reaction_from_id(
    model_original: cobra.Model,
    modelseed_rxn_db: dict,
    modelseed_cpd_db: dict,
    rxn_id: str,
    rxn_compartment: str = "c0",
    internal_compartment: str = "c0",
    external_compartment: str = "e0",
) -> cobra.Model:
    """
    Create a COBRApy Reaction object (and any new metabolites) from a ModelSEED
    database entry and add it to a model.

    Parameters:
        model_original (cobra.Model): Model to which the reaction and new metabolites will be added.
        modelseed_rxn_db (dict): Dictionary of ModelSEED reaction entries.
        modelseed_cpd_db (dict): Dictionary of ModelSEED compound entries.
        rxn_id (str): Reaction ID to add.
        rxn_compartment (str): Compartment code for the reaction.
        internal_compartment (str): Compartment code for internal metabolites. Default is "c0".
        external_compartment (str): Compartment code for external metabolites. Default is "e0".

    Returns:
        cobra.Reaction: The created reaction.
        list: List of new metabolites added to the model.
    """
    # Make a copy of the model to edit
    model = model_original.copy()
    # Get the reaction entry from the modelSEED database
    rxn = modelseed_rxn_db[rxn_id]
    # TODO: Check if the reaction is a transport reaction and handle it accordingly
    reaction_id = rxn["id"] + "_" + rxn_compartment
    reaction = cobra.Reaction(reaction_id, name=rxn["name"])
    mets_to_add = []
    # Set reaction bounds based on reversibility
    if rxn["reversibility"] == ">":
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
    elif rxn["reversibility"] == "<":
        reaction.lower_bound = -1000
        reaction.upper_bound = 0
    else:  # '=' or any unexpected value
        reaction.lower_bound = -1000
        reaction.upper_bound = 1000
    # Add metabolites to the reaction
    for met_str in rxn["stoichiometry"].split(";"):
        parts = met_str.split(":")
        if len(parts) < 3:
            continue  # Skip malformed entries
        coeff = float(parts[0])
        met_id_base = parts[1]
        comp_code = parts[2]
        # Map compartment code to tag
        if comp_code == "0":
            compartment = internal_compartment
        elif comp_code == "1":
            compartment = external_compartment
        else:
            compartment = comp_code
        met_id = f"{met_id_base}_{compartment}"
        # Add metabolite if not present
        # TODO: Make the metabolite more detailed (e.g. add formula, charge, etc.)
        if met_id not in model.metabolites:
            met_obj = create_cobra_metabolite(
                modelseed_cpd_db, met_id_base, compartment
            )
            mets_to_add.append(met_obj)
        else:
            met_obj = model.metabolites.get_by_id(met_id)
        reaction.add_metabolites({met_obj: coeff})
    # Add the reaction to the model
    model.add_reactions([reaction])
    # Add the new metabolites to the model
    model.add_metabolites(mets_to_add)
    # Handle tranport reactions: Add corresponding exchange reactions
    if rxn["is_transport"]:
        # Get all of the external metabolites that have been added
        external_mets = [
            met for met in mets_to_add if met.compartment == external_compartment
        ]
        # Add exchange reactions for each external metabolite
        # Set the lower bound to 0 to allow only for release
        for met in external_mets:
            model.add_boundary(met, type="exchange", lb=0, ub=1000)
    # Return the model
    return model


def create_cobra_metabolite(
    modelseed_cpd_db: dict, cpd_id: str, compartment: str
) -> cobra.Metabolite:
    """
    Create a COBRApy Metabolite object from a ModelSEED database entry.

    Parameters:
        modelseed_cpd_db (dict): Dictionary of ModelSEED compound entries.
        cpd_id (str): Compound ID to add.

    Returns:
        cobra.Metabolite: The created metabolite.
    """
    # Get the compound entry from the modelSEED database
    cpd = modelseed_cpd_db[cpd_id]

    # Make the metabolite ID compartment-specific and retrieve the name
    met_id = f"{cpd_id}_{compartment}"
    met_name = cpd["name"]

    # Make a new metabolite object with the ID and name
    met = cobra.Metabolite(id=met_id, name=met_name, compartment=compartment)

    # Add the rest of the information to the metabolite object
    met.formula = cpd["formula"] if "formula" in cpd else ""
    met.charge = cpd["charge"] if "charge" in cpd else ""
    met.annotation = parse_aliases(cpd["aliases"]) if "aliases" in cpd else {}

    return met


def parse_aliases(aliases: list) -> dict:
    """
    Parse the aliases list to extract relevant annotations.

    Parameters:
        aliases (list): List of alias strings.

    Returns:
        dict: Dictionary of parsed annotations.
    """
    annotations = {}
    for alias in aliases:
        # Split the string into a key and value around the first ":"
        parts = alias.split(":")
        if len(parts) == 2:
            key = parts[0].lower().strip()
            value = parts[1].strip()
        else:
            continue
        # Skip if the key is "name"
        if key == "name":
            continue
        # Split the value by ";" to get a list of values
        value = [v.strip() for v in value.split(";")]
        # Add the key-value pair to the annotations dictionary
        annotations[key] = value
    return annotations
