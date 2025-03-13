import cobra
import pandas as pd


def process_intervention(
    model: cobra.Model, df: pd.DataFrame, template_rxn_db: dict
) -> cobra.Model:
    """
    Process interventions on a metabolic model based on a DataFrame of changes.

    Parameters:
        model (cobra.Model): The original metabolic model.
        df (pandas.DataFrame): DataFrame containing intervention details with columns 'reaction', 'change', and 'gene'.
        template_rxn_db (dict): Dictionary of template reactions for adding new reactions.

    Returns:
        cobra.Model: The modified metabolic model after applying interventions.
    """
    model_new = model.copy()
    report = []  # to store a report of interventions
    for idx, row in df.iterrows():
        rxn_id = row["reaction"].strip()  # remove any extra whitespace
        change = row["change"].strip() if " change" in row else row["change"].strip()
        gene = row["gene"].strip() if " gene" in row else row["gene"].strip()
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
                new_rxn, new_mets = create_cobra_reaction(
                    model_new, template_rxn_db, rxn_id
                )
                # Only add gene information if the gene is not "unknown"
                if gene.lower() != "unknown":
                    # Here you might update the gene_reaction_rule or add gene objects,
                    # e.g., new_rxn.gene_reaction_rule = gene
                    new_rxn.gene_reaction_rule = gene
                    report.append(f"Added reaction {rxn_id} with gene {gene}")
                else:
                    report.append(
                        f"Added reaction {rxn_id} without gene info (gene is unknown)"
                    )
                model_new.add_reactions([new_rxn])
                model_new.add_metabolites(new_mets)
        else:
            report.append(f"Unrecognized change '{change}' for reaction {rxn_id}")
    # --- Print the report ---
    print("Intervention Report:")
    for line in report:
        print(line)
    return model_new


def create_cobra_reaction(
    model: cobra.Model,
    modelseed_rxn_db: dict,
    modelseed_cpd_db: dict,
    rxn_id: str,
    compartment: str,
) -> tuple:
    """
    Create a COBRApy Reaction object from a ModelSEED database entry.

    Parameters:
        model (cobra.Model): Model to which the reaction and new metabolites will be added.
        modelseed_rxn_db (dict): Dictionary of ModelSEED reaction entries.
        modelseed_cpd_db (dict): Dictionary of ModelSEED compound entries.
        rxn_id (str): Reaction ID to add.
        compartment (str): Compartment code for the reaction.

    Returns:
        cobra.Reaction: The created reaction.
        list: List of new metabolites added to the model.
    """
    rxn = modelseed_rxn_db[rxn_id]
    # TODO: Check if the reaction is a transport reaction and handle it accordingly
    reaction_id = (
        rxn["id"] + "_c0"
    )  # Assume cytosolic compartment  # FIXME: Hardcoded compartment
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
            compartment = "c0"
        elif comp_code == "1":
            compartment = "e0"
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
    model.add_reactions([reaction])
    return reaction, mets_to_add


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
