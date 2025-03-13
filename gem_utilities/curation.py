import cobra


def process_intervention(model, df, template_rxn_db):
    # --- Process the interventions __-
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


def create_cobra_reaction(model, modelseed_db, rxn_id):
    """
    Create a COBRApy Reaction object from a ModelSEED database entry.
    Parameters:
        model (cobra.Model): Model to which the reaction and new metabolites will be added.
        modelseed_db (dict): Dictionary of ModelSEED reaction entries.
        rxn_id (str): Reaction ID to add.
    Returns:
        cobra.Reaction: The created reaction.
    """
    rxn = modelseed_db[rxn_id]
    reaction_id = rxn["id"] + "_c0"  # Assume cytosolic compartment
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
        met_name = parts[4].strip('"') if len(parts) >= 5 else met_id
        # Add metabolite if not present
        if met_id not in model.metabolites:
            met_obj = cobra.Metabolite(
                id=met_id, name=met_name, compartment=compartment
            )
            mets_to_add.append(met_obj)
        else:
            met_obj = model.metabolites.get_by_id(met_id)
        reaction.add_metabolites({met_obj: coeff})
    model.add_reactions([reaction])
    return reaction, mets_to_add
