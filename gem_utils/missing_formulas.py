import logging


def copy_formulas(model, source_id, target_id):
    """Copy the chemical formula from one metabolite to another.

    Args:
        model (cobra.Model): A cobra model.
        source_id (str): The ID of the metabolite to copy the chemical formula from.
        target_id (str): The ID of the metabolite to copy the chemical formula to.

    Returns:
        None

    """
    # Get the metabolite to copy the chemical formula from
    source_met = model.metabolites.get_by_id(source_id)

    # Get the metabolite to copy the chemical formula to
    target_met = model.metabolites.get_by_id(target_id)

    # Copy the chemical formula
    target_met.formula = source_met.formula

    return None


def find_matching_metabolite(model, metabolite):
    """Find a metabolite with the same name in a different compartment.

    Args:
        model (cobra.Model): A cobra model.
        metabolite (cobra.Metabolite): A cobra metabolite.

    Returns:
        cobra.Metabolite: A cobra metabolite with the same name in a different compartment.

    """
    # Get the name of the metabolite
    name = metabolite.name

    # Get the compartment of the metabolite
    compartment = metabolite.compartment

    # Find metabolites with the same name in a different compartment
    matching_met = model.metabolites.query(lambda x: x.name == name and x.compartment != compartment)

    # Check that there is only one matching metabolite
    if len(matching_met) != 1:
        logging.warning('There is not exactly one matching metabolite.')
    else:
        # Return the matching metabolite
        return matching_met[0]
