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
