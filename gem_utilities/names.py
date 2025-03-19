import cobra


def fix_names_w_compartment_suffix(model):
    """Fix metabolite names with a compartment suffix.

    Args:
        model (cobra.Model): A cobra model.

    Returns:
        None

    """
    # Find metabolites with a compartment suffix
    names_w_compartment_suffix = find_names_w_compartment_suffix(model)

    # Fix the metabolite names
    for id in names_w_compartment_suffix:
        # Trim the compartment suffix
        trimmed_name = trim_name(names_w_compartment_suffix[id])

        # Find the metabolite with the compartment suffix
        met = model.metabolites.get_by_id(id)

        # Change the metabolite name
        met.name = trimmed_name

    return None


def trim_name(name):
    """Trim the compartment ID to remove the compartment suffix.

    Args:
        name (str): A cobra metabolite name.

    Returns:
        str: The cobra metabolite name without the compartment suffix.

    """
    # Get the compartment suffix
    suffix = name[-3:]

    # Check that the compartment suffix is valid
    if suffix not in ["[c]", "[e]"]:
        raise ValueError("The compartment suffix is not valid.")

    # Trim the compartment suffix
    trimmed_name = name[:-3]

    return trimmed_name


def find_names_w_compartment_suffix(model):
    """Find metabolites with a compartment suffix.

    Args:
        model (cobra.Model): A cobra model.

    Returns:
        dict: A dictionary of metabolite IDs and names.

    """
    # Find metabolites with a compartment suffix
    names_w_compartment_suffix = {
        met.id: met.name for met in model.metabolites if met.name[-3:] in ["[c]", "[e]"]
    }

    return names_w_compartment_suffix


# Write a function to go from the stoichiometry dictionary to a string of the reaction with human-friendly metabolite names
def build_reaction_string(cobra_reaction: cobra.Reaction):
    """
    Given a cobra reaction, return a string of the reaction with
    human-friendly metabolite names, and an arrow indicating the allowable
    direction of the reaction (determined from the reaction bounds).

    Args:
        cobra_reaction (cobra.Reaction): A cobra reaction.

    Returns:
        str: A string of the reaction with human-friendly metabolite names.
    """
    reactants_side = " + ".join(
        [
            f"{abs(coeff)} {met.name}"
            for met, coeff in cobra_reaction.metabolites.items()
            if coeff < 0
        ]  # TODO: Change this into another function, so I can reuse it below
    )
    products_side = " + ".join(
        [
            f"{abs(coeff)} {met.name}"
            for met, coeff in cobra_reaction.metabolites.items()
            if coeff > 0
        ]
    )
    forward_arrow = ">" if cobra_reaction.upper_bound > 0 else ""
    reverse_arrow = "<" if cobra_reaction.lower_bound < 0 else ""
    return f"{reactants_side} {reverse_arrow}--{forward_arrow} {products_side}"
