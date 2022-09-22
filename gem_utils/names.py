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
    for name in names_w_compartment_suffix:
        # Trim the compartment suffix
        trimmed_name = trim_name(name)

        # Find the metabolite with the compartment suffix
        met = model.metabolites.get_by_id(name)

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
    if suffix not in ['[c]', '[e]']:
        raise ValueError('The compartment suffix is not valid.')

    # Trim the compartment suffix
    trimmed_name = name[:-3]

    return trimmed_name


def find_names_w_compartment_suffix(model):
    """Find metabolites with a compartment suffix.

    Args:
        model (cobra.Model): A cobra model.

    Returns:
        list: A list of metabolite names with a compartment suffix.

    """
    # Find metabolites with a compartment suffix
    names_w_compartment_suffix = [met.name for met in model.metabolites if met.name[-3:] in ['[c]', '[e]']]

    return names_w_compartment_suffix
