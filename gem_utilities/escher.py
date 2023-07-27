import cobra
import warnings


def create_comparison_data_file(*models):
    """A function that takes any number of COBRA models and creates a
    JSON file to use as 'reaction data' in Escher, where the reaction
    line will be weighted based on the number of models that the reaction
    is present in. This is useful for visually comparing models.
    Ideally, there would be a function in Escher for comparing models,
    but this is a workaround for now.

    Parameters
    ----------
    models : cobra.Model
        Any number of cobra model objects

    Returns
    -------
    A dictionary that can be saved to make a reaction data JSON file
    """
    # Check that all arguments are cobra models
    for model in models:
        if not isinstance(model, cobra.Model):
            raise TypeError('All arguments must be cobra models.')

    # If there is only one model, give an error
    if len(models) == 1:
        raise ValueError('This function requires at least two models.')

    # TODO: Check that the models use the same comparement nomenclature
    # (e.g. 'c' vs 'c0') since the compartment ID is ususally in the
    # reaction IDs

    # Create a dictionary of reactions and the number of models they are
    # present in and return it
    reaction_dict = {}
    for model in models:
        for reaction in model.reactions:
            if reaction.id not in reaction_dict:
                reaction_dict[reaction.id] = 1
            else:
                reaction_dict[reaction.id] += 1
    return reaction_dict
