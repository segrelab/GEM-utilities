import warnings

import cobra


def clean_media(model: cobra.Model, media: dict) -> dict:
    """clean_media
    Removes exchange reactions from the media that are not present in the model

    Parameters
    ----------
    model : cobra.Model
        The model to set the media for.
    media : dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction.

    Returns
    -------
    dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction
    """
    # Make an empty dictionary for the media
    clean_medium = {}
    # Loop through the media and set the exchange reactions that are present
    for ex_rxn, lb in media.items():
        if ex_rxn in [r.id for r in model.reactions]:
            clean_medium[ex_rxn] = lb
        else:
            warnings.warn(
                "Model does not have the exchange reaction "
                + ex_rxn
                + ", so it was not set in the media."
            )

    # Return the clean medium
    return clean_medium
