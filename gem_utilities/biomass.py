from typing import List

import cobra


def unlump_biomass(
    model: cobra.Model,
    biomass_compounds: List[str],
    lumped_metabolites=["cpd11461_c0", "cpd11463_c0", "cpd11462_c0"],
):
    """
    "Un-lump" the biomass so that any lumped biomass component (e.g. DNA) is
    separated into its constituent metabolites (e.g. dAMP, dCMP, dGMP, dTMP)

    Args:
    model (cobra.Model): The model to test
    biomass_compounds (list): List of biomass compounds to test
    lumped_metabolites (list): List of metabolites that are lumped in the
        biomass reaction (defaults to
        ['cpd11461_c0', 'cpd11463_c0', 'cpd11462_c0'] which is DNA, protein,
        and RNA respectively)

    Returns:
    list: List of individual metabolites that make up the biomass
    """
    # Make a copy of the biomass compounds
    unlumped_compounds = biomass_compounds.copy()

    # Loop through the lumped metabolites
    for lumped_metabolite in lumped_metabolites:
        # Remove the metabolite from the list of biomass compounds
        unlumped_compounds.remove(lumped_metabolite)
        # Find the reaction in the model that makes the lumped metabolite
        synth_rxn = [
            r
            for r in model.reactions
            if model.metabolites.get_by_id(lumped_metabolite) in r.products
        ]
        # Throw an error if there is not exactly one reaction that makes the lumped metabolite
        if len(synth_rxn) != 1:
            raise ValueError(
                "There should be exactly one reaction that makes the lumped metabolite"
            )
        # Get the reaction
        synth_rxn = synth_rxn[0]
        # Get the metabolites that are consumed in the reaction
        substrates = [m.id for m in synth_rxn.reactants]
        # Add the metabolites to the list of biomass compounds
        unlumped_compounds += substrates

    # Return the list of individual metabolites that make up the biomass
    return unlumped_compounds
