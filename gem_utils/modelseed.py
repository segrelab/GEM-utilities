from cobra import Reaction


def generate_rxn(ec_number, database):
    """Generate a reaction from an EC number using the ModelSEED database."""
    # Get the dictionary from the list of dictionaries with the correct entry in the 'ec_numbers' key
    rxn_dict = next((rxn for rxn in database if ec_number in rxn['ec_numbers']), None)

    # Create COBRA reaction
    rxn = Reaction(ec_number)
    rxn.name = rxn_dict['name']
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    # TODO: Add metabolites
    # for met_id, stoich in rxn_dict['stoichiometry'].items():
    #     rxn.add_metabolites({database['compounds'].find_one({'id': met_id})['abbreviation']: stoich})

    return rxn
