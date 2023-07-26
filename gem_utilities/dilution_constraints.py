'''
dilution_constraints.py
Written by Devlin Moyer; last updated on 2023-07-26

Adds a "dilution reaction" for each metabolite in the model (or a given list of
metabolites in the given model) and constrains the flux through that reaction to
be 1/dil_factor * the sum of the absolute values of the fluxes through all
reactions that involve that metabolite weighted by that metabolite's
stoichiometric coefficients in each reaction.

Requires the model to be capable of net production of all metabolites that
participate in at least one reaction with nonzero predicted flux; this prevents
"perfect" recycling of cofactors by requiring at least some net production of
them. This can improve predictions of knockout phenotypes by requiring fluxes
through biosynthetic pathways that might otherwise not be required to have flux.
Mitigates but does not completely solve the problem of unbounded loop fluxes.

This concept was introduced in https://doi.org/10.1186/gb-2010-11-4-r43 and this
implementation of the idea is based on the implementation presented in
https://doi.org/10.1371/journal.pcbi.1003126
'''

from optlang.symbolics import Zero
import cobra

def add_dilution_constraints(
    given_model, mets_to_dilute = None, split_rxns = False, dil_factor = 0,
    threads = 1, verbose = 1
):
    '''
    Given a Cobrapy Model object:
        - add a "dilution" reaction for all mets_to_dilute (if mets_to_dilute
          wasn't specified, do all metabolites in the model)
        - make all reversible reactions that involve any of the mets_to_dilute
          irreversible
        - add one new reaction for each formerly reversible reaction that goes
          in the opposite direction but is otherwise identical
        - create a constraint that sets each metabolite in mets_to_dilute's new
          dilution reaction's flux * dil_factor = sum of fluxes through all
          other reactions that metabolite participates in
    If mets_to_dilute is not specified, will use all metabolites in the GSMM
    except those that have "tRNA" in their names, because tRNA metabolites tend
    to only be capable of being recycled and lack biosynthesis reactions
    (because tRNA biosynthesis is extremely complex), and dilution constraints
    require some de novo biosynthesis of all metabolites that participate in any
    reactions that carry flux, which would make all reactions that involve tRNAs
    incapable of having flux (and they usually primarily participate in biomass
    reactions, which tend to be things people consider important)
    Return the new model object
    '''
    # modify a copy of the given model
    model = given_model.copy()
    if (dil_factor > 0) and (mets_to_dilute is None):
        # if no list of metabolites was given but a non-zero dilution factor was
        # given, find all metabolites that don't seem like tRNA metabolites,
        # since diluting those tends to cause problems
        mets_to_dilute = [
            m.id for m in model.metabolites
            # be case-insensitive
            if ('trna' not in m.name.lower()) and ('trna' not in m.id.lower())
        ]
        if verbose > 0:
            msg = f'Creating dilution constraints for {len(mets_to_dilute)} '
            msg += f'metabolites (all metabolites that don\'t have "tRNA" in '
            msg += 'their names)'
            print(msg)
    # only do this if the dilution factor isn't 0 cuz then there's no point
    if (dil_factor > 0) and (len(mets_to_dilute) > 0): 
        # raise ValueError if mets_to_dilute has Metabolite objects and not
        # IDs, cuz for unclear reasons that tends to cause problems
        if not all(isinstance(m, str) for m in mets_to_dilute):
            msg = 'mets_to_dilute can only contain the IDs of cobra.Metabolite'
            msg += 'objects and not the Metabolite objects themselves'
            raise ValueError(msg)
        # make a dilution reaction for everything in mets_to_dilute
        dil_rxns = [
            make_dilution_reaction(model.metabolites.get_by_id(m))
            for m in mets_to_dilute
        ]
        # if we're gonna predict fluxes from this model with Cobrapy, we don't
        # have to separate the reversible reactions into irreversible halves,
        # cuz that's automatically done under the hood, but if we're gonna pass
        # this model to Matlab (e.g. to run CRHMC on it), we need to separate
        # the reversible reactions in order to properly impose dilution
        # constraints in Matlab
        if split_rxns:
            # only worry about reversible reactions that involve one of the
            # mets_to_dilute
            rev_copies = [
                make_irrev(model, r) for r in model.reactions
                # don't split exchange reactions; that'd get weird
                if r.reversibility and not r.boundary and
                any(m.id in mets_to_dilute for m in r.metabolites)
            ]
            # make_irrev altered the IDs of some existing reactions; that will
            # cause problems if we don't call model.repair() right now
            model.repair()
            # add all the new reactions to the model
            model.add_reactions(dil_rxns + rev_copies)
            if verbose > 0:
                msg = f'Adding dilution reactions for {len(mets_to_dilute)} '
                msg += f'metabolites increased the total number of reactions in '
                msg += f'the model from {len(given_model.reactions)} to '
                msg += f'{len(model.reactions)}'
                print(msg)
        else:
            # if we split the reversible reactions, Cobrapy/optlang will get
            # weird about the dilution constraints, and we probably only did it
            # cuz we're about to write this model as a Matlab file, and Cobrapy
            # doesn't keep any extra constraints you add when writing models to
            # files regardless of the format, so only make these constraints if
            # we haven't split the reversible reactions
            model.add_reactions(dil_rxns)
            dil_consts = [
                make_dilution_constraint(model, met_id, dil_factor)
                for met_id in mets_to_dilute
            ]
            # Cobrapy works a lot faster if you add a list of constraints all at
            # once than if you add each constraint individually
            model.add_cons_vars(dil_consts)
    else:
        if verbose > 0:
            msg = 'Since dil_factor was 0 and/or mets_to_dilute was '
            msg += 'an empty list, no dilution constraints were added'
            print(msg)
    return(model)

def make_dilution_reaction(met):
    '''
    Make a reaction that irreversibly consumes the given metabolite and make
    the reaction ID <metabolite ID>_dilution so they're easy to find later
    '''
    dil_rxn = cobra.Reaction(f'{met.id}_dilution')
    dil_rxn.name = f'{met.name} Dilution'
    dil_rxn.lower_bound = 0
    dil_rxn.upper_bound = float('Inf')
    dil_rxn.add_metabolites({met : -1.0})
    return(dil_rxn)

def make_irrev(model, given_rxn):
    '''
    Given a reversible reaction and the model that contains it, make two new
    reaction objects with identical stoichiometry and GPRs, but set their
    bounds so that one can only go forwards and one can only go backwards
    '''
    # constrain the given reaction to only go forwards
    given_rxn.lower_bound = 0
    # make a new reaction with the same metabolites but on opposite sides of
    # the equation that can also only go forwards (also the same GPR and bound)
    rev_copy = cobra.Reaction(given_rxn.id + '__rev')
    rev_copy.subtract_metabolites(given_rxn.metabolites)
    rev_copy.gene_reaction_rule = given_rxn.gene_reaction_rule
    rev_copy.upper_bound = given_rxn.upper_bound
    # add a __fwd to the given reaction ID so that you know it started out
    # reversible later
    given_rxn.id = given_rxn.id + '__fwd'
    return(rev_copy)

def make_dilution_constraint(model, met_id, dil_factor):
    # start out with an optlang/SymPy expression that's just zero, and add in
    # the optlang/SymPy variables associated with the reactions that this
    # metabolite participates in
    expression = Zero
    # passing around Metabolite objects directly has frequently given us weird
    # errors, so we're passing around the metabolite IDs instead
    met_obj = model.metabolites.get_by_id(met_id)
    for r in met_obj.reactions:
        if 'dilution' not in r.id:
            # add both the forward and reverse variable so this is always
            # positive, regardless of which direction the reaction is going in
            expression += r.forward_variable + r.reverse_variable
        else:
            # if this is the metabolite's dilution reaction, subtract its
            # flux times the dilution factor (since it's irreversible, reverse
            # variable should always be 0, but subtracted just in case)
            expression -= dil_factor * (r.forward_variable - r.reverse_variable)
    # set the upper and lower bounds on this constraint to 0 so that the
    # dilution flux (scaled by the dilution factor) must equal the sum of fluxes
    # through all other reactions involving this metabolite
    dilution_constraint = model.problem.Constraint(
        expression, lb = 0, ub = 0, name = f'{met_id}_dilution_constraint'
    )
    return(dilution_constraint)
