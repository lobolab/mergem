"""
    Uses mergem dictionaries to translate models to common namespace
    and then merges the models to remove metabolite and reaction
    duplicates.
    Resulting model can be FBA simulated.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""
from . import __modelHandling
from . import __version
import cobra


# merges models in a list to the template/first model
# set_objective can be an integer for model obj or 'merge'
def merge(input_models, set_objective='merge'):
    """
    Takes a list of cobra models or file names as input and merges them into a single model with the chosen objective. \n
    :param input_models: list of cobra models or file names
    :param set_objective: objective reaction from one of the models or merge (default) all model objectives
    :return: A dictionary of the merged model, met & reac jaccard distances, num of mets and reacs merged,
            and met & reac sources.
    """
    __modelHandling.load_met_univ_id_dict()
    models = []

    for input_model in input_models:
        if isinstance(input_model, str):
            try:
                input_model = __modelHandling.load_model(input_model)
                models.append(input_model)

            except Exception as e:
                print("Error loading model: ", e)

        else:
            models.append(input_model)

    objective_reactions = []
    met_model_id_dict, met_sources_dict, reac_sources_dict, merged_model_reactions_dict = {}, {}, {}, {}

    merged_model_id = 'mergem'
    merged_model_name = 'Mergem of '
    merged_model_metabolites = []
    merged_model_reactions = []

    # Add first model
    model = models[0]
    merged_model_id += '_' + model.id
    merged_model_name += model.name if model.name else model.id
    model_objectives = []

    for metabolite in model.metabolites:
        merged_model_metabolites.append(metabolite)
        old_met_id = metabolite.id
        new_met_id = map_metabolite_to_mergem_id(metabolite)

        if (new_met_id is None) or (new_met_id in met_sources_dict):
            met_sources_dict[old_met_id] = {0}
        else:
            met_sources_dict[new_met_id] = {0}
            if new_met_id in met_model_id_dict:
                met_model_id_dict[new_met_id].append(old_met_id)
            else:
                met_model_id_dict[new_met_id] = [old_met_id]
            metabolite.id = new_met_id

    for reaction in model.reactions:
        reac_id = reaction.id
        if reac_id in str(model.objective): # processing objective reactions
            model_objectives.append(reaction)
        else:
            merged_model_reactions.append(reaction)
            reac_sources_dict[reac_id] = {0}
            reaction_key, rev_reaction_key = create_reaction_key(reaction)
            if not reaction_key in merged_model_reactions_dict:
                merged_model_reactions_dict[reaction_key] = reac_id

    objective_reactions.append(model_objectives)

    # Merge rest of models
    for model_index in range(1, len(models)):
        model = models[model_index]
        merged_model_id += '_' + model.id
        merged_model_name += '; ' + (model.name if model.name else model.id)
        model_objectives = []
        model_metabolite_ids = {m.id for m in model.metabolites}

        for metabolite in model.metabolites:
            old_met_id = metabolite.id
            new_met_id = map_metabolite_to_mergem_id(metabolite)

            if new_met_id is None:
                if old_met_id in met_sources_dict:
                    met_sources_dict[old_met_id].add(model_index)
                else:
                    met_sources_dict[old_met_id] = {model_index}
                    merged_model_metabolites.append(metabolite)
            elif old_met_id in met_sources_dict: # priority is given to original metabolite ids
                met_sources_dict[old_met_id].add(model_index)
            elif new_met_id in met_sources_dict: # new metabolite id previously found
                old_met_ids = met_model_id_dict[new_met_id]
                if (old_met_id not in old_met_ids) and any(id in old_met_ids for id in model_metabolite_ids): # model has a better match
                    met_sources_dict[old_met_id] = {model_index}
                    merged_model_metabolites.append(metabolite)
                elif model_index in met_sources_dict[new_met_id]: # model already had a metabolite for this mergem id
                    for reaction in metabolite.reactions: # replace id in its reactions
                        if new_met_id in [met.id for met in reaction.metabolites]: # new metabolite id conflict, keep it
                            if old_met_id not in met_sources_dict: # first reaction with conflict
                                met_sources_dict[old_met_id] = {model_index}
                                merged_model_metabolites.append(metabolite)
                        else: # substitute metabolite in reaction
                            st_coeff = reaction.metabolites[metabolite]
                            reaction.add_metabolites({old_met_id: -st_coeff})
                            reaction.add_metabolites({new_met_id: st_coeff})
                else:
                    metabolite.id = new_met_id
                    met_sources_dict[new_met_id].add(model_index)
            else:
                metabolite.id = new_met_id
                met_sources_dict[new_met_id] = {model_index}
                merged_model_metabolites.append(metabolite)
                if new_met_id in met_model_id_dict:
                    met_model_id_dict[new_met_id].append(old_met_id)
                else:
                    met_model_id_dict[new_met_id] = [old_met_id]

        for reaction in model.reactions:
            reac_id = reaction.id
            if reac_id in str(model.objective): # processing objective reactions
                model_objectives.append(reaction)
            else:
                reaction_key, rev_reaction_key = create_reaction_key(reaction)
                if reaction_key in merged_model_reactions_dict:
                    existing_reac_id = merged_model_reactions_dict[reaction_key]
                    reac_sources_dict[existing_reac_id].add(model_index)
                    if reac_id in reac_sources_dict:
                        reac_sources_dict[reac_id].add(model_index)
                elif rev_reaction_key in merged_model_reactions_dict:
                    rev_existing_reac_id = merged_model_reactions_dict[rev_reaction_key]
                    reac_sources_dict[rev_existing_reac_id].add(model_index)
                    if reac_id in reac_sources_dict:
                        reac_sources_dict[reac_id].add(model_index)
                else:
                    if reac_id in reac_sources_dict:
                        reac_id += '~'
                        while reac_id in reac_sources_dict:
                            reac_id += '~'
                        reaction.id = reac_id

                    merged_model_reactions.append(reaction)
                    merged_model_reactions_dict[reaction_key] = reac_id
                    reac_sources_dict[reac_id] = {model_index}

        objective_reactions.append(model_objectives)

    merged_model_name += ' (mergem v' + __version._version + ')'

    merged_model = cobra.Model(merged_model_id, merged_model_name)
    merged_model.add_metabolites(merged_model_metabolites)
    merged_model.add_reactions(merged_model_reactions)

    jacc_matrix = compute_jaccard_matrix(len(models), met_sources_dict, reac_sources_dict)

    num_mets_merged = sum([len(model.metabolites) for model in models]) - len(met_sources_dict)
    num_reacs_merged = sum([len(model.reactions) for model in models]) - len(reac_sources_dict)
    num_obj_reactions = sum([len(model_objectives) for model_objectives in objective_reactions])
    if num_obj_reactions:
        num_reacs_merged -= num_obj_reactions
        if set_objective == 'merge':
            num_reacs_merged += 1

    merged_model, reac_sources_dict = set_objective_expression(merged_model, reac_sources_dict, models,
                                                              objective_reactions, set_objective)

    merged_model.repair()

    # Revert mergem ids and convert sources to lists
    met_sources_dict = {met_model_id_dict.get(met_id, [met_id])[0]: list(sources) for (met_id, sources) in met_sources_dict.items()}
    reac_sources_dict = {reac_id : list(sources) for (reac_id, sources) in reac_sources_dict.items()}
    for metabolite in merged_model.metabolites:
        old_met_id = met_model_id_dict.get(metabolite.id, [metabolite.id])[0]
        if old_met_id != metabolite.id:
            if old_met_id in merged_model.metabolites:
                alt_old_met_id = old_met_id
                while alt_old_met_id in merged_model.metabolites:
                    if (idx := alt_old_met_id.rfind('@')) > 0 or (idx := alt_old_met_id.rfind('_')) > 0:
                        alt_old_met_id = alt_old_met_id[:idx] + '~' + alt_old_met_id[idx:]
                    else:
                        alt_old_met_id += '~'
                met_sources_dict[alt_old_met_id] = met_sources_dict[old_met_id]
                old_met_id = alt_old_met_id
            metabolite.id = old_met_id

    merged_model.repair()

    result = {}
    result['merged_model'] = merged_model
    result['jacc_matrix'] = jacc_matrix
    result['num_met_merged'] = num_mets_merged
    result['num_reac_merged'] = num_reacs_merged
    result['met_sources'] = met_sources_dict
    result['reac_sources'] = reac_sources_dict

    return result


# returns a metabolite id in mergem namespace with cellular localization
def map_metabolite_to_mergem_id(metabolite):
    """
    Takes a metabolite object as input and returns mergem_id notation for metabolite
    :param metabolite: Cobra metabolite object
    :return: mergem_id notation for the metabolite or None if there is no mapping
    """
    met_id = metabolite.id.lstrip('_') # remove leading underscores

    met_univ_id = __modelHandling.met_univ_id_dict.get(met_id)

    split = None
    if met_univ_id is None: # Re-check mapping without compartment
        if '@' in met_id:
            split = met_id.rsplit('@', 1)
        else:
            split = met_id.rsplit("_", 1)
        met_univ_id = __modelHandling.met_univ_id_dict.get(split[0])
        if met_univ_id is None:
            return None

    met_compartment = __modelHandling.map_localization(metabolite.compartment)
    if met_compartment == '':
        if split is None:
            if '@' in met_id:
                split = met_id.rsplit('@', 1)
            else:
                split = met_id.rsplit("_", 1)
        if len(split) > 1:
            met_compartment = __modelHandling.map_localization(split[1])
            if met_compartment == '':
                met_compartment = split[1]

    return "mergem_" + str(met_univ_id) + "_" + met_compartment


# reaction key is a frozenset of tuples of participating mets with their stoichiometric coeffs
def create_reaction_key(reaction):
    """
    Takes a reaction object as input and creates a key(frozen set) of all pairs of metabolite ID and stoichiometric
    coefficients. \n
    :param reaction: Cobra reaction object
    :return: frozen set of pairs of IDs of participating metabolite and their stoichiometric coefficients
    """
    reac_metabolite_set = set()
    reac_rev_met_set = set()
    for reactant in reaction.reactants:
        if (not reactant.id.startswith(__modelHandling.proton_mergem_id)) and (reactant.id[-1] != 'b') and (reactant.name != "PMF"):
            id = reactant.id if reactant.id.startswith('mergem_') else map_metabolite_to_mergem_id(reactant)
            metabolite_set = (id, -1)
            rev_met_set = (id, 1)
            reac_metabolite_set.add(metabolite_set)
            reac_rev_met_set.add(rev_met_set)

    for product in reaction.products:
        if (not product.id.startswith(__modelHandling.proton_mergem_id)) and (product.id[-1] != 'b') and (product.name != "PMF"):
            id = product.id if product.id.startswith('mergem_') else map_metabolite_to_mergem_id(product)
            metabolite_set = (id, 1)
            rev_met_set = (id, -1)
            reac_metabolite_set.add(metabolite_set)
            reac_rev_met_set.add(rev_met_set)

    reac_metabolite_set = frozenset(reac_metabolite_set)
    reac_rev_met_set = frozenset(reac_rev_met_set)

    return reac_metabolite_set, reac_rev_met_set


# sets the objective for merged model
def set_objective_expression(merged_model, reac_sources_dict, models, objective_reactions, set_objective):
    """
    Sets the objective expression for merged model.
    :param merged_model: Model whose objective is to be set.
    :param models: List of all input models from which objective is chosen.
    :param objective_reactions: List of (lists of) all objective reactions in input models.
    :param set_objective: 'merge' all input model objective reacs or from one of the input models.
    :return: model with its objective expression set.
    """
    if set_objective == 'merge':
        merged_model, reac_sources_dict = create_merged_objective(merged_model, reac_sources_dict, objective_reactions)
    else:
        model_num = int(set_objective) - 1
        reactions = objective_reactions[model_num]
        if len(reactions):
            merged_model.add_reactions(reactions)
            merged_model.objective = models[model_num].objective.expression
            for reaction in reactions:
                reac_sources_dict[reaction.id] = {model_num}

    return merged_model, reac_sources_dict


# create a reaction that contains input model objective reacs merged together
def create_merged_objective(merged_model, reac_sources_dict, objective_reactions):
    """
    Merges objective reactions in list and sets as objective in input merged model.
    :param merged_model: Merged model to set objective of
    :param objective_reactions: list of objective reaction lists to merge
    :return: Merged model with its objective set to the merged objective reactions
    """
    merged_obj_reaction_id = 'merged-objectives'
    merged_obj_reaction_name = 'Merged all objectives'
    st_dict = {}
    metabolite_dict = {}

    for reaction_list in objective_reactions:
        for reaction in reaction_list:
            for metabolite in reaction.metabolites:
                if metabolite.id not in st_dict:
                    st_dict[metabolite.id] = set()
                    metabolite_dict[metabolite.id] = metabolite

                st_dict[metabolite.id].add(reaction.metabolites[metabolite])

            merged_obj_reaction_name += '; ' + reaction.name

    if len(st_dict):
        merged_obj_reaction = cobra.Reaction(merged_obj_reaction_id, merged_obj_reaction_name)

        for metabolite_id, met_stoichiometries in st_dict.items():
            avg_stoichiometry = sum(met_stoichiometries)/len(met_stoichiometries)
            merged_obj_reaction.add_metabolites({metabolite_dict[metabolite_id]: avg_stoichiometry})

        merged_model.add_reactions([merged_obj_reaction])
        merged_model.objective = merged_obj_reaction_id

        reac_sources_dict[merged_obj_reaction_id] = list(range(0, len(objective_reactions)))

    return merged_model, reac_sources_dict


def compute_jaccard_matrix(num_models, met_source_dict, reac_source_dict):
    """
    Creates Jaccard distances matrix using dictionaries with source of each metabolite
    and reaction in merged model.
    :param num_models: number of input models
    :param met_source_dict: dictionary with source of each metabolite
    :param reac_source_dict: dictionary with source of each reaction
    :return: Jaccard distance matrices for metabolites and reactions
    """
    jacc_matrix = []

    model_met_ids = []
    model_reac_ids = []
    for i in range(0, num_models):
        model_met_ids.append({id for id,model_set in met_source_dict.items() if i in model_set})
        model_reac_ids.append({id for id,model_set in reac_source_dict.items() if i in model_set})

    for i in range(0, num_models):
        jd_row = []
        for j in range(0, num_models):
            if i > j:
                union = len(model_reac_ids[i] | model_reac_ids[j])
                if union == 0:
                    jd_row += [0]
                else:
                    jd_row +=  [1 - (len(model_reac_ids[i] & model_reac_ids[j]) / union)]

            elif i < j:
                union = len(model_met_ids[i] | model_met_ids[j])
                if union == 0:
                    jd_row += [0]
                else:
                    jd_row +=  [1 - (len(model_met_ids[i] & model_met_ids[j]) / union)]

            else:
                jd_row += [0]

        jacc_matrix += [jd_row]

    return jacc_matrix
