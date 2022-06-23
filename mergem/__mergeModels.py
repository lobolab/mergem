"""
    Uses mergem dictionaries to translate models to common namespace
    and then merges the models to remove metabolite and reaction
    duplicates.
    Resulting model can be FBA simulated.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""
from . import __modelHandling
from .__database_id_merger import __return_mapping_and_info_dicts
import cobra


# update metabolite id mapping pickle
def update_id_mapper():
    """
    Calls the update id mapping pickles function, which downloads most recent database identifier and
    cross-reference files and merges identifiers based on common properties.
    """
    __modelHandling.__update_id_mapping_pickles()


# loads reaction and metabolite ID mapping and info dictionaries from pickles
def get_mapping_and_info_dicts():
    """
    Returns dictionaries loaded from pickles. New ID mapping and info dictionaries and pickles are created
    if they do not exist in data folder.
    """
    reaction_id_mapper, reaction_info, metabolite_id_mapper, metabolite_info = __return_mapping_and_info_dicts()
    return reaction_id_mapper, reaction_info, metabolite_id_mapper, metabolite_info


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
    merged_model_name = 'Mergem'
    merged_model_reactions = []

    # Add first model
    model = models[0]
    merged_model_id += '-' + model.id
    merged_model_name += '; ' + (model.name if model.name else model.id)
    model_objectives = []

    for metabolite in model.metabolites:
        old_met_id = metabolite.id
        new_met_id = __modelHandling.map_to_metabolite_mergem_id(metabolite)

        if (new_met_id is None) or (new_met_id in met_sources_dict):
            met_sources_dict[old_met_id] = {0}
        else:
            met_sources_dict[new_met_id] = {0}
            met_model_id_dict[new_met_id] = old_met_id
            metabolite.id = new_met_id

    for reaction in model.reactions:
        reac_id = reaction.id
        if reac_id in str(model.objective): # processing objective reactions
            model_objectives.append(reaction)
        else:
            reaction_key, rev_reaction_key = __modelHandling.__create_reaction_key(reaction, False)
            merged_model_reactions.append(reaction)
            merged_model_reactions_dict[reaction_key] = reac_id
            merged_model_reactions_dict[rev_reaction_key] = reac_id
            reac_sources_dict[reac_id] = {0}

    objective_reactions.append(model_objectives)

    # Merge rest of models
    for model_index in range(1, len(models)):
        model = models[model_index]
        merged_model_id += '-' + model.id
        merged_model_name += '; ' + (model.name if model.name else model.id)
        model_objectives = []

        for metabolite in model.metabolites:
            old_met_id = metabolite.id

            if old_met_id in met_sources_dict:
                met_sources_dict[old_met_id].add(model_index)
            else:
                new_met_id = __modelHandling.map_to_metabolite_mergem_id(metabolite)

                if new_met_id is None:
                    met_sources_dict[old_met_id] = {model_index}
                else:
                    if new_met_id in met_sources_dict:
                        if model_index in met_sources_dict[new_met_id]:
                            met_sources_dict[old_met_id] = {model_index}
                        else:
                            metabolite.id = new_met_id
                            met_sources_dict[new_met_id].add(model_index)
                    else:
                        metabolite.id = new_met_id
                        met_sources_dict[new_met_id] = {model_index}
                        met_model_id_dict[new_met_id] = old_met_id

        for reaction in model.reactions:
            reac_id = reaction.id
            if reac_id in str(model.objective): # processing objective reactions
                model_objectives.append(reaction)
            else:
                reaction_key, rev_reaction_key = __modelHandling.__create_reaction_key(reaction, False)
                if (reaction_key in merged_model_reactions_dict):
                    existing_reac_id = merged_model_reactions_dict[reaction_key]
                    reac_sources_dict[existing_reac_id].add(model_index)
                elif (rev_reaction_key in merged_model_reactions_dict):
                    rev_existing_reac_id = merged_model_reactions_dict[rev_reaction_key]
                    reac_sources_dict[rev_existing_reac_id].add(model_index)
                else:
                    merged_model_reactions.append(reaction)
                    merged_model_reactions_dict[reaction_key] = reac_id
                    merged_model_reactions_dict[rev_reaction_key] = reac_id
                    reac_sources_dict[reac_id] = {model_index}

        objective_reactions.append(model_objectives)

    merged_model = cobra.Model(merged_model_id, merged_model_name)
    merged_model.add_reactions(merged_model_reactions)

    jacc_matrix = __compute_jaccard_matrix(len(models), met_sources_dict, reac_sources_dict)

    num_mets_merged = sum([len(model.metabolites) for model in models]) - len(met_sources_dict)
    num_reacs_merged = sum([len(model.reactions) for model in models]) - len(reac_sources_dict)

    merged_model, reac_sources_dict = __modelHandling.__set_objective_expression(merged_model, reac_sources_dict, models,
                                                              objective_reactions, set_objective)

    # Revert mergem ids and convert sources to lists
    met_sources_dict = {met_model_id_dict.get(met_id, met_id): list(sources) for (met_id, sources) in met_sources_dict.items()}
    reac_sources_dict = {reac_id : list(sources) for (reac_id, sources) in reac_sources_dict.items()}
    for metabolite in merged_model.metabolites:
        metabolite.id = met_model_id_dict.get(metabolite.id, metabolite.id)

    merged_model.repair()

    result = {}
    result['merged_model'] = merged_model
    result['jacc_matrix'] = jacc_matrix
    result['num_met_merged'] = num_mets_merged
    result['num_reac_merged'] = num_reacs_merged
    result['met_sources'] = met_sources_dict
    result['reac_sources'] = reac_sources_dict

    return result


def __compute_jaccard_matrix(num_models, met_source_dict, reac_source_dict):
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
