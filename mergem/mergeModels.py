"""
    Uses Fluxer dictionaries to translate models to common namespace
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
def merge(list_of_inputs, set_objective='merge'):
    """
    Takes a list of cobra models as input and merges them into a single model with the chosen objective. \n
    :param list_of_inputs: list of cobra models
    :param set_objective: objective reaction from one of the models or merge (default) all model objectives
    :return: A dictionary of the merged model, met & reac jaccard distances, num of mets and reacs merged,
            and met & reac sources.
    """
    list_of_models = []

    for mod_input in list_of_inputs:
        if isinstance(mod_input, str):
            try:
                input_model = __modelHandling.__load_model(mod_input)
                list_of_models.append(input_model)

            except Exception as e:
                print("Error loading model: ", e)

        else:
            list_of_models.append(mod_input)

    result = {}
    list_of_objective_reactions = []
    met_source_dict, reac_source_dict, met_model_id_dict = {}, {}, {}
    merged_model_reactions = {}
    total_mets_merged, total_reacs_merged = 0, 0
    met_and_reac_sources = []

    merged_model_id = 'merged_'
    template_model = cobra.Model(merged_model_id)
    __modelHandling.__load_or_create_id_mapper()
    cleaned_models = []

    # model 'cleanup'
    for model in list_of_models:
        model_mets = {}
        update_reactions = []
        for met in model.metabolites:
            met_fluxer_id = __modelHandling.create_fluxer_metabolite_id(met)
            if met_fluxer_id not in model_mets:
                model_mets[met_fluxer_id] = met.copy()
            else:
                existing_met = model.metabolites.get_by_id(model_mets[met_fluxer_id].id)
                if existing_met.reactions != met.reactions:
                    for reaction in met.reactions:
                        if met in reaction.metabolites:
                            update_reactions += [(reaction.id, existing_met.id, met.id)]

        for reaction_id, new_met_id, old_met_id in update_reactions:
            reaction = model.reactions.get_by_id(reaction_id)
            reaction_copy = reaction.copy()
            for met in reaction_copy.metabolites:
                if met.id == old_met_id:
                    st_coeff = reaction_copy.metabolites[met]
                    reaction.add_metabolites({new_met_id: st_coeff})
                    reaction.add_metabolites({met.id: -st_coeff})

        model.repair()
        cleaned_models.append(model)

    for model_index in range(0, len(cleaned_models)):
        model = list_of_models[model_index]
        merged_model_id += model.id + "_"
        list_model_objectives = []
        met_counted = set()

        for reaction in model.reactions:
            if reaction.id in str(model.objective): # processing objective reactions
                list_model_objectives.append(reaction)
                if reaction.id in reac_source_dict:
                    reac_source_dict[reaction.id] |= {model_index}
                else:
                    reac_source_dict[reaction.id] = {model_index}

            else:
                for reaction_met in reaction.metabolites: # processing metabolic reactions that are not in objective
                    new_met_id = __modelHandling.create_fluxer_metabolite_id(reaction_met)

                    if new_met_id is not None:
                        if new_met_id in met_source_dict:
                            met_source_dict[new_met_id] |= {model_index}
                            met_model_id_dict[new_met_id] += [reaction_met.id]
                        else:
                            met_source_dict[new_met_id] = {model_index}
                            met_model_id_dict[new_met_id] = [reaction_met.id]
                        met_counted |= {new_met_id}
                    else:
                        if reaction_met.id in met_source_dict:
                            met_source_dict[reaction_met.id] |= {model_index}
                            met_counted |= {reaction_met.id}
                        else:
                            met_source_dict[reaction_met.id] = {model_index}

                    if (new_met_id not in model.metabolites) and (new_met_id is not None):
                        reaction_met.id = new_met_id

                reaction_key, rev_reaction_key = __modelHandling.__create_reaction_key(reaction, False)
                if (reaction_key not in merged_model_reactions) and (rev_reaction_key not in merged_model_reactions):
                    template_model.add_reactions([reaction.copy()])
                    merged_model_reactions[reaction_key] = reaction.copy()
                    merged_model_reactions[rev_reaction_key] = reaction.copy()

                    reac_source_dict[reaction.id] = {model_index}

                else:
                    existing_reac = merged_model_reactions.get(reaction_key)
                    if existing_reac is not None:
                        reac_source_dict[existing_reac.id] |= {model_index}

                    else:
                        rev_existing_reac = merged_model_reactions.get(rev_reaction_key)
                        if rev_existing_reac is not None:
                            reac_source_dict[rev_existing_reac.id] |= {model_index}

        list_of_objective_reactions.append(list_model_objectives)

    jacc_matrix_met, jacc_matrix_reac = __modelHandling.__create_jaccard_matrix(len(cleaned_models), met_source_dict,
                                                                                reac_source_dict)

    met_source_cleaned = {}
    for metabolite_id in met_source_dict.keys():
        source_set = met_source_dict[metabolite_id]
        if "fluxer" in metabolite_id:  # Only revert fluxer IDs
            met_source_cleaned[met_model_id_dict[metabolite_id][0]] = source_set
        else:
            met_source_cleaned[metabolite_id] = source_set

    if set_objective == 'merge':
        reac_source_dict['merged_objectives'] = set(model_num for model_num in range(0, len(list_of_models)))

    met_and_reac_sources += (__modelHandling.__set_to_list_source_dict(met_source_cleaned),
                             __modelHandling.__set_to_list_source_dict(reac_source_dict))
    template_model = __modelHandling.__set_objective_expression(template_model, list_of_models,
                                                                list_of_objective_reactions, set_objective)
    merged_model_id = merged_model_id[:-1]
    merged_model = __modelHandling.__convert_template_to_merged_model(template_model, merged_model_id,
                                                                      met_model_id_dict)

    merged_model.repair()

    result['merged_model'] = merged_model
    result['jacc_d'] = [jacc_matrix_met, jacc_matrix_reac]
    result['Met_merged'] = total_mets_merged
    result['Reac_merged'] = total_reacs_merged
    result['Met_sources'] = met_and_reac_sources[0]
    result['Reac_sources'] = met_and_reac_sources[1]

    return result
