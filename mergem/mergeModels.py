"""
    Uses Fluxer dictionaries to translate models to common namespace
    and then merges the models to remove metabolite and reaction
    duplicates.
    Resulting model can be FBA simulated.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""
from . import __modelHandling
import cobra


# update metabolite id mapping pickle
def update_id_mapper():
    """
    Calls the update id mapping pickles function, which downloads most recent database identifier and
    cross-reference files and merges identifiers based on common properties.
    """
    __modelHandling.__update_id_mapping_pickles()


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
    fluxer_ids_created = {}
    model_jaccard_distances = []
    ref_model_mets = set()
    met_and_reac_sources = []

    merged_model_id = 'merged_'
    template_model = cobra.Model(merged_model_id)
    num_ref_model_reacs = len(list_of_models[0].reactions)
    __modelHandling.__load_or_create_id_mapper()
    for model_index in range(0, len(list_of_models)):
        model = list_of_models[model_index]
        merged_model_id += model.id + "_"
        list_model_objectives = []
        num_met_merged = 0
        num_reac_merged = 0

        for metabolite in model.metabolites:
            met_fluxer_id = __modelHandling.__create_fluxer_metabolite_id(metabolite)
            met_copy = metabolite.copy()

            if met_fluxer_id is None:
                curr_id = metabolite.id
            else:
                fluxer_ids_created[metabolite.id] = met_fluxer_id
                curr_id = met_fluxer_id

            if curr_id in met_model_id_dict:
                met_model_id_dict[curr_id] += [metabolite.id]
                met_source_dict[curr_id] |= {model_index}
            else:
                met_model_id_dict[curr_id] = [metabolite.id]
                met_copy.id = curr_id
                template_model.add_metabolites(met_copy)
                met_source_dict[curr_id] = {model_index}

            if model_index == 0:
                ref_model_mets |= {curr_id}
            else:
                if curr_id in ref_model_mets:
                    num_met_merged += 1
                    total_mets_merged += 1

        for reaction in model.reactions:
            for reaction_met in reaction.metabolites:
                new_met_id = fluxer_ids_created.get(reaction_met.id)
                if (new_met_id not in model.metabolites) and (new_met_id is not None):
                    reaction_met.id = new_met_id

                elif new_met_id in model.metabolites:
                    st_coeff = reaction.metabolites[reaction_met]
                    reaction.add_metabolites({new_met_id: st_coeff})
                    reaction.add_metabolites({reaction_met: 0})

            if reaction.id in str(model.objective):
                list_model_objectives.append(reaction)
                if reaction.id in reac_source_dict:
                    reac_source_dict[reaction.id] |= {model_index}
                else:
                    reac_source_dict[reaction.id] = {model_index}

            else:
                reaction_key = __modelHandling.__create_reaction_key(reaction, False)
                if reaction_key not in merged_model_reactions:
                    template_model.add_reactions([reaction.copy()])
                    merged_model_reactions[reaction_key] = reaction.copy()

                    if reaction.id not in reac_source_dict:
                        reac_source_dict[reaction.id] = {model_index}
                    else:
                        reac_source_dict[reaction.id] |= {model_index}

                    if model_index == 0:
                        num_ref_model_reacs += 1

                else:
                    existing_reac = merged_model_reactions.get(reaction_key)
                    reac_source_dict[existing_reac.id] |= {model_index}  # source of reactions by metabolite key

                    if model_index != 0:
                        if 0 in reac_source_dict[existing_reac.id]:
                            num_reac_merged += 1
                            total_reacs_merged += 1

                    existing_reac.gene_reaction_rule = __modelHandling.__update_gene_rule(existing_reac, reaction)

        if model_index == 0:
            model_jaccard_distances += [(0, 0)]
        else:
            model_jaccard_distances += [(__modelHandling.__calculate_jaccard_distance(len(ref_model_mets), num_met_merged),
                                        __modelHandling.__calculate_jaccard_distance(num_ref_model_reacs, num_reac_merged))]

        list_of_objective_reactions.append(list_model_objectives)

    met_source_cleaned = {}
    for metabolite_id in met_source_dict.keys():
        source_set = met_source_dict[metabolite_id]
        if "fluxer" in metabolite_id:  # Only revert fluxer IDs
            met_source_cleaned[met_model_id_dict[metabolite_id][0]] = source_set
        else:
            met_source_cleaned[metabolite_id] = source_set

    met_and_reac_sources += (__modelHandling.__set_to_list_source_dict(met_source_cleaned),
                             __modelHandling.__set_to_list_source_dict(reac_source_dict))
    template_model = __modelHandling.__set_objective_expression(template_model, list_of_models,
                                                                list_of_objective_reactions, set_objective)
    merged_model_id = merged_model_id[:-1]
    merged_model = __modelHandling.__convert_template_to_merged_model(template_model, merged_model_id, met_model_id_dict)

    result['merged_model'] = merged_model
    result['jacc_d'] = model_jaccard_distances
    result['Met_merged'] = total_mets_merged
    result['Reac_merged'] = total_reacs_merged
    result['Met_sources'] = met_and_reac_sources[0]
    result['Reac_sources'] = met_and_reac_sources[1]

    return result
