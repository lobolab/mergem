"""
    Uses mergem dictionaries to translate models to common namespace
    and then merges the models to remove metabolite and reaction
    duplicates.
    Resulting model can be FBA simulated.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""

from . import __model_handling
from . import __version
import cobra
from collections import defaultdict

# translate all metabolite and reaction IDs to a target namespace
def translate(input_model, trans_to_db=None):
    """
    Translates metabolite and reaction IDs to a target namespace
    :param input_model: a cobra model or file name
    :param trans_to_db: target database to be translated to
    :return: model with translated metabolite and reaction IDs.
    """
    return merge([input_model], trans_to_db=trans_to_db)['merged_model']


# merges models in a list to the template/first model
# set_objective can be an integer for model obj or 'merge'
def merge(input_models, set_objective='merge', exact_sto=False, use_prot=False, extend_annot=False, trans_to_db=None):
    """
    Takes a list of cobra models or file names as input and merges them into a single model with the chosen objective. \n
    :param input_models: list of cobr+a models or file names
    :param set_objective: objective reaction from one of the models or merge (default) all model objectives
    :param exact_sto: Boolean which determines whether exact stoichiometry of metabolites is used during merging
    :param use_prot: Boolean to consider hydrogen and proton when merging reactions
    :param add_annot: Boolean to add additional metabolite and reaction annotations from mergem dictionaries
    :param trans_to_db: target database to be translated to
    :return: a dictionary of the merged model, met & reac jaccard distances, num of mets and reacs merged,
            and met & reac sources.
    """
    __model_handling.load_met_univ_id_dict()
    models = []

    for input_model in input_models:
        if isinstance(input_model, str):
            try:
                input_model = __model_handling.load_model(input_model)
                models.append(input_model)

            except Exception as e:
                print("Error loading model: ", e)

        else:
            models.append(input_model)

    objective_reactions = []
    met_model_id_dict, met_sources_dict, merged_model_reactions_dict = {}, {}, {}
    met_sources_dict = defaultdict(lambda:defaultdict(list))
    reac_sources_dict = defaultdict(lambda:defaultdict(list))

    merged_model_id = 'mergem'
    merged_model_name = 'Mergem of '
    merged_model_metabolites = []
    merged_model_reactions = []

    # Add first model
    model = models[0]
    merged_model_id += '_' + model.id
    merged_model_name += model.name if model.name else model.id
    merged_compartments = model.compartments
    model_objectives = []

    dict_met_annot, dict_reac_annot, dict_gprs = {}, {}, {}

    for metabolite in model.metabolites:
        merged_model_metabolites.append(metabolite)
        old_met_id = metabolite.id
        new_met_id = map_metabolite_to_mergem_id(metabolite)

        if (new_met_id is None) or (new_met_id in met_sources_dict):
            met_sources_dict[old_met_id][0].append(old_met_id)
            dict_met_annot[old_met_id] = metabolite.annotation
        else:
            met_sources_dict[new_met_id][0].append(old_met_id)
            dict_met_annot[new_met_id] = metabolite.annotation
            if new_met_id in met_model_id_dict:
                met_model_id_dict[new_met_id].append(old_met_id)
            else:
                met_model_id_dict[new_met_id] = [old_met_id]
            metabolite.id = new_met_id

    for reaction in model.reactions:
        reac_id = reaction.id
        if reac_id in str(model.objective):  # processing objective reactions
            model_objectives.append(reaction)
        else:
            merged_model_reactions.append(reaction)
            reac_sources_dict[reac_id][0].append(reac_id)
            reaction_key, rev_reaction_key = create_reaction_key(reaction, exact_sto, use_prot)
            if not reaction_key in merged_model_reactions_dict:
                merged_model_reactions_dict[reaction_key] = reac_id

            dict_reac_annot[reac_id] = reaction.annotation
            dict_gprs[reac_id] = reaction.gpr
    objective_reactions.append(model_objectives)

    # Merge rest of models
    for model_index in range(1, len(models)):
        model = models[model_index]
        merged_model_id += '_' + model.id
        merged_model_name += '; ' + (model.name if model.name else model.id)
        merged_compartments = model.compartments | merged_compartments
        model_objectives = []
        model_metabolite_ids = {m.id for m in model.metabolites}

        for metabolite in model.metabolites:
            old_met_id = metabolite.id
            new_met_id = map_metabolite_to_mergem_id(metabolite)

            if new_met_id is None:
                if old_met_id in met_sources_dict:
                    met_sources_dict[old_met_id][model_index].append(old_met_id)
                    __model_handling.add_annotations(old_met_id, dict_met_annot, metabolite)
                else:
                    met_sources_dict[old_met_id][model_index].append(old_met_id)
                    merged_model_metabolites.append(metabolite)
                    dict_met_annot[old_met_id] = metabolite.annotation

            elif old_met_id in met_sources_dict:  # priority is given to original metabolite ids
                met_sources_dict[old_met_id][model_index].append(old_met_id)
                __model_handling.add_annotations(old_met_id, dict_met_annot, metabolite)

            elif new_met_id in met_sources_dict:  # new metabolite id previously found
                old_met_ids = met_model_id_dict[new_met_id]
                if (old_met_id not in old_met_ids) and any(id in old_met_ids for id in model_metabolite_ids):  # model has a better match
                    met_sources_dict[old_met_id][model_index].append(old_met_id)
                    merged_model_metabolites.append(metabolite)
                    dict_met_annot[old_met_id] = metabolite.annotation

                elif model_index in met_sources_dict[new_met_id]:  # model already had a metabolite for this mergem id
                    for reaction in metabolite.reactions:  # replace id in its reactions
                        if new_met_id in [met.id for met in reaction.metabolites]: # new metabolite id conflict, keep it
                            if old_met_id not in met_sources_dict:  # first reaction with conflict
                                met_sources_dict[old_met_id][model_index].append(old_met_id)
                                merged_model_metabolites.append(metabolite)
                                dict_met_annot[old_met_id] = metabolite.annotation

                        else:  # substitute metabolite in reaction
                            st_coeff = reaction.metabolites[metabolite]
                            reaction.add_metabolites({old_met_id: -st_coeff})
                            reaction.add_metabolites({new_met_id: st_coeff})
                else:
                    met_sources_dict[new_met_id][model_index].append(old_met_id)
                    metabolite.id = new_met_id
                    __model_handling.add_annotations(new_met_id, dict_met_annot, metabolite)
            else:
                metabolite.id = new_met_id
                met_sources_dict[new_met_id][model_index].append(old_met_id)
                merged_model_metabolites.append(metabolite)
                dict_met_annot[new_met_id] = metabolite.annotation
                if new_met_id in met_model_id_dict:
                    met_model_id_dict[new_met_id].append(old_met_id)
                else:
                    met_model_id_dict[new_met_id] = [old_met_id]

        for reaction in model.reactions:
            reac_id = reaction.id
            if reac_id in str(model.objective):  # processing objective reactions
                model_objectives.append(reaction)
            else:
                reaction_key, rev_reaction_key = create_reaction_key(reaction, exact_sto, use_prot)
                if reaction_key in merged_model_reactions_dict:
                    existing_reac_id = merged_model_reactions_dict[reaction_key]
                    reac_sources_dict[existing_reac_id][model_index].append(reac_id)
                    __model_handling.add_annotations(existing_reac_id, dict_reac_annot, reaction)
                    __model_handling.add_gpr(existing_reac_id, dict_gprs, reaction)
                    if reac_id in reac_sources_dict:
                        if reac_id != existing_reac_id:
                            reac_sources_dict[reac_id][model_index].append(reac_id)
                        __model_handling.add_annotations(reac_id, dict_reac_annot, reaction)
                        __model_handling.add_gpr(reac_id, dict_gprs, reaction)

                elif rev_reaction_key in merged_model_reactions_dict:
                    rev_existing_reac_id = merged_model_reactions_dict[rev_reaction_key]
                    reac_sources_dict[rev_existing_reac_id][model_index].append(reac_id)
                    __model_handling.add_annotations(rev_existing_reac_id, dict_reac_annot, reaction)
                    __model_handling.add_gpr(rev_existing_reac_id, dict_gprs, reaction)
                    if reac_id in reac_sources_dict:
                        if reac_id != rev_existing_reac_id:
                            reac_sources_dict[reac_id][model_index].append(reac_id)
                        __model_handling.add_annotations(reac_id, dict_reac_annot, reaction)
                        __model_handling.add_gpr(reac_id, dict_gprs, reaction)
                else:
                    orig_reac_id = reac_id
                    if reac_id in reac_sources_dict:
                        reac_id += '~'
                        while reac_id in reac_sources_dict:
                            reac_id += '~'
                        reaction.id = reac_id

                    merged_model_reactions.append(reaction)
                    merged_model_reactions_dict[reaction_key] = reac_id
                    reac_sources_dict[reac_id][model_index].append(orig_reac_id)
                    dict_reac_annot[reac_id] = reaction.annotation
                    dict_gprs[reac_id] = reaction.gpr

        objective_reactions.append(model_objectives)

    if exact_sto:
        merged_model_id += '_exactsto'
        merged_model_name += ' with exact stoichiometry'

    if use_prot:
        merged_model_id += '_useprot'
        merged_model_name += ' with protonation'

    if trans_to_db:
        merged_model_id += '_trans_' + trans_to_db
        merged_model_name += ' translated to ' + trans_to_db

    merged_model_name += ' (mergem v' + __version._version + ')'

    merged_model = cobra.Model(merged_model_id, merged_model_name)
    merged_model.add_metabolites(merged_model_metabolites)
    merged_model.add_reactions(merged_model_reactions)
    merged_model.compartments = {(k,v) for k,v in merged_compartments.items() if k in merged_model.compartments}

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

    # Convert IDs back to originals
    met_sources_dict = {met_model_id_dict.get(met_id, [met_id])[0]: sources 
                        for (met_id, sources) in met_sources_dict.items()}
    reac_sources_dict = {reac_id: sources 
                         for (reac_id, sources) in reac_sources_dict.items()}
    
    if trans_to_db:
        trans_to_db += ':'

    # Post-processing metabolites
    for metabolite in merged_model.metabolites:
        metabolite.annotation = dict_met_annot.get(metabolite.id, {})
        old_met_id = None
        if trans_to_db or extend_annot:
            met_id_array = metabolite.id.split("_")
            met_univ_id = int(met_id_array[1]) if met_id_array[0] == "mergem" else __model_handling.map_metabolite_univ_id(metabolite.id)

            if met_univ_id:
                met_props = __model_handling.get_metabolite_properties(met_univ_id)

                if extend_annot:
                    __model_handling.extend_metabolite_annotations(metabolite, met_props)

                if trans_to_db:
                    trans_met_id = next(
                        ((s[len(trans_to_db):] + (('_' + met_id_array[-1]) if len(met_id_array) > 1 else '')) for s in met_props['ids'] if s.startswith(trans_to_db)), 
                        None)
                    if trans_met_id:
                        met_sources_dict[trans_met_id] = met_sources_dict[met_model_id_dict.get(metabolite.id, [metabolite.id])[0]]
                        old_met_id = trans_met_id
        
        if not old_met_id:
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

    # Post-processing reactions
    for reaction in merged_model.reactions:
        reaction.annotation = dict_reac_annot.get(reaction.id, {})
        if reaction.id in dict_gprs: reaction.gpr = dict_gprs[reaction.id]

        if trans_to_db or extend_annot:
            if reac_mergem_id := __model_handling.map_reaction_univ_id(reaction.id):
                reac_props = __model_handling.get_reaction_properties(reac_mergem_id)

                if extend_annot:
                    __model_handling.extend_reaction_annotations(reaction, reac_props)

                if trans_to_db:
                    new_reac_id = next(
                        ((s[len(trans_to_db):]) for s in reac_props['ids'] if s.startswith(trans_to_db)),
                        None)

                    if new_reac_id:
                        while new_reac_id in merged_model.reactions:
                            new_reac_id += '~'

                        reac_sources_dict[new_reac_id] = reac_sources_dict[reaction.id]
                        reaction.id = new_reac_id
            elif trans_to_db: # Translate reactions with a metabolite id (e.g., exchange reactions)
                reac_id_array = reaction.id.split('_')
                met_univ_id = None if len(reac_id_array) < 2 else __model_handling.map_metabolite_univ_id('_'.join(reac_id_array[1:]))
                if met_univ_id:
                    met_props = __model_handling.get_metabolite_properties(met_univ_id)
                    new_reac_id = next(
                        ((reac_id_array[0] + '_' + s[len(trans_to_db):] + (('_' + reac_id_array[-1]) if len(reac_id_array) > 2 else '')) for s in met_props['ids'] if s.startswith(trans_to_db)), 
                        None)
                    if new_reac_id:
                        while new_reac_id in merged_model.reactions:
                            new_reac_id += '~'

                        reac_sources_dict[new_reac_id] = reac_sources_dict[reaction.id]
                        reaction.id = new_reac_id

    # Post-processing genes
    if extend_annot:
        for gene in merged_model.genes:
            gene.annotation['sbo'] = 'SBO:0000243'

    merged_model.repair()
    
    # Clean source dicts
    met_sources_dict = {m.id: {model_index:(ids if len(ids) > 1 else ids[0]) 
                               for model_index,ids in met_sources_dict[m.id].items()} 
                               for m in merged_model.metabolites}
    reac_sources_dict = {r.id: {model_index:(ids if len(ids) > 1 else ids[0]) 
                                for model_index,ids in reac_sources_dict[r.id].items()} 
                                for r in merged_model.reactions}

    results = {}
    results['merged_model'] = merged_model
    results['jacc_matrix'] = jacc_matrix
    results['num_met_merged'] = num_mets_merged
    results['num_reac_merged'] = num_reacs_merged
    results['met_sources'] = met_sources_dict
    results['reac_sources'] = reac_sources_dict

    return results


# returns a metabolite id in mergem namespace with cellular localization
def map_metabolite_to_mergem_id(metabolite):
    """
    Takes a metabolite object as input and returns mergem_id notation for metabolite
    :param metabolite: Cobra metabolite object
    :return: mergem_id notation for the metabolite or None if there is no mapping
    """
    met_id = metabolite.id.lstrip('_')  # remove leading underscores

    met_univ_id = __model_handling.met_univ_id_dict.get(met_id)

    split = None
    if met_univ_id is None:  # Re-check mapping without compartment
        if '@' in met_id:
            split = met_id.rsplit('@', 1)
        else:
            split = met_id.rsplit("_", 1)
        met_univ_id = __model_handling.met_univ_id_dict.get(split[0])

        if (met_univ_id is None) and ('mergem' not in met_id):  # no mapping for metabolite ID
            for annot in metabolite.annotation:
                if annot != 'sbo':
                    met_id_from_annot = metabolite.annotation[annot]  # get xref from annotation
                    if type(met_id_from_annot) == str:
                        met_univ_id = __model_handling.met_univ_id_dict.get(met_id_from_annot)
                    elif type(met_id_from_annot) == list:
                        for annot_met_id in met_id_from_annot:
                            if ':' in annot_met_id:
                                split_annot_id = annot_met_id.split(':', 1)[1]
                                met_univ_id = __model_handling.met_univ_id_dict.get(split_annot_id)
                            else:
                                met_univ_id = __model_handling.met_univ_id_dict.get(annot_met_id)
                            if met_univ_id:
                                break
                    if met_univ_id:
                        break
        if met_univ_id is None:
            return None

    met_compartment = __model_handling.map_localization(metabolite.compartment)
    if met_compartment == '':
        if split is None:
            if '@' in met_id:
                split = met_id.rsplit('@', 1)
            else:
                split = met_id.rsplit("_", 1)
        if len(split) > 1:
            met_compartment = __model_handling.map_localization(split[1])
            if met_compartment == '':
                met_compartment = split[1]

    return "mergem_" + str(met_univ_id) + "_" + met_compartment


# reaction key is a frozenset of tuples of participating mets with their stoichiometric coeffs
def create_reaction_key(reaction, exact_sto, use_prot):
    """
    Takes a reaction object as input and creates a key(frozen set) of all pairs of metabolite ID and stoichiometric
    coefficients. \n
    :param reaction: Cobra reaction object
    :param exact_sto: Reaction stoichiometric coefficient
    :param use_prot: Inclue hydrogen and protons
    :return: frozen set of pairs of IDs of participating metabolite and their stoichiometric coefficients
    """
    reac_metabolite_set = set()
    reac_rev_met_set = set()
    for reactant in reaction.reactants:
        if (not use_prot) and (reactant.id.startswith(__model_handling.proton_mergem_id) or reactant.name == "PMF"):
            continue

        elif reactant.id[-1] != 'b':
            id = reactant.id if reactant.id.startswith('mergem_') else map_metabolite_to_mergem_id(reactant)
            stoc = reaction.metabolites[reactant] if exact_sto else 1
            metabolite_set = (id, -stoc)
            rev_met_set = (id, stoc)
            reac_metabolite_set.add(metabolite_set)
            reac_rev_met_set.add(rev_met_set)

    for product in reaction.products:
        if (not use_prot) and (product.id.startswith(__model_handling.proton_mergem_id) or product.name == "PMF"):
            continue
        elif product.id[-1] != 'b':
            id = product.id if product.id.startswith('mergem_') else map_metabolite_to_mergem_id(product)
            stoc = reaction.metabolites[product] if exact_sto else 1
            metabolite_set = (id, stoc)
            rev_met_set = (id, -stoc)
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
    if len(models) > 1 and set_objective == 'merge':
        merged_model, reac_sources_dict = create_merged_objective(merged_model, reac_sources_dict, objective_reactions)
    else:
        model_index = int(1 if set_objective == 'merge' else set_objective) - 1
        reactions = objective_reactions[model_index]
        if len(reactions):
            merged_model.add_reactions(reactions)
            merged_model.objective = models[model_index].objective.expression
            for reaction in reactions:
                reac_sources_dict[reaction.id][model_index].append(reaction.id)

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

    for model_index, reaction_list in enumerate(objective_reactions):
        for reaction in reaction_list:
            reac_sources_dict[merged_obj_reaction_id][model_index].append(reaction.id)
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
