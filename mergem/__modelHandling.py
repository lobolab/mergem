"""
    Functions to handle cobra models including transforming metabolite IDs to mergem namespace,
    creating reaction keys to compare them better, calculating jaccard distances, and
    loading and exporting cobra models.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""

import cobra.io
from cobra import Model, Reaction
from pickle import load
import os
from . import __database_id_merger

curr_dir = os.path.dirname(__file__)
__mergem_met_id_dict, __mergem_met_info_dict = {}, {}

localization_dict = {'p': 'p', 'p0': 'p', 'periplasm': 'p', 'periplasm_0': 'p', 'mnxc19': 'p',
                     'c': 'c', 'c0': 'c', 'cytosol': 'c', 'cytosol_0': 'c', 'cytoplasm': 'c', 'mnxc3': 'c',
                     'cytoplasm_0': 'c',
                     'e': 'e', 'e0': 'e', 'extracellular': 'e', 'extracellular_0': 'e',
                     'extracellular space': 'e', 'mnxc2': 'e',
                     'm': 'm', 'mitochondria': 'm', 'mitochondria_0': 'm', 'mnxc4': 'm',
                     'b': 'b', 'boundary': 'b',
                     'x': 'p/glyoxysome', 'mnxc24': 'p/glyoxysome', 'mnxc13': 'p/glyoxysome',
                     'h': 'h', 'choloroplast': 'h', 'mnxc8': 'h',
                     'v': 'v', 'vacuole': 'v', 'mnxc9': 'v',
                     'n': 'n', 'nucleus': 'n', 'mnxc6': 'n'}

# To convert model met ID to mergem ID
def __load_or_create_id_mapper():
    """
    Checks if metabolite id mapper exists in current directory, else downloads the latest database files,
    merges the database identifiers based on common properties and saves the mapping table as a pickle.
    """
    global __mergem_met_id_dict, __mergem_met_info_dict
    met_conv_file = "metaboliteIdMapper.p"
    met_info_file = 'metaboliteInfo.p'
    met_conv_pickle_path = os.path.join(curr_dir, "data", met_conv_file)
    met_info_pickle_path = os.path.join(curr_dir, "data", met_info_file)

    if (os.path.exists(met_conv_pickle_path)) and (os.path.exists(met_info_pickle_path)):
        f = open(met_conv_pickle_path, "rb")
        __mergem_met_id_dict = load(f)
        f.close()

        f = open(met_info_pickle_path, "rb")
        __mergem_met_info_dict = load(f)
        f.close()

    else:
        print("Building ID mapping and metabolite Info table.")
        __database_id_merger.__create_id_mapping_pickle()
        __mergem_met_id_dict = __database_id_merger.__return_mapping_and_info_dicts()[2]
        __mergem_met_info_dict = __database_id_merger.__return_mapping_and_info_dicts()[3]


def __update_id_mapping_pickles():
    """
    Downloads latest versions of database files, merges the metabolite identifiers that have common
    properties and saves the mapping table as a pickle.
    """
    __database_id_merger.__create_id_mapping_pickle()


# returns mergem identifier for a metabolite
def __get_mergem_id(db_id):
    """
    Returns mergem_id from metabolite ID mapping dictionary.
    :param db_id: External database id from model
    :return: metabolite mergem_id, if exists else None
    """
    mergem_id = __mergem_met_id_dict.get(db_id)

    if mergem_id is None:
        if '@' in db_id:
            split_db_id = db_id.rsplit('@', 1)[0]
        else:
            split_db_id = db_id.rsplit("_", 1)[0]
        mergem_id = __mergem_met_id_dict.get(split_db_id)

    return mergem_id


# convert cellular localization to single namespace
def map_localization(id_or_model_localization):
    """
    Converts localization suffixes into common notation.
    :param id_or_model_localization: cellular localization of entity in model
    :return: single letter cellular localization
    """
    localization = localization_dict.get(id_or_model_localization.lower(), '')

    return localization


# reaction key is a frozenset of tuples of participating mets with their stoichiometric coeffs
def __create_reaction_key(reaction, reverse=False):
    """
    Takes a reaction object as input and creates a key(frozen set) of all pairs of metabolite ID and stoichiometric
    coefficients. \n
    :param reaction: Cobra reaction object
    :param reverse: reverse the reaction before creating key
    :return: frozen set of pairs of IDs of participating metabolite and their stoichiometric coefficients
    """
    reac_metabolite_set = set()
    reac_rev_met_set = set()
    for reactant in reaction.reactants:
        if ('mergem_78_' not in reactant.id) and (reactant.id[-1] != 'b') and (reactant.name != "PMF"):
            metabolite_set = (reactant.id, -1)
            rev_met_set = (reactant.id, 1)
            reac_metabolite_set.add(metabolite_set)
            reac_rev_met_set.add(rev_met_set)

    for product in reaction.products:
        if ('mergem_78_' not in product.id) and (product.id[-1] != 'b') and (product.name != "PMF"):
            metabolite_set = (product.id, 1)
            rev_met_set = (product.id, -1)
            reac_metabolite_set.add(metabolite_set)
            reac_rev_met_set.add(rev_met_set)

    reac_metabolite_set = frozenset(reac_metabolite_set)
    reac_rev_met_set = frozenset(reac_rev_met_set)

    return reac_metabolite_set, reac_rev_met_set


# returns a metabolite id in mergem namespace with cellular localization
def map_to_metabolite_mergem_id(metabolite):
    """
    Takes a metabolite object as input and returns mergem_id notation for metabolite
    :param metabolite: Cobra metabolite object
    :return: mergem_id notation for the metabolite
    """
    met_mergem_id = __get_mergem_id(metabolite.id)
    met_compartment = map_localization(metabolite.compartment)

    if met_mergem_id is not None:
        if (met_compartment == '') & ('_' in metabolite.id):
            if '@' in metabolite.id:
                met_compartment = map_localization(metabolite.id.rsplit('@', 1)[1])
            else:
                met_compartment = map_localization(metabolite.id.rsplit("_", 1)[1])

            if met_compartment != '':
                met_id = "mergem_" + str(met_mergem_id) + "_" + met_compartment
                return met_id

        elif met_compartment != '':
            met_id = "mergem_" + str(met_mergem_id) + "_" + met_compartment
            return met_id

    return None


# create a reaction that contains input model objective reacs merged together
def __create_merged_objective(input_merged_model, list_of_obj_reac_lists):
    """
    Merges objective reactions in list and sets as objective in input merged model.
    :param input_merged_model: Merged model to set objective of
    :param list_of_obj_reac_lists: list of objective reaction lists to merge
    :return: Merged model with its objective set to the merged objective reactions
    """
    merged_reaction = Reaction('merged_objectives')
    st_dict = {}
    metabolite_dict, obj_mets_dict = {}, {}

    for reaction_list in list_of_obj_reac_lists:
        for reaction in reaction_list:
            for metabolite in reaction.metabolites:
                if metabolite.id not in st_dict:
                    st_dict[metabolite.id] = set()

                if metabolite.id not in metabolite_dict:
                    metabolite_dict[metabolite.id] = metabolite.copy()

                st_dict[metabolite.id] |= {reaction.metabolites[metabolite]}

            merged_reaction.name += reaction.id + "; "

    for metabolite_id, met_stoichiometries in st_dict.items():
        avg_stoichiometry = sum(met_stoichiometries)/len(met_stoichiometries)
        merged_reaction.add_metabolites({metabolite_dict[metabolite_id]: avg_stoichiometry})

    if len(merged_reaction.metabolites) > 0:
        input_merged_model.add_reactions({merged_reaction})
        input_merged_model.objective = merged_reaction.id

    return input_merged_model


# Creates a model containing the merged reactions and metabolites in original namespaces
def __convert_template_to_merged_model(model, merged_model_id, dict_met_id_conv):
    """
    Restore model identifiers from model id conv dict
    :param model: merged model in mergem namespace
    :param merged_model_id: model id for merged model
    :param dict_met_id_conv: dictionary mapping mergem ids to original model ids
    :return: merged model in original model id namespace
    """
    global __mergem_met_info_dict
    merged_model = Model(merged_model_id)

    for reaction in model.reactions:
        reaction_copy = reaction.copy()
        for metabolite in reaction_copy.metabolites:
            if "mergem" in metabolite.id:  # Only revert mergem IDs
                new_met_id = dict_met_id_conv[metabolite.id][0]
                if new_met_id not in reaction_copy.metabolites:
                    metabolite.id  = new_met_id
                else:
                    uniq_ids = [model_id for model_id in dict_met_id_conv[metabolite.id]
                                if model.metabolites.get_by_id(model_id) not in reaction.metabolites]
                    if len(uniq_ids) > 0:
                        metabolite.id = uniq_ids[0]
        merged_model.add_reactions([reaction_copy])

    merged_model.objective = model.objective

    return merged_model


# sets the objective for merged model
def __set_objective_expression(merged_model, list_of_models, list_of_obj_reacs, set_objective):
    """
    Sets the objective expression for merged model.
    :param merged_model: Model whose objective is to be set.
    :param list_of_models: List of all input models from which objective is chosen.
    :param list_of_obj_reacs: List of (lists of) all objective reactions in input models.
    :param set_objective: 'merge' all input model objective reacs or from one of the input models.
    :return: model with its objective expression set.
    """
    if set_objective == 'merge':
        merged_model = __create_merged_objective(merged_model, list_of_obj_reacs)
        return merged_model

    else:
        model_num = int(set_objective) - 1
        for reaction in list_of_obj_reacs[model_num]:
            merged_model.add_reactions({reaction.copy()})
        merged_model.objective = list_of_models[model_num].objective.expression
        return merged_model


# converts source model sets into lists
def __set_to_list_source_dict(dict_source):
    """
    Converts set values in dictionary to lists.
    :param dict_source: Dictionary with set values
    :return: Dictionary with list values
    """
    new_dict_of_sources = {}
    for key, value in dict_source.items():
        new_dict_of_sources[key] = list(value)

    return new_dict_of_sources


# loads and returns cobra model based on file format
def __load_model(filename):
    """
    Loads a model from the given filename/path.
    :param filename: Name of file to load model from.
    :return: Cobra model loaded from file.
    """
    if not os.path.exists(filename):
        raise IOError('File {} not found.'.format(filename))

    file_format = os.path.splitext(filename)[1][1:].strip().lower()

    if file_format in ["sbml", "xml"]:
        cobra_model = cobra.io.read_sbml_model(filename)
    elif file_format in ["mat", "m", "matlab"]:
        cobra_model = cobra.io.load_matlab_model(filename)
    elif file_format == "json":
        cobra_model = cobra.io.load_json_model(filename)
    elif file_format == "yaml":
        cobra_model = cobra.io.load_yaml_model(filename)
    else:
        raise TypeError('Cannot load file of {} format'.format(file_format))

    return cobra_model


# saves the merged model in xml format
def __export_merged_model(cobra_model, file_name):
    """
    Takes a cobra model as input and exports as file_name.
    :param cobra_model: cobra model to be saved.
    :param file_name: filename with format extension to save model as.
    """
    file_format = os.path.splitext(file_name)[1][1:].strip().lower()
    if file_format in ["sbml", "xml"]:
        cobra.io.write_sbml_model(cobra_model, file_name)
    elif file_format in ["mat", "m", "matlab"]:
        cobra.io.save_matlab_model(cobra_model, file_name)
    elif file_format == "json":
        cobra.io.save_json_model(cobra_model, file_name)
    elif file_format == "yaml":
        cobra.io.save_yaml_model(cobra_model, file_name)
    else:
        raise IOError('Unable to save merged model. Check file format {}'.format(file_name))


def __create_jaccard_matrix(num_models, met_source_dict, reac_source_dict):
    """
    Creates Jaccard distances matrix using dictionaries with source of each metabolite
    and reaction in merged model.
    :param num_models: number of input models
    :param met_source_dict: dictionary with source of each metabolite
    :param reac_source_dict: dictionary with source of each reaction
    :return: Jaccard distance matrices for metabolites and reactions
    """
    met_uniq_num = [0] * num_models
    reac_uniq_num = [0] * num_models
    jacc_matrix = []
    total_mets_merged, total_reacs_merged = 0, 0

    for met_sources in met_source_dict.values():
        if len(met_sources) > 1:
            total_mets_merged += 1
        else:
            met_uniq_num[list(met_sources)[0]] += 1

    for reac_sources in reac_source_dict.values():
        if len(reac_sources) > 1:
            total_reacs_merged += 1
        else:
            reac_uniq_num[list(reac_sources)[0]] += 1

    for i in range(0, num_models):
        jd_row = []
        for j in range(0, num_models):
            num_met_merged = 0
            num_reac_merged = 0

            if i > j:
                for reac_source in reac_source_dict.values():
                    if (i in reac_source) and (j in reac_source):
                        num_reac_merged += 1

                jd_row += [__calculate_jaccard_distance(reac_uniq_num[i], num_reac_merged, reac_uniq_num[j])]

            elif i < j:
                for met_source in met_source_dict.values():
                    if (i in met_source) and (j in met_source):
                        num_met_merged += 1

                jd_row += [__calculate_jaccard_distance(met_uniq_num[i], num_met_merged, met_uniq_num[j])]

            else:
                jd_row += [0]

        jacc_matrix += [jd_row]

    return jacc_matrix


def __calculate_jaccard_distance(num_reference, num_merged, num_unique):
    """
    Calculates jaccard distance given the number of entities in reference and
    number of entities common between reference and other set.
    :param num_reference: number of entities (met/reacs) in reference model
    :param num_merged: number of entities (met/reacs) merged
    :return: jaccard distance of model from reference model.
    """
    jaccard_distance = round(float(1 - num_merged/(num_unique + num_reference + num_merged)), 2)
    return jaccard_distance


# adding genes to merged reactions
def __update_gene_rule(existing_reaction, new_reaction):
    """
    Updates the gene rule for existing reaction from new reaction.
    :param existing_reaction: cobra reaction object that exists in model
    :param new_reaction: cobra reaction model that is being merged
    :return: a gene rule that includes rule from new reaction
    """
    if (existing_reaction.gene_reaction_rule == '') and (new_reaction.gene_reaction_rule != ''):
        updated_gene_rule = new_reaction.gene_reaction_rule
    elif (existing_reaction.gene_reaction_rule != '') and (new_reaction.gene_reaction_rule == ''):
        updated_gene_rule = existing_reaction.gene_reaction_rule
    elif (existing_reaction.gene_reaction_rule != '') and (new_reaction.gene_reaction_rule != ''):
        updated_gene_rule = existing_reaction.gene_reaction_rule + ' or ' + new_reaction.gene_reaction_rule
    else:
        updated_gene_rule = ''

    return updated_gene_rule
