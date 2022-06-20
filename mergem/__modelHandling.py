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
    met_mergem_id = __mergem_met_id_dict.get(metabolite.id)

    if met_mergem_id is None:
        if '@' in metabolite.id:
            split = metabolite.id.rsplit('@', 1)
        else:
            split = metabolite.id.rsplit("_", 1)
        met_mergem_id =  __mergem_met_id_dict.get(split[0])


    if met_mergem_id is not None:
        met_compartment = map_localization(metabolite.compartment)
        if (met_compartment == ''):
            met_compartment = map_localization(split[1])

        met_mergem_id = "mergem_" + str(met_mergem_id) + "_" + met_compartment

    return met_mergem_id


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
        merged_model.add_reactions(list_of_obj_reacs[model_num])
        merged_model.objective = list_of_models[model_num].objective.expression
        return merged_model


# loads and returns cobra model based on file format
def load_model(filename):
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


# saves a cobra model in a file with appropriate format
def save_model(cobra_model, file_name):
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


# Initialize mapper
__load_or_create_id_mapper()