"""
    Functions to handle cobra models including transforming metabolite IDs to mergem namespace,
    creating reaction keys to compare them better, calculating jaccard distances, and
    loading and exporting cobra models.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""

import cobra.io
from pickle import dump, load
import os
from .__database_id_merger import build_id_mapping

met_univ_id_dict, met_univ_id_prop_dict, reac_univ_id_dict, reac_univ_id_prop_dict = {}, {}, {}, {}

curr_dir = os.path.dirname(__file__)
met_univ_id_dict_file = os.path.join(curr_dir, 'data', 'metaboliteIdMapper.p')
met_univ_id_prop_dict_file = os.path.join(curr_dir, 'data', 'metaboliteInfo.p')
reac_univ_id_dict_file = os.path.join(curr_dir, 'data', 'reactionIdMapper.p')
reac_univ_id_prop_dict_file = os.path.join(curr_dir, 'data', 'reactionInfo.p')

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

proton_mergem_id = ''


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


def load_dict(dict, file):
    if not dict:
        if os.path.exists(file):
            f = open(file, "rb")
            dict = load(f)
            f.close()
        else:
            update_id_mapper()

    return dict


def load_met_univ_id_dict():
    global met_univ_id_dict, proton_mergem_id
    met_univ_id_dict = load_dict(met_univ_id_dict, met_univ_id_dict_file)
    proton_mergem_id = 'mergem_' + str(met_univ_id_dict['C00080']) + '_'


def load_met_univ_id_prop_dict():
    global met_univ_id_prop_dict
    met_univ_id_prop_dict = load_dict(met_univ_id_prop_dict, met_univ_id_prop_dict_file)


def load_reac_univ_id_dict():
    global reac_univ_id_dict
    reac_univ_id_dict = load_dict(reac_univ_id_dict, reac_univ_id_dict_file)


def load_reac_univ_id_prop_dict():
    global reac_univ_id_prop_dict
    reac_univ_id_prop_dict = load_dict(reac_univ_id_prop_dict, reac_univ_id_prop_dict_file)


def update_id_mapper():
    """
    Downloads the latest database files,
    merges the database identifiers based on common properties and saves the mapping table as a pickle.
    """
    global met_univ_id_dict, met_univ_id_prop_dict, reac_univ_id_dict, reac_univ_id_prop_dict

    print("Building ID mapping and properties information.")
    met_univ_id_dict, met_univ_id_prop_dict, reac_univ_id_dict, reac_univ_id_prop_dict = build_id_mapping()

    log("Creating reaction id pickle")
    with open(pickle_dir + 'reactionIdMapper.p', 'wb') as file:
        dump(dict_any_reac_id_to_univ_id, file)

    log("Creating reaction info pickle")
    with open(pickle_dir + 'reactionInfo.p', 'wb') as file:
        dump(dict_univ_id_to_reac_prop, file)

    log("Creating metabolite id pickle")
    with open(pickle_dir + 'metaboliteIdMapper.p', 'wb') as file:
        dump(dict_any_met_id_to_univ_id, file)

    log("Creating metabolite info pickle")
    with open(pickle_dir + 'metaboliteInfo.p', 'wb') as file:
        dump(dict_univ_id_to_met_prop, file)


# convert cellular localization to single namespace
def map_localization(id_or_model_localization):
    """
    Converts localization suffixes into common notation.
    :param id_or_model_localization: cellular localization of entity in model
    :return: single letter cellular localization
    """
    localization = localization_dict.get(id_or_model_localization.lower(), '')

    return localization


def map_metabolite_univ_id(met_id):
    """
    Maps metabolite id to metabolite universal id
    """
    if not met_univ_id_dict:
        load_met_univ_id_dict()

    met_univ_id = met_univ_id_dict.get(met_id)

    if met_univ_id is None:
        met_univ_id = met_univ_id_dict.get(remove_localization(met_id))

    return met_univ_id


def map_reaction_univ_id(reac_id):
    """
    Maps reaction id to metabolite universal id
    """
    if not reac_univ_id_dict:
        load_reac_univ_id_dict()

    reac_univ_id = reac_univ_id_dict.get(reac_id)

    if reac_univ_id is None:
        reac_univ_id = reac_univ_id_dict.get(remove_localization(reac_id))

    return reac_univ_id


def get_metabolite_properties(met_univ_id):
    if not met_univ_id_prop_dict:
        load_met_univ_id_prop_dict()

    return met_univ_id_prop_dict.get(met_univ_id)


def get_reaction_properties(reac_univ_id):
    if not reac_univ_id_prop_dict:
        load_reac_univ_id_prop_dict()

    return reac_univ_id_prop_dict.get(reac_univ_id)


def remove_localization(id):
    if '@' in id:
        return id.rsplit('@', 1)[0]
    else:
        return id.rsplit("_", 1)[0]
