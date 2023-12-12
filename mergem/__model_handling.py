"""
    Functions to handle cobra models including transforming metabolite IDs to mergem namespace,
    creating reaction keys to compare them better, calculating jaccard distances, and
    loading and exporting cobra models.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""

from .__database_processing import build_id_mapping
import cobra
# This hack solves the problem of cobrapy replacements introducing control ASCII characters in ids,
# which breaks the glpk solver and crashes the Python kernel
# See _f_reaction in sbml.py in cobrapy: __(NUMBER)__ replaced with the character value of NUMBER
from typing import Match
def _number_to_chr_safe(numberStr: Match) -> str:
    ascii = int(numberStr.group(1))
    return numberStr.group(1) if (ascii <= 31 or ascii == 127) else chr(ascii)
cobra.io.sbml._number_to_chr = _number_to_chr_safe

from pickle import dump, load
import os
import csv

met_univ_id_dict, met_univ_id_prop_dict, reac_univ_id_dict, reac_univ_id_prop_dict = {}, {}, {}, {}

curr_dir = os.path.dirname(__file__)
data_dir = os.path.join(curr_dir, 'data')
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

met_univ_id_dict_file = os.path.join(data_dir, 'metaboliteIdMapper.p')
met_univ_id_prop_dict_file = os.path.join(data_dir, 'metaboliteInfo.p')
reac_univ_id_dict_file = os.path.join(data_dir, 'reactionIdMapper.p')
reac_univ_id_prop_dict_file = os.path.join(data_dir, 'reactionInfo.p')

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
        if not os.path.exists(file):
            print("Dictionary not found. Updating mapping dictionaries...")
            update_id_mapper()

        f = open(file, "rb")
        dict = load(f)
        f.close()

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


def update_id_mapper(delete_database_files = True):
    """
    Downloads the latest database files,
    merges the database identifiers based on common properties and saves the mapping tables as pickles.
    """
    global met_univ_id_dict, met_univ_id_prop_dict, reac_univ_id_dict, reac_univ_id_prop_dict

    met_univ_id_dict, met_univ_id_prop_dict, reac_univ_id_dict, reac_univ_id_prop_dict = \
        build_id_mapping(delete_database_files)

    with open(met_univ_id_dict_file, 'wb') as file:
        dump(met_univ_id_dict, file)

    with open(met_univ_id_prop_dict_file, 'wb') as file:
        dump(met_univ_id_prop_dict, file)

    with open(reac_univ_id_dict_file, 'wb') as file:
        dump(reac_univ_id_dict, file)

    with open(reac_univ_id_prop_dict_file, 'wb') as file:
        dump(reac_univ_id_prop_dict, file)


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

    met_id = met_id.replace('~', '')
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

    reac_id = reac_id.replace('~', '')
    reac_univ_id = reac_univ_id_dict.get(reac_id)

    if reac_univ_id is None:
        reac_univ_id = reac_univ_id_dict.get(remove_localization(reac_id))

    return reac_univ_id


def get_metabolite_properties(met_univ_id):
    """
    Retrieves the properties of a metabolite using its universal id
    """
    if not met_univ_id_prop_dict:
        load_met_univ_id_prop_dict()
    #mergem dictionary
    return met_univ_id_prop_dict.get(met_univ_id)


def get_reaction_properties(reac_univ_id):
    """
    Retrieves the properties of a reaction using its universal id
    """
    if not reac_univ_id_prop_dict:
        load_reac_univ_id_prop_dict()

    return reac_univ_id_prop_dict.get(reac_univ_id)


def remove_localization(id):
    if '@' in id:
        return id.rsplit('@', 1)[0]
    else:
        return id.rsplit("_", 1)[0]


def save_mapping_tables(metabolites_file_name = 'mergem_univ_id_mapper_metabolites.csv',
                        reactions_file_name = 'mergem_univ_id_mapper_reactions.csv'):
    """
    Saves database id mapping tables as CSV files
    """
    if not met_univ_id_prop_dict:
        load_met_univ_id_prop_dict()
    if not reac_univ_id_prop_dict:
        load_reac_univ_id_prop_dict()

    property_dicts = {metabolites_file_name: met_univ_id_prop_dict,
                      reactions_file_name: reac_univ_id_prop_dict}

    for filename, property_dict in property_dicts.items():
        list_univ_ext_ids = []
        for univ_id, prop in property_dict.items():
            list_ids = [univ_id]
            list_ids += sorted(prop['ids'])
            list_univ_ext_ids += [list_ids]

        with open(filename, 'w', newline='') as f:
            write = csv.writer(f)
            write.writerows(list_univ_ext_ids)


def merge_unique(a,b):
    '''
    Merge strings and list of strings without duplicates
    :param a: a string or list of strings
    :param b: a string or list of strings
    :return: a string or unique list of strings
    '''
    if not a:
        return b
    
    if type(a) != list: a = [a]
    if type(b) != list: b = [b]
    a += [x for x in b if x not in a]
    if len(a) == 1:
        return a[0]
    else:
        return a


def add_annotations(id, dict_annot, source):
    '''
    Add annotations from a metabolite or reaction without duplicates
    :param id: id of metabolite or reaction to add annotations
    :param dict_annot: dictionary with all annotations
    :param source: reaction or metabolite source
    '''
    annotations =  dict_annot[id]
    for source_anot_key, source_anot_value in source.annotation.items():
        if anot_value := annotations.get(source_anot_key):
            if anot_value != source_anot_value:
                annotations[source_anot_key] = merge_unique(anot_value, source_anot_value)
        else:
            annotations[source_anot_key] = source_anot_value


def extend_metabolite_annotations(metabolite, props):
    annotation = metabolite.annotation
    for db_ids in props.get('ids', []):
        db_name, db_id = db_ids.split(':', 1)

        if 'bigg' in db_name:
            annotation['bigg.metabolite'] = merge_unique(annotation.get('bigg.metabolite'), db_id)

        elif 'chebi' in db_name:
            annotation['chebi'] = merge_unique(annotation.get('chebi'), 'CHEBI:' + db_id)

        elif 'metanetx' in db_name:
            annotation['metanetx.chemical'] = merge_unique(annotation.get('metanetx.chemical'), db_id)

        elif 'seed' in db_name:
            annotation['seed.compound'] = merge_unique(annotation.get('seed.compound'), db_id)
        
        elif 'kegg' in db_name:
            annotation['kegg.compound'] = merge_unique(annotation.get('kegg.compound'), db_id)

        elif 'reactome' in db_name:
            annotation['reactome.compound'] = merge_unique(annotation.get('reactome.compound'), db_id)

        else:
            annotation[db_name] = merge_unique(annotation.get(db_name), db_id)

    if (not annotation.get('inchikey')) and (not annotation.get('inchi_key')) and props.get('inchikey'):
        annotation['inchikey'] = props['inchikey']

    if not annotation.get('sbo'):
        annotation['sbo'] = 'SBO:0000247'

    if not metabolite.formula and len(props.get('formula', [])) > 0:
        metabolite.formula = props['formula'][0]


def extend_reaction_annotations(reaction, props):
    annotation = reaction.annotation
    for db_ids in props.get('ids', []):
        db_name, db_id = db_ids.split(':', 1)

        if 'bigg' in db_name:
            annotation['bigg.reaction'] = merge_unique(annotation.get('bigg.reaction'), db_id)

        elif 'metanetx' in db_name:
            annotation['metanetx.reaction'] = merge_unique(annotation.get('metanetx.reaction'), db_id)

        elif 'seed' in db_name:
            annotation['seed.reaction'] = merge_unique(annotation.get('seed.reaction'), db_id)

        else:
            annotation[db_name] = merge_unique(annotation.get(db_name), db_id)

    if not annotation.get('ec-code') and props.get('EC_num'):
        annotation['ec-code'] = props['EC_num']


def add_gpr(id, dict_gprs, source):
    '''
    :param id: id of reaction for which gpr is being updated
    :param dict_gprs: dictionary with all reaction gprs
    :param source: reaction with gpr to add to dict_gprs
    '''
    if source_gpr_str := cobra.core.gene.GPR.to_string(source.gpr):
        if gpr := dict_gprs.get(id):
            if gpr_str := cobra.core.gene.GPR.to_string(gpr):
                if source_gpr_str not in gpr_str:
                    dict_gprs[id] = cobra.core.gene.GPR.from_string(gpr_str + ' or ' + source_gpr_str)
            else:
                dict_gprs[id] = source.gpr
        else:
            dict_gprs[id] = source.gpr