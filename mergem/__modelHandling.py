"""
    Functions to handle cobra models including transforming metabolite IDs to fluxer namespace,
    creating reaction keys to compare them better, calculating jaccard distances, and
    loading and exporting cobra models.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""

import cobra.io
from cobra import Model, Reaction
from pickle import load
import os
from . import __database_met_id_merger

curr_dir = os.path.dirname(os.path.abspath(__file__))
__fluxer_met_id_dict = {}


# To convert model met ID to fluxer ID
def __load_or_create_id_mapper():
    """
    Checks if metabolite id mapper exists in current directory, else downloads the latest database files,
    merges the database identifiers based on common properties and saves the mapping table as a pickle.
    """
    global __fluxer_met_id_dict
    met_conv_file = "metaboliteIdMapper.p"
    pickle_path = os.path.join(curr_dir, "data", met_conv_file)
    if os.path.exists(pickle_path):
        f = open(pickle_path, "rb")
        __fluxer_met_id_dict = load(f)
        f.close()
    else:
        print("Building ID mapping table.")
        __database_met_id_merger.__create_id_mapping_pickle()


def __update_id_mapping_pickles():
    """
    Downloads latest versions of database files, merges the metabolite identifiers that have common
    properties and saves the mapping table as a pickle.
    """
    __database_met_id_merger.__create_id_mapping_pickle()


# returns fluxer identifier for a metabolite
def __get_fluxer_id(db_id):
    """
    Returns fluxer_id from metabolite ID mapping dictionary.
    :param db_id: External database id from model
    :return: metabolite fluxer_id, if exists else None
    """
    split_db_id = db_id.rsplit("_", 1)[0]
    fluxer_id = __fluxer_met_id_dict.get(db_id)

    if fluxer_id is None:
        fluxer_id = __fluxer_met_id_dict.get(split_db_id)

    return fluxer_id


# convert cellular localization to single namespace
def __get_localization(id_or_model_localization):
    """
    Converts localization suffixes into common notation.
    :param id_or_model_localization: cellular localization of entity in model
    :return: single letter cellular localization
    """
    localization_dict = {'p': 'p', 'p0': 'p', 'periplasm': 'p', 'periplasm_0': 'p',
                         'c': 'c', 'c0': 'c', 'cytosol': 'c', 'cytosol_0': 'c', 'cytoplasm': 'c',
                         'cytoplasm_0': 'c',
                         'e': 'e', 'e0': 'e', 'extracellular': 'e', 'extracellular_0': 'e',
                         'extracellular space': 'e',
                         'm': 'm', 'mitochondria': 'm', 'mitochondria_0': 'm',
                         'b': 'b', 'boundary': 'b',
                         'x': 'p/glyoxysome',
                         'h': 'h', 'choloroplast': 'h',
                         'v': 'v', 'vacuole': 'v',
                         'n': 'n', 'nucleus': 'n'}

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
    reac_metabolite_set = set()  # sets are unordered -> easy to compare to one another than tuple
    for metabolite in reaction.metabolites:
        if reverse:
            stoichiometric_coeff = - reaction.metabolites[metabolite]
        else:
            stoichiometric_coeff = reaction.metabolites[metabolite]

        metabolite_set = (metabolite.id, stoichiometric_coeff)  # since models are already translated
        reac_metabolite_set.add(metabolite_set)

    reac_metabolite_set = frozenset(reac_metabolite_set)

    return reac_metabolite_set


# returns a metabolite id in fluxer namespace with cellular localization
def __create_fluxer_metabolite_id(metabolite):
    """
    Takes a metabolite object as input and returns fluxer_id notation for metabolite
    :param metabolite: Cobra metabolite object
    :return: Fluxer_id notation for the metabolite
    """
    met_fluxer_id = __get_fluxer_id(metabolite.id)
    met_compartment = __get_localization(metabolite.compartment)

    if met_fluxer_id is not None:
        if (met_compartment == '') & ('_' in metabolite.id):
            met_compartment = __get_localization(metabolite.id.rsplit("_", 1)[1])
            if met_compartment != '':
                met_id = "fluxer_" + str(met_fluxer_id) + "_" + met_compartment
                return met_id

        elif met_compartment != '':
            met_id = "fluxer_" + str(met_fluxer_id) + "_" + met_compartment
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

    for reaction_list in list_of_obj_reac_lists:
        for reaction in reaction_list:
            merged_reaction = __merge_obj_reactions(merged_reaction, reaction)
            merged_reaction.name += reaction.id + "; "

    if len(merged_reaction.metabolites) > 1:
        input_merged_model.add_reactions({merged_reaction})
        input_merged_model.objective = merged_reaction.id

    return input_merged_model


# merges the right reaction into merged reaction
def __merge_obj_reactions(merged_r, right_r, operation='avg'):
    """
    Merge right reaction into the first reaction.
    :param merged_r: Reaction to which other reaction is merged
    :param right_r: The reaction to be merged
    :param operation: 'Avg' the stoichiometries of common metabolites
    :return: merged reaction
    """
    merged_metabolites = {}

    for met in merged_r.metabolites:
        merged_metabolites[met.id] = merged_r.metabolites[met]  # store st coefficients

    for met in right_r.metabolites:
        right_met_stoichiometry = right_r.metabolites[met]
        merged_met_stoichiometry = merged_metabolites.get(met.id, 0)
        copy_met = met.copy()

        if (merged_met_stoichiometry != 0) & (merged_met_stoichiometry != right_met_stoichiometry):
            if operation == 'avg':
                met_stoichiometry = ((merged_met_stoichiometry + right_met_stoichiometry) / 2.0)
            elif operation == 'max':
                met_stoichiometry = max(merged_met_stoichiometry, right_met_stoichiometry)
            else:
                met_stoichiometry = min(merged_met_stoichiometry, right_met_stoichiometry)

        elif (merged_met_stoichiometry != 0) & (merged_met_stoichiometry == right_met_stoichiometry):
            continue

        elif merged_met_stoichiometry == 0:
            met_stoichiometry = right_met_stoichiometry

        merged_r.add_metabolites({copy_met: met_stoichiometry})

    return merged_r


# Creates a model containing the merged reactions and metabolites in original namespaces
def __convert_template_to_merged_model(model, merged_model_id, dict_met_id_conv):
    """
    Restore model identifiers from model id conv dict
    :param model: merged model in fluxer namespace
    :param merged_model_id: model id for merged model
    :param dict_met_id_conv: dictionary mapping fluxer ids to original model ids
    :return: merged model in original model id namespace
    """
    merged_model = Model(merged_model_id)

    for reaction in model.reactions:
        reaction_copy = reaction.copy()
        for metabolite in reaction_copy.metabolites:
            if "fluxer" in metabolite.id:  # Only revert fluxer IDs
                metabolite.id = dict_met_id_conv[metabolite.id][0]  # converted_mets[metabolite.id]
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


def __calculate_jaccard_distance(num_reference, num_merged):
    """
    Calculates jaccard distance given the number of entities in reference and
    number of entities common between reference and other set.
    :param num_reference: number of entities (met/reacs) in reference model
    :param num_merged: number of entities (met/reacs) merged
    :return: jaccard distance of model from reference model.
    """
    jaccard_distance = round(float(1 - num_merged/num_reference), 2)
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
