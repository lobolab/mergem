"""
    Downloads metabolite information files from Kegg, Chebi, ModelSeed, MetaNetX and Bigg
    databases. Creates a dictionary mapping external database metabolite IDs
    to a Fluxer ID and another dictionary mapping Fluxer ID to metabolite properties.
    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""

import urllib.request
import shutil
from contextlib import closing
from datetime import datetime
from gzip import open as gzip_open
from sys import maxsize
from time import perf_counter
from pickle import dump, load
import os
import requests
import ssl

curr_dir = os.path.dirname(os.path.abspath(__file__))
files_dir = os.path.join(curr_dir, "downloads/")
pickle_dir = os.path.join(curr_dir, "/data/")
log_dir = os.path.join(curr_dir, "/logs/")

# Database URLs
modelSeed_met_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv"
modelSeed_met_aliases_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt"
metaNetX_compounds_url = "https://ftp.vital-it.ch/databases/metanetx/MNXref/latest/chem_prop.tsv"
metaNetX_met_xref_url = "https://ftp.vital-it.ch/databases/metanetx/MNXref/latest/chem_xref.tsv"
bigg_metabolites_url = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
chebi_compound_structure_url = "http://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/structures.csv.gz"
chebi_compounds_url = "http://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz"

# Metabolite filenames
ms_met_filename = files_dir + "modelSeed_compounds.tsv"
ms_met_aliases_filename = files_dir + "modelseed_compound_aliases.txt"
mx_chem_prop_filename = files_dir + "metanetx_chem_prop.tsv"
mx_met_xref_filename = files_dir + "metanetx_chem_xref.tsv"
bigg_metabolites_filename = files_dir + "bigg_models_metabolites.txt"
kegg_metabolites_filename = files_dir
chebi_compound_st_zipped_filename = files_dir + "chebi_structures.csv.gz"
chebi_compound_structure_filename = files_dir + "chebi_structures.csv"
chebi_compounds_zipped_filename = files_dir + "chebi_compounds.tsv.gz"
chebi_compounds_filename = files_dir + "chebi_compounds.tsv"

url_dictionary = {ms_met_filename: modelSeed_met_url,
                  ms_met_aliases_filename: modelSeed_met_aliases_url,
                  mx_chem_prop_filename: metaNetX_compounds_url,
                  mx_met_xref_filename: metaNetX_met_xref_url,
                  bigg_metabolites_filename: bigg_metabolites_url,
                  chebi_compounds_zipped_filename: chebi_compounds_url,
                  chebi_compound_st_zipped_filename: chebi_compound_structure_url
                  }

url_dictionary_chebi = {chebi_compounds_filename: chebi_compounds_zipped_filename,
                        chebi_compound_structure_filename: chebi_compound_st_zipped_filename
                        }

# Dictionaries
dict_any_met_id_to_fluxer_id = {}
dict_fluxer_id_to_met_prop = {}
list_primary_ids = set()
last_fluxer_id = 0
start_time = datetime.now().strftime("%Y%m%d_%HH%MM")
database_names_dict = {'seed': set(), 'metanetx': set(), 'bigg': set(), 'kegg': set()}
primary_dbs = ['seed', 'metanetx', 'bigg', 'kegg', 'chebi']
met_prop_collected = {'seed': set(), 'metanetx': set(), 'bigg': set(), 'kegg': set(), 'chebi': set()}


def __log(message):
    dt_string = datetime.now().strftime("%Y%m%d_%H:%M:%S")
    log_line = dt_string + " " + message
    print(log_line)
    log_file = open(log_dir + start_time + "_Metabolites_DatabasesDownloadLog.txt", "a")
    log_file.write(log_line + "\n")
    log_file.close()


def __create_directories():
    if not os.path.exists(files_dir):
        os.makedirs(files_dir)

    if not os.path.exists(pickle_dir):
        os.makedirs(pickle_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)


def __download_database_files():
    global kegg_metabolites_filename
    if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
            getattr(ssl, '_create_unverified_context', None)):
        ssl._create_default_https_context = ssl._create_unverified_context

    for filename in url_dictionary.keys():
        __log("Downloading " + filename)
        url = url_dictionary[filename]
        with closing(urllib.request.urlopen(url)) as r:
            with open(filename, 'wb') as f:
                shutil.copyfileobj(r, f)

        __log(filename + " downloaded.")

    for filename, zipped_filename in url_dictionary_chebi.items():
        with gzip_open(zipped_filename, 'rb') as file_in:
            with open(filename, 'wb') as file_out:
                shutil.copyfileobj(file_in, file_out)

    kegg_metabolites_filename += __check_for_kegg_update()


def __check_for_kegg_update():
    """
    Checks for a KEGG update and creates a filename for latest version.\n
    :return: filename for latest version
    """
    __log("Checking Kegg stats. ")
    kegg_cpd_stats_response = requests.get("http://rest.kegg.jp/info/cpd")
    release_version = kegg_cpd_stats_response.text.split('\n')[1].split('Release ', 1)[1]
    filename = "kegg_" + \
               release_version.translate({ord(c): "-" for c in "/!@#$%^&*()[]{};:,.<>?\\|`~-=+"}).replace(" ", "-")

    filename = filename.split("--", 1)[0] + ".p"

    if not os.path.isfile(files_dir + filename):
        __log("New Kegg version found: {}".format(release_version))
        __log("Downloading new Kegg version.")
        __download_kegg_compounds(filename)

    return filename


def __download_kegg_compounds(filename):
    """
    Uses API to get information on each compound from KEGG database.\n
    :param filename: filename to save kegg compound information as
    """
    kegg_compounds_request = requests.get("http://rest.kegg.jp/list/compound")
    kegg_compounds_list = kegg_compounds_request.text.split('\n')
    kegg_compounds_dict = {}

    for compound in kegg_compounds_list:
        if "\t" in compound:
            split_text = compound.split('\t', 1)
            kegg_id = (split_text[0]).split(':', 1)[1]
            kegg_compounds_dict[kegg_id] = {'Name': [], 'mass': [], 'formula': [], 'chebi': []}

            compound_info_request = requests.get("http://rest.kegg.jp/get/" + kegg_id)
            compound_info = compound_info_request.text.split('\n')

            for info_line in compound_info:
                line = info_line.replace(" ", "")
                if 'FORMULA' in line:
                    kegg_compounds_dict[kegg_id]['formula'].append(line.split('FORMULA', 1)[1])
                if 'NAME' in line:
                    kegg_compounds_dict[kegg_id]['Name'].append(line.split('NAME', 1)[1][:-1])
                elif 'MOL_WEIGHT' in line:
                    kegg_compounds_dict[kegg_id]['mass'].append(line.split('MOL_WEIGHT', 1)[1])
                elif 'EXACT_MASS' in line:
                    kegg_compounds_dict[kegg_id]['mass'].append(line.split('EXACT_MASS', 1)[1])
                elif 'chebi:' in line.lower():
                    kegg_compounds_dict[kegg_id]['chebi'].append(line.lower())
    __log(f"Kegg ids processed: {len(kegg_compounds_dict)}")

    with open(files_dir + filename, 'wb') as kegg_file:
        dump(kegg_compounds_dict, kegg_file)


def __process_met_file(file_name, file_line_reader):
    """
    Uses reader function to read lines of file and append informaiton to met properties dictionary.\n
    :param file_name: name of database file
    :param file_line_reader: reader function for database file
    """
    __log("Processing file " + file_name)
    if ".csv" in file_name:
        with open(file_name, "r") as db_file:
            next(db_file)  # skip header
            for line in db_file:
                met_properties = file_line_reader(line.strip().split(','))
                if met_properties is not None:
                    __append_met_properties(met_properties)
    else:
        with open(file_name, "r") as db_file:
            next(db_file)  # skip header
            for line in db_file:
                met_properties = file_line_reader(line.strip().split('\t'))
                if met_properties is not None:
                    __append_met_properties(met_properties)

    __log("Done processing file " + file_name)
    __log(f"Number of metabolite ids: {len(dict_any_met_id_to_fluxer_id)}")
    __log(f"Number of fluxer ids: {len(dict_fluxer_id_to_met_prop)}")
    __log(f"List of primary ids: {len(list_primary_ids)}")
    __log("")


def __process_cross_ref_info(file_name, xref_line_reader):
    """
    Read cross-reference information from file and maps database identifiers
    :param file_name: name of database file containing cross-reference information
    :param xref_line_reader: reader function for database file
    """
    __log("Processing file " + file_name)

    xref_dict = xref_line_reader(file_name)

    for source_id, xref_set in xref_dict.items():
        xref_list = list(xref_set)
        source_fluxer_id = dict_any_met_id_to_fluxer_id[source_id]
        source_id_prop = dict_fluxer_id_to_met_prop[source_fluxer_id]

        for other_id in xref_list:
            if other_id in list_primary_ids:
                __merge_identifiers(source_id, other_id)
            elif (other_id not in dict_any_met_id_to_fluxer_id) and (other_id.split(":", 1)[0] not in primary_dbs):
                dict_any_met_id_to_fluxer_id[other_id] = source_fluxer_id
                source_id_prop['ids'] += [other_id]

    __log("Done processing file " + file_name)
    __log(f"Number of metabolite ids: {len(dict_any_met_id_to_fluxer_id)}")
    __log(f"Number of fluxer ids: {len(dict_fluxer_id_to_met_prop)}")
    __log("")


def __append_met_properties(met_properties):
    global list_primary_ids
    fluxer_id = dict_any_met_id_to_fluxer_id.get(met_properties['ids'][0], maxsize)

    if fluxer_id == maxsize:
        global last_fluxer_id
        last_fluxer_id = last_fluxer_id + 1
        fluxer_id = last_fluxer_id

        dict_fluxer_id_to_met_prop[fluxer_id] = {'Name': [], 'ids': [], 'formula': [],
                                                 'mass': [], 'inchikey': [], 'xref_links': []}

    dict_any_met_id_to_fluxer_id[met_properties['ids'][0]] = fluxer_id
    fluxer_met_properties = dict_fluxer_id_to_met_prop[fluxer_id]
    list_primary_ids |= {db_id for db_id in met_properties['ids']}

    for key, value in met_properties.items():
        if (key == 'ids') and (len(value) > 1):
            for met_id in value:
                dict_any_met_id_to_fluxer_id[met_id] = fluxer_id
        __add_values_to_property_list(fluxer_met_properties[key], value)


def __process_kegg_compounds(filename):
    __log(f"Processing file {filename}")
    kegg_file = open(filename, "rb")
    kegg_compounds_dictionary = load(kegg_file)
    kegg_file.close()

    for kegg_id, properties in kegg_compounds_dictionary.items():
        ids = ["kegg:" + kegg_id]  # .lower()]

        if len(properties['chebi']) > 0:
            ids += [properties['chebi'][0]]

        property_dict = {'ids': ids,
                         'Name': properties['Name'],
                         'formula': properties['formula'],
                         'mass': properties['mass'],
                         }
        __append_met_properties(property_dict)

    met_prop_collected['kegg'] |= {'Name', 'formula', 'mass'}
    database_names_dict['kegg'] |= {'chebi', 'kegg'}
    __log(f"Done processing file {filename}")
    __log(f"Number of metabolite ids: {len(dict_any_met_id_to_fluxer_id)}")
    __log(f"Number of fluxer ids: {len(dict_fluxer_id_to_met_prop)}")
    __log("")


def __chebi_compounds_inchi_reader(line):
    """
    Reader function for chebi compounds. \n
    :param line: line from file
    :return: dictionary with properties from line
    """
    if 'inchikey' in (item.lower() for item in line):
        return {'ids': ['chebi:' + str(line[1])],
                'inchikey': [line[2]]
                }

    met_prop_collected['chebi'] |= {'inchikey'}


def __chebi_compounds_names(line):
    """
    Reader function for chebi compounds. \n
    :param line: line from file
    :return: dictionary with properties from line
    """
    if line[5] != 'null':
        return {'ids': [line[2].lower()],
                'Name': [line[5]]
                }

    met_prop_collected['chebi'] |= {'Name'}


def __metanetx_chem_prop_line_reader(line):
    """
    Reader function for metanetx compounds. \n
    :param line: line from file
    :return: dictionary with properties from line
    """
    if line[0][0] != "M":
        return None

    ids = ["metanetx:" + line[0]]
    met_prop_collected['metanetx'] |= {'Name', 'formula', 'mass', 'inchikey'}

    if len(line) > 7:
        if len(line[5]) > 0:
            mass = [line[5]]
        else:
            mass = []

        if len(line[7]) > 0:
            inchikey = [line[7].split('=', 1)[1]]  # .lower().split('=', 1)[1]]
        else:
            inchikey = []

        return {'ids': ids,
                'Name': [line[1]],
                'formula': [line[3]],
                'mass': mass,
                'inchikey': inchikey
                }
    else:
        return {'ids': ids,
                'Name': [line[1]],
                'mass': [0],
                'inchikey': []}


def __metanetx_chem_prop_xref_reader(file_name):
    """
    Reader function for cross references from metanetx compounds file. \n
    :param file_name: line from file
    :return: dictionary with properties from line
    """
    xref_dict = {}
    with open(file_name, "r") as db_file:
        next(db_file)  # skip header
        for line in db_file:
            if (line[0][0] != "M") or (":" not in line[2]):
                continue
            [other_db, other_db_id] = line[2].split(":", 1)
            source_id = "metanetx:" + line[0]
            other_id = None

            if other_db.lower() == "mnx":
                other_id = "metanetx:" + other_db_id
                other_db_name = "metanetx"

            elif other_db.lower() in {"kegge", "envipathm", "envipath", "biggm", "seedm", "reactomem", "lipidmapsm"}:
                continue

            elif other_db[:-1].lower() in {'bigg', 'kegg', 'seed', 'metacyc'}:
                other_db_name = other_db[:-1].lower()
                other_id = other_db_name + ":" + other_db_id

            elif other_db == "slm":
                other_db_name = "slm"
                other_id = other_db_name + ":" + other_db_id

            if other_id is not None:
                database_names_dict['metanetx'] |= {other_db_name}
                if source_id in xref_dict:
                    xref_dict[source_id] |= {other_id}
                else:
                    xref_dict[source_id] = {other_id}
    database_names_dict['metanetx'] |= {'metanetx'}
    return xref_dict


def __metanetx_chem_xref_reader(file_name):
    """
        Reader function for metanetx cross reference file. \n
        :param file_name: name of database file
        :return: dictionary with properties from line
    """
    xref_dict = {}
    with open(file_name, "r") as db_file:
        next(db_file)  # skip header
        for line in db_file:
            line = line.strip().split('\t')

            if (line[0][0] == "#") or (line[0][0:3] == "MNX") or (line[0][0:3] == "mnx") or (line[1] == "MNXM0") or \
                    (line[1] == "BIOMASS") or (line[1] == "WATER") or (":" not in line[0]):
                continue

            if ("unknown" in line[2].lower()) or ("no description" in line[2].lower()) or ("obsolete" in line[2].lower()) \
                    or ("molecular entity" in line[2].lower()):
                continue

            source_id = "metanetx:" + line[1]  # .lower()
            [other_db, other_db_id] = line[0].split(":", 1)  # .lower().split(":", 1)
            other_id = None

            if "." in other_db:
                other_db = other_db.split(".")[0]
                other_id = other_db + ":" + other_db_id

            elif other_db not in {"keggC", "envipathM", "envipath", "seedM", "CHEBI", "biggM", "keggD", "SLM", "keggE",
                                  "keggG", "reactomeM", "sabiorkM", "rheaP", "rheaG", "lipidmapsM", "metacycM"}:
                other_id = other_db + ":" + other_db_id

            if other_id is not None:  # and (other_db_name in primary_dbs):
                database_names_dict['metanetx'] |= {other_db}
                if source_id in xref_dict:
                    xref_dict[source_id] |= {other_id}
                else:
                    xref_dict[source_id] = {other_id}

    return xref_dict


def __modelseed_metabolites_line_reader(line):
    """
    Line reader function for modelseed metabolites. \n
    :param line: line from file
    :return: dictionary of metabolite information from line
    """
    ids = ["seed:" + line[0]]  # .lower()]
    met_prop_collected['seed'] |= {'Name', 'formula', 'mass', 'inchikey'}

    if (len(line[4]) > 0) and (line[4].lower() != 'none') and (line[4].lower() != 'null'):
        mass = [line[4]]
    else:
        mass = []

    return {'ids': ids,
            'Name': [line[2]],
            'formula': [line[3]],
            'mass': mass,
            'inchikey': [line[6]]
            }


def __modelseed_metabolites_xref_reader(file_name):
    """
    Reader function for cross reference information from modelseed metabolites file. \n
    :param file_name: modelseed metabolites file name
    :return: dictionary with cross reference information
    """
    xref_dict = {}
    with open(file_name, "r") as db_file:
        next(db_file)
        for line in db_file:
            line = line.strip().split('\t')
            source_id = "seed:" + line[0]
            ids = []
            if line[10] != 'null':
                ids += [('seed:' + ele) for ele in (line[10].split(";"))]  # .lower().split(";"))]

            if len(ids) > 0:
                if source_id in xref_dict:
                    xref_dict[source_id] |= {seed_id for seed_id in ids}
                else:
                    xref_dict[source_id] = {seed_id for seed_id in ids}
    database_names_dict['seed'] |= {'seed'}
    return xref_dict


def __modelseed_met_aliases_reader(file_name):
    """
    Reader function for modelseed metabolite aliases file.\n
    :param file_name: name of metabolite aliases file from modelseed
    :return: dictionary of cross reference information
    """

    xref_dict = {}
    with open(file_name, "r") as db_file:
        next(db_file)  # skip header
        for line in db_file:
            line = line.strip().split('\t')
            source_id = "seed:" + line[0].lower()
            other_db = line[2].lower()

            if other_db == 'metanetx.chemical':
                other_db_name = 'metanetx'
                database_names_dict['seed'] |= {other_db_name}
                other_id = other_db_name + ":" + line[1]  # .lower()

            elif other_db in primary_dbs:
                other_db_name = other_db
                database_names_dict['seed'] |= {other_db_name}
                other_id = other_db_name + ":" + line[1]  #.lower()

            else:
                continue

            if source_id in xref_dict:
                xref_dict[source_id] |= {other_id}
            else:
                xref_dict[source_id] = {other_id}

    return xref_dict


def __bigg_metabolites_line_reader(line):
    """
    Line reader for metabolites from BiGG database.\n
    :param line: line in file
    :return: dictionary with metabolite information from line
    """
    if "recon" in line[1].lower():
        return None

    ids = []
    inchikeys = []
    xref_links = []
    names = []
    met_prop_collected['bigg'] |= {'Name', 'inchikey', 'xref_links'}

    ids += ['bigg:' + line[1]]

    if len(line) > 2:
        names += [line[2]]

    database_names_dict['bigg'] |= {'bigg'}

    if (len(line) > 4) and "http" in line[4]:
        for link in line[4].split('; '):
            link_part = link.split('/')
            other_db = link_part[3].lower()
            other_db_id = link_part[4]

            if other_db == 'inchikey':
                inchikeys += [other_db_id]

            if "kegg" in link:
                xref_links += ["http" + link.split("http", 1)[1]]
            elif ("metanetx" in link) or ("seed" in link):
                xref_links += ["http" + link.split("http", 1)[1]]

    prop = {}
    if ids:
        prop['ids'] = ids
    if inchikeys:
        prop['inchikey'] = inchikeys
    if xref_links:
        prop['xref_links'] = xref_links
    if names:
        prop['Name'] = names

    if prop:
        return prop
    else:
        return None


def __bigg_models_xref_reader(file_name):
    """
    Reader function for cross reference information from BiGG metabolites file.
    :param file_name: name of file
    :return: dictionary of cross-reference information
    """
    xref_dict = {}
    with open(file_name, "r") as db_file:
        next(db_file)  # skip header
        for line in db_file:
            line = line.strip().split('\t')

            if len(line) < 5:
                continue

            ids = []
            source_id = 'bigg:' + line[1]

            if "http" in line[4]:
                for link in line[4].split('; '):
                    link_part = link.split('/')
                    other_db = link_part[3].lower()
                    other_db_id = link_part[4]

                    if other_db == 'inchikey':
                        continue
                    elif other_db == 'chebi':
                        database_names_dict['bigg'] |= {'chebi'}
                        ids += [other_db_id.lower()]
                    else:
                        other_db_name = other_db.split('.')[0].lower()
                        database_names_dict['bigg'] |= {other_db_name}
                        if other_db_name in primary_dbs:  # get only primary identifiers from bigg - incorrect chebi
                            ids += [other_db_name + ':' + other_db_id]

            if source_id in xref_dict:
                xref_dict[source_id] |= {db_id for db_id in ids}
            else:
                xref_dict[source_id] = {db_id for db_id in ids}
    return xref_dict


def __add_values_to_property_list(property_list, prop_value):
    invalid_values_list = ['\'\'', '\"\"', 'null', '-', '']
    for value in prop_value:
        if (value not in invalid_values_list) and (value not in property_list):
            property_list.append(value)


def __merge_identifiers(source_id, other_id):
    """
    Merges the two met IDs into lowest fluxer ID only if at least two properties match.
    :param source_id: primary metabolite ID from database being processed
    :param other_id: cross referenced metabolite ID to be mapped to primary met ID
    """
    source_fluxer_id = dict_any_met_id_to_fluxer_id.get(source_id)
    source_id_properties = dict_fluxer_id_to_met_prop[source_fluxer_id]

    matched_properties = []
    check_prop = []
    other_fluxer_id = dict_any_met_id_to_fluxer_id.get(other_id)

    if other_fluxer_id != source_fluxer_id:
        other_properties = dict_fluxer_id_to_met_prop[other_fluxer_id]

        check_prop += [prop for prop in ['mass', 'inchikey', 'Name', 'formula']
                       if len(other_properties[prop]) > 0 if len(source_id_properties[prop]) > 0]

        property_mismatch = False
        for prop in check_prop:
            if prop == 'mass':
                if abs(float(other_properties[prop][0]) - float(source_id_properties[prop][0])) < 0.5:
                    matched_properties.append(prop)
                else:
                    property_mismatch = True
            else:
                if prop != 'inchikey':
                    for string_prop in other_properties[prop]:
                        if string_prop.lower() in \
                                (source_string_prop.lower() for source_string_prop in source_id_properties[prop]):
                            matched_properties.append(prop)
                            break
                else:
                    if other_properties[prop][0].lower() == source_id_properties[prop][0].lower():
                        matched_properties.append(prop)
                    else:
                        property_mismatch = True
            if (len(matched_properties) > 1) or property_mismatch:
                break

        if len(matched_properties) > 1:
            minimum_fluxer_id = min([source_fluxer_id, other_fluxer_id])

            if minimum_fluxer_id == source_fluxer_id:
                for xref_id in other_properties['ids']:
                    dict_any_met_id_to_fluxer_id[xref_id] = minimum_fluxer_id
                for key, value in other_properties.items():
                    __add_values_to_property_list(source_id_properties[key], value)
                del other_properties

            else:
                for xref_id in source_id_properties['ids']:
                    dict_any_met_id_to_fluxer_id[xref_id] = minimum_fluxer_id
                for key, value in source_id_properties.items():
                    __add_values_to_property_list(other_properties[key], value)
                del source_id_properties


def __clean_id_mapping_dictionary(id_dictionary, info_dictionary):
    dict_copy_id_converter = {}
    dict_copy_met_info = {}
    db_name_dict = {}
    db_preference = {'kegg': 0, 'chebi': 1, 'metanetx': 2, 'bigg': 3, 'seed': 4, 'sabiork': 5}
    for db_id, fl_id in id_dictionary.items():
        [db_name, new_key] = db_id.split(":", 1)
        conflict_fl_id = dict_copy_id_converter.get(new_key)
        if conflict_fl_id:
            if db_preference.get(db_name, maxsize) > db_preference.get(db_name_dict[conflict_fl_id], maxsize):
                fl_id = conflict_fl_id
            else:
                continue
        dict_copy_id_converter[new_key] = fl_id
        if db_name_dict.get(fl_id) is None:
            db_name_dict[fl_id] = db_name

    for fluxer_id in dict_copy_id_converter.values():
        dict_copy_met_info[fluxer_id] = info_dictionary[fluxer_id]

    __log(f"Number of metabolite ids: {len(dict_copy_id_converter)}")
    __log(f"Number of fluxer ids: {len(dict_copy_met_info)}")
    __log("")

    return dict_copy_id_converter, dict_copy_met_info


def __dict_values_set_to_list(dict_db_names_table):
    cleaned_dict = {}
    for source_db, mapped_dbs in dict_db_names_table.items():
        cleaned_dict[source_db] = list(mapped_dbs)

    return cleaned_dict


# Main program
def __create_id_mapping_pickle():
    """
    Main function that downloads database files and processes them to merge identifiers into a mapping dictionary.
    Mapping dictionary is serialized and saved.
    """
    global dict_any_met_id_to_fluxer_id, dict_fluxer_id_to_met_prop, database_names_dict, met_prop_collected

    print("Creating directories")
    __create_directories()

    __log("Downloading files")
    tic = perf_counter()

    __download_database_files()

    toc = perf_counter()
    __log("")
    __log(f"All metabolite files downloaded in {(toc - tic) / 60:0.3f} min")

    __log("Processing metabolites")
    tic = perf_counter()

    __process_kegg_compounds(kegg_metabolites_filename)
    __process_met_file(chebi_compound_structure_filename, __chebi_compounds_inchi_reader)
    __process_met_file(chebi_compounds_filename, __chebi_compounds_names)
    __process_met_file(mx_chem_prop_filename, __metanetx_chem_prop_line_reader)
    __process_met_file(ms_met_filename, __modelseed_metabolites_line_reader)
    __process_met_file(bigg_metabolites_filename, __bigg_metabolites_line_reader)

    __process_cross_ref_info(mx_met_xref_filename, __metanetx_chem_xref_reader)
    __process_cross_ref_info(ms_met_filename, __modelseed_metabolites_xref_reader)
    __process_cross_ref_info(ms_met_aliases_filename, __modelseed_met_aliases_reader)
    __process_cross_ref_info(mx_chem_prop_filename, __metanetx_chem_prop_xref_reader)
    __process_cross_ref_info(bigg_metabolites_filename, __bigg_models_xref_reader)

    __log("Cleaning metabolite id dictionary")
    dict_any_met_id_to_fluxer_id, dict_fluxer_id_to_met_prop = __clean_id_mapping_dictionary(dict_any_met_id_to_fluxer_id,
                                                                                             dict_fluxer_id_to_met_prop)
    __log("Creating metabolite id pickle")
    with open(pickle_dir + 'metaboliteIdMapper.p', 'wb') as file:
        dump(dict_any_met_id_to_fluxer_id, file)

    toc = perf_counter()
    __log("")
    __log(f"All metabolites processed in {(toc - tic) / 60:0.3f} min")
    __log("")
    print(f"New metabolite ID mapping table created.")

