"""
    Downloads metabolite information files from Kegg, Chebi, ModelSeed, MetaNetX and Bigg
    databases. Creates a dictionary mapping external database metabolite IDs
    to a univ ID and another dictionary mapping univ ID to metabolite properties.
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

curr_dir = os.path.dirname(__file__)
files_dir = os.path.join(curr_dir, "downloads/")
log_dir = os.path.join(curr_dir, "logs/")

# Database URLs
modelSeed_met_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv"
modelSeed_met_aliases_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt"
metaNetX_compounds_url = "https://ftp.vital-it.ch/databases/metanetx/MNXref/latest/chem_prop.tsv"
metaNetX_met_xref_url = "https://ftp.vital-it.ch/databases/metanetx/MNXref/latest/chem_xref.tsv"
metaNetX_met_depr_url = "https://ftp.vital-it.ch/databases/metanetx/MNXref/latest/chem_depr.tsv"
bigg_metabolites_url = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
chebi_compound_structure_url = "http://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/structures.csv.gz"
chebi_compounds_url = "http://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz"

# Metabolite filenames
ms_met_filename = files_dir + "modelSeed_compounds.tsv"
ms_met_aliases_filename = files_dir + "modelseed_compound_aliases.txt"
mx_chem_prop_filename = files_dir + "metanetx_chem_prop.tsv"
mx_met_xref_filename = files_dir + "metanetx_chem_xref.tsv"
mx_met_depr_filename = files_dir + "metanetx_chem_depr.tsv"
bigg_metabolites_filename = files_dir + "bigg_models_metabolites.txt"
kegg_metabolites_filename = files_dir
chebi_compound_st_zipped_filename = files_dir + "chebi_structures.csv.gz"
chebi_compound_structure_filename = files_dir + "chebi_structures.csv"
chebi_compounds_zipped_filename = files_dir + "chebi_compounds.tsv.gz"
chebi_compounds_filename = files_dir + "chebi_compounds.tsv"

modelSeed_reac_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/reactions.tsv"
modelSeed_reac_aliases_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
modelSeed_reac_pathways_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Pathways.txt"

metaNetX_reactions_url = "https://ftp.vital-it.ch/databases/metanetx/MNXref/latest/reac_prop.tsv"
metaNetX_reac_xref_url = "https://ftp.vital-it.ch/databases/metanetx/MNXref/latest/reac_xref.tsv"
bigg_reactions_url = "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt"

# Reactions
ms_reac_filename = files_dir + "modelSeed_reactions.tsv"
ms_reac_aliases_filename = files_dir + "modelSeed_reaction_aliases.txt"
mx_reac_prop_filename = files_dir + "metaNetX_reac_prop.tsv"
mx_reac_xref_filename = files_dir + "metaNetX_reac_xref.tsv"
bigg_reactions_filename = files_dir + "bigg_models_reactions.txt"
ms_reac_pathways_filename = files_dir + "modelSeed_reaction_pathways.txt"

url_dictionary = {ms_met_filename: modelSeed_met_url,
                  ms_met_aliases_filename: modelSeed_met_aliases_url,
                  mx_chem_prop_filename: metaNetX_compounds_url,
                  mx_met_xref_filename: metaNetX_met_xref_url,
                  mx_met_depr_filename: metaNetX_met_depr_url,
                  bigg_metabolites_filename: bigg_metabolites_url,
                  chebi_compounds_zipped_filename: chebi_compounds_url,
                  chebi_compound_st_zipped_filename: chebi_compound_structure_url,

                  ms_reac_filename: modelSeed_reac_url,
                  ms_reac_aliases_filename: modelSeed_reac_aliases_url,
                  ms_reac_pathways_filename: modelSeed_reac_pathways_url,
                  mx_reac_prop_filename: metaNetX_reactions_url,
                  mx_reac_xref_filename: metaNetX_reac_xref_url,
                  bigg_reactions_filename: bigg_reactions_url
                  }


url_dictionary_chebi = {chebi_compounds_filename: chebi_compounds_zipped_filename,
                        chebi_compound_structure_filename: chebi_compound_st_zipped_filename
                        }

# Dictionaries
met_univ_id_dict = {}
met_univ_id_prop_dict = {}
reac_univ_id_dict = {}
reac_univ_id_prop_dict = {}
list_primary_ids = set()
met_last_univ_id, reac_last_univ_id = 0, 0
start_time = datetime.now().strftime("%Y%m%d_%HH%MM")
primary_dbs = ['seed', 'metanetx', 'bigg', 'kegg', 'chebi']


def log(message):
    dt_string = datetime.now().strftime("%Y%m%d_%H:%M:%S")
    log_line = dt_string + " " + message
    print(log_line)
    log_file = open(log_dir + start_time + "_DatabasesDownloadLog.txt", "a")

    log_file.write(log_line + "\n")
    log_file.close()


def create_directories():
    if not os.path.exists(files_dir):
        os.makedirs(files_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)


def download_database_files():
    global kegg_metabolites_filename
    if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
            getattr(ssl, '_create_unverified_context', None)):
        ssl._create_default_https_context = ssl._create_unverified_context

    for filename in url_dictionary.keys():
        log("Downloading " + filename)
        url = url_dictionary[filename]
        with closing(urllib.request.urlopen(url)) as r:
            with open(filename, 'wb') as f:
                shutil.copyfileobj(r, f)

        log(filename + " downloaded.")

    for filename, zipped_filename in url_dictionary_chebi.items():
        with gzip_open(zipped_filename, 'rb') as file_in:
            with open(filename, 'wb') as file_out:
                shutil.copyfileobj(file_in, file_out)

    kegg_metabolites_filename += check_for_kegg_update()


def check_for_kegg_update():
    """
    Checks for a KEGG update and creates a filename for latest version.\n
    :return: filename for latest version
    """
    log("Checking Kegg stats. ")
    kegg_cpd_stats_response = requests.get("http://rest.kegg.jp/info/cpd")
    release_version = kegg_cpd_stats_response.text.split('\n')[1].split('Release ', 1)[1]
    filename = "kegg_" + \
               release_version.translate({ord(c): "-" for c in "/!@#$%^&*()[]{};:,.<>?\\|`~-=+"}).replace(" ", "-")

    filename = filename.split("--", 1)[0] + ".p"

    if not os.path.isfile(files_dir + filename):
        log("New Kegg version found: {}".format(release_version))
        log("Downloading new Kegg version.")
        download_kegg_compounds(filename)

    return filename


def download_kegg_compounds(filename):
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
            spl = (split_text[0]).split(':', 1)
            kegg_id = spl[1] if len(spl) > 1 else split_text[0]
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
    log(f"Kegg ids processed: {len(kegg_compounds_dict)}")

    with open(files_dir + filename, 'wb') as kegg_file:
        dump(kegg_compounds_dict, kegg_file)


def process_reac_file(file_name, file_line_reader):
    log("Processing file " + file_name)
    with open(file_name, "r") as file:
        next(file)  # skip header
        for line in file:
            reac_properties = file_line_reader(line.strip().split('\t'))
            if reac_properties is not None:
                append_reac_properties(reac_properties)

    log("Done processing file " + file_name)
    log(f"Number of reaction ids: {len(reac_univ_id_dict)}")
    log(f"Number of univ ids: {len(reac_univ_id_prop_dict)}")
    log("")


def process_met_file(file_name, file_line_reader):
    """
    Uses reader function to read lines of file and append informaiton to met properties dictionary.\n
    :param file_name: name of database file
    :param file_line_reader: reader function for database file
    """
    log("Processing file " + file_name)
    if ".csv" in file_name:
        with open(file_name, "r") as db_file:
            next(db_file)  # skip header
            for line in db_file:
                met_properties = file_line_reader(line.strip().split(','))
                if met_properties is not None:
                    append_met_properties(met_properties)
    else:
        with open(file_name, "r") as db_file:
            next(db_file)  # skip header
            for line in db_file:
                met_properties = file_line_reader(line.strip().split('\t'))
                if met_properties is not None:
                    append_met_properties(met_properties)

    log("Done processing file " + file_name)
    log(f"Number of metabolite ids: {len(met_univ_id_dict)}")
    log(f"Number of univ ids: {len(met_univ_id_prop_dict)}")
    log(f"List of primary ids: {len(list_primary_ids)}")
    log("")


def process_cross_ref_info(file_name, xref_line_reader):
    """
    Read cross-reference information from file and maps database identifiers
    :param file_name: name of database file containing cross-reference information
    :param xref_line_reader: reader function for database file
    """
    log("Processing file " + file_name)

    global met_univ_id_dict, met_univ_id_prop_dict
    xref_dict = xref_line_reader(file_name)

    for source_id, xref_list in xref_dict.items():
        for other_id in xref_list:
            source_univ_id = met_univ_id_dict[source_id]
            if other_id not in met_univ_id_dict:
                met_univ_id_dict[other_id] = source_univ_id
                met_univ_id_prop_dict[source_univ_id]['ids'] += [other_id]
            else:
                other_prop = met_univ_id_prop_dict[met_univ_id_dict[other_id]]
                source_db = source_id.rsplit(':', 1)[0]
                existing_db_mappings = [db_id.rsplit(':', 1)[0] for db_id in other_prop['ids']]
                if source_db not in existing_db_mappings:
                    merge_identifiers(source_id, other_id)

    log("Done processing file " + file_name)
    log(f"Number of metabolite ids: {len(met_univ_id_dict)}")
    log(f"Number of univ ids: {len(met_univ_id_prop_dict)}")
    log("")


def append_met_properties(met_properties):
    global list_primary_ids, met_univ_id_dict, met_univ_id_prop_dict
    univ_id = met_univ_id_dict.get(met_properties['ids'][0], maxsize)

    if univ_id == maxsize:
        global met_last_univ_id
        met_last_univ_id = met_last_univ_id + 1
        univ_id = met_last_univ_id

        met_univ_id_prop_dict[univ_id] = {'Name': [], 'ids': [], 'formula': [],
                                                 'mass': [], 'inchikey': [], 'xref_links': []}

    met_univ_id_dict[met_properties['ids'][0]] = univ_id
    univ_met_properties = met_univ_id_prop_dict[univ_id]
    list_primary_ids |= {db_id for db_id in met_properties['ids']}
    univ_met_properties['ids'] += [met_properties['ids'][0]]

    for key, value in met_properties.items():
        if (key == 'ids') and (len(value) > 1):
            for met_id in value:
                if met_id not in met_univ_id_dict:
                    met_univ_id_dict[met_id] = univ_id
                    univ_met_properties[key] += [met_id]
        else:
            add_values_to_property_list(univ_met_properties[key], value)


def append_reac_properties(reac_properties):
    global reac_univ_id_dict, reac_univ_id_prop_dict
    univ_id = min(fl_id for fl_id in (reac_univ_id_dict.get(reac_id, maxsize)
                                      for reac_id in reac_properties['ids']))

    if univ_id == maxsize:
        global reac_last_univ_id
        reac_last_univ_id = reac_last_univ_id + 1
        univ_id = reac_last_univ_id
        reac_univ_id_prop_dict[univ_id] = {'ids': [], 'Name': [],
                                                  'EC_num': [], 'Pathways': [], 'xref_links': []}
    unadded_db_ids = []
    for key, value in reac_properties.items():
        if key == 'ids':
            for other_id in value:
                other_univ_id = reac_univ_id_dict.get(other_id)
                if other_univ_id is None:
                    reac_univ_id_dict[other_id] = univ_id
                    reac_univ_id_prop_dict[univ_id]['ids'] += [other_id]
                else:
                    unadded_db_ids.append(other_id)
        else:
            add_values_to_property_list(reac_univ_id_prop_dict[univ_id][key], value)

    for db_id in unadded_db_ids:
        other_univ_id = reac_univ_id_dict.get(db_id)
        merge_identifiers(univ_id, other_univ_id, True)


def process_kegg_compounds(filename):
    log(f"Processing file {filename}")
    kegg_file = open(filename, "rb")
    kegg_compounds_dictionary = load(kegg_file)
    kegg_file.close()

    for kegg_id, properties in kegg_compounds_dictionary.items():
        ids = ["kegg:" + kegg_id]

        if len(properties['chebi']) > 0:
            ids += properties['chebi']

        property_dict = {'ids': ids,
                         'Name': properties['Name'],
                         'formula': properties['formula'],
                         'mass': properties['mass'],
                         }
        append_met_properties(property_dict)

    log(f"Done processing file {filename}")
    log(f"Number of metabolite ids: {len(met_univ_id_dict)}")
    log(f"Number of univ ids: {len(met_univ_id_prop_dict)}")
    log("")


def chebi_compounds_inchi_reader(line):
    """
    Reader function for chebi compounds. \n
    :param line: line from file
    :return: dictionary with properties from line
    """
    if 'inchikey' in (item.lower() for item in line):
        return {'ids': ['chebi:' + str(line[1])],
                'inchikey': [line[2]]
                }


def chebi_compounds_names(line):
    """
    Reader function for chebi compounds. \n
    :param line: line from file
    :return: dictionary with properties from line
    """
    if line[5] != 'null':
        return {'ids': [line[2].lower()],
                'Name': [line[5]]
                }


def metanetx_chem_prop_line_reader(line):
    """
    Reader function for metanetx compounds. \n
    :param line: line from file
    :return: dictionary with properties from line
    """
    if line[0][0] != "M":
        return None

    ids = ["metanetx:" + line[0]]

    if len(line) > 7:
        if len(line[5]) > 0:
            mass = [line[5]]
        else:
            mass = []

        if len(line[7]) > 0:
            inchikey = [line[7].split('=', 1)[1]]
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


def metanetx_chem_prop_xref_reader(file_name):
    """
    Reader function for cross references from metanetx compounds file. \n
    :param file_name: line from file
    :return: dictionary with properties from line
    """
    xref_dict = {}
    with open(file_name, "r") as db_file:
        next(db_file)  # skip header
        for line in db_file:
            if (line[0][0] == '#') or (':' not in line[2]):
                continue
            [other_db, other_db_id] = line[2].split(":", 1)
            source_id = "metanetx:" + line[0]
            other_id = None

            if other_db.lower() == "mnx":
                other_id = "metanetx:" + other_db_id

            elif other_db.lower() in {"kegge", "envipathm", "envipath", "biggm", "seedm", "reactomem", "lipidmapsm"}:
                continue

            elif other_db[:-1].lower() in {'bigg', 'kegg', 'seed', 'metacyc'}:
                other_db_name = other_db[:-1].lower()
                other_id = other_db_name + ":" + other_db_id

            elif other_db == "slm":
                other_db_name = "slm"
                other_id = other_db_name + ":" + other_db_id

            if other_id is not None:
                if source_id in xref_dict:
                    xref_dict[source_id] += [other_id]
                else:
                    xref_dict[source_id] = [other_id]
    return xref_dict


def metanetx_chem_depr_reader(file_name):
    """
    Reader function for deprecated metabolite IDs in MetaNetX. \n
    :param file_name: name of deprecated IDs file
    :return: dictionary with latest MetaNetX IDs mapped to deprecated IDs
    """
    xref_dict = {}
    version = ''
    with open(file_name, "r") as db_file:
        next(db_file)  # skip header
        for line in db_file:
            if ('#VERSION' in line) and (version == ''):
                version = line.split('\n')[0][-3:]

            if line[0][0] == '#':
                continue

            line = line.strip().split('\t')
            if line[2] == version:
                mnx_id = 'metanetx:' + line[1]
                mnx_old = 'metanetx:' + line[0]

                if mnx_id in xref_dict:
                    xref_dict[mnx_id] += [mnx_old]
                else:
                    xref_dict[mnx_id] = [mnx_old]

    return xref_dict


def metanetx_chem_xref_line_reader(line):
    """
        Reader function for metanetx compounds. \n
        :param line: line from file
        :return: dictionary with properties from line
    """
    if line[0][0] == '#':
        return None

    ids = ["metanetx:" + line[1]]
    name = [line[2]]

    if 'obsolete' not in name[0]:
        return {'ids': ids
                }


def metanetx_chem_xref_reader(file_name):
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
                    (":" not in line[0]) or ("unknown" in line[2].lower()) or ("no description" in line[2].lower()) or \
                    ("obsolete" in line[2].lower()) \
                    or ("molecular entity" in line[2].lower()):
                continue

            source_id = "metanetx:" + line[1]
            [other_db, other_db_id] = line[0].split(":", 1)
            other_id = None

            if "." in other_db:
                other_db = other_db.split(".")[0]
                other_id = other_db + ":" + other_db_id

            elif other_db not in {"keggC", "envipathM", "envipath", "seedM", "CHEBI", "biggM", "keggD", "SLM", "keggE",
                                  "keggG", "reactomeM", "sabiorkM", "rheaP", "rheaG", "lipidmapsM", "metacycM"}:
                other_id = other_db + ":" + other_db_id

            if other_id is not None:
                if source_id in xref_dict:
                    xref_dict[source_id] += [other_id]
                else:
                    xref_dict[source_id] = [other_id]

    return xref_dict


def modelseed_metabolites_line_reader(line):
    """
    Line reader function for modelseed metabolites. \n
    :param line: line from file
    :return: dictionary of metabolite information from line
    """
    ids = ["seed:" + line[0]]

    if (len(line[4]) > 0) and (line[4].lower() != 'none') and (line[4].lower() != 'null')\
            and (float(line[4]) != 10000000):
        mass = [line[4]]
    else:
        mass = []

    return {'ids': ids,
            'Name': [line[2]],
            'formula': [line[3]],
            'mass': mass,
            'inchikey': [line[6]]
            }


def modelseed_metabolites_xref_reader(file_name):
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
                ids += [('seed:' + ele) for ele in (line[10].split(";"))]

            if len(ids) > 0:
                if source_id in xref_dict:
                    xref_dict[source_id] += [seed_id for seed_id in ids]
                else:
                    xref_dict[source_id] = [seed_id for seed_id in ids]
    return xref_dict


def modelseed_met_aliases_reader(file_name):
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
                other_id = other_db_name + ":" + line[1]

            elif other_db in primary_dbs:
                other_db_name = other_db
                other_id = other_db_name + ":" + line[1]

            else:
                continue

            if source_id in xref_dict:
                xref_dict[source_id] += [other_id]
            else:
                xref_dict[source_id] = [other_id]

    return xref_dict


def bigg_metabolites_line_reader(line):
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
    ids += ['bigg:' + line[1]]

    if len(line) > 2:
        names += [line[2]]

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


def bigg_models_xref_reader(file_name):
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
                        ids += [other_db_id.lower()]
                    else:
                        other_db_name = other_db.split('.')[0].lower()
                        ids += [other_db_name + ':' + other_db_id]

            if source_id in xref_dict:
                xref_dict[source_id] += [db_id for db_id in ids]
            else:
                xref_dict[source_id] = [db_id for db_id in ids]
    return xref_dict


def modelSeed_reactions_line_reader(line):
    """
    Reader function for collecting reaction info from modelseed file. \n
    :param line: line from file
    :return: dictionary mapping id to properties
    """
    ids = ["seed:" + line[0]]

    if line[19] != 'null':
        ids += [('seed:' + ele) for ele in (line[19].split(";"))]

    return {'ids': ids,
            'Name': [line[2]],
            'EC_num': [line[13]]}


def modelSeed_reaction_aliases_line_reader(line):
    """
    Reader function for collecting identifiers from modelseed aliases file. \n
    :param line: line from file
    :return: dictionary with database ids
    """
    other_db = line[2].lower()
    if other_db == 'metanetx.reaction':
        return {'ids': ["seed:" + line[0], "metanetx:" + line[1]]}

    elif other_db in {'bigg', 'bigg1'}:
        return None

    elif other_db.startswith("i"):
        return None

    else:
        return {'ids': ["seed:" + line[0], other_db + ":" + line[1]]}


def modelSeed_reaction_pathways_line_reader(line):
    """
    Reader function for collecting pathway information. \n
    :param line: line from modelseed pathways file
    :return: dictionary mapping seed id to pathway
    """
    return {'ids': ["seed:" + line[0]],
            'Pathways': [line[1]]}


def metanetx_reaction_prop_line_reader(line):
    """
    Reader function for cross references from metanetx reaction file. \n
    :param file_name: line from file
    :return: dictionary with properties from line
    """

    if line[0][0] != "M":
        return None

    [other_db, other_db_id] = line[2].split(":")
    if other_db.lower() == "mnx":
        other_id = "metanetx:" + other_db_id
    else:
        if line[2].split(":")[0][-1] == "R":
            other_db_name = other_db[:-1].lower()
        else:
            other_db_name = other_db
        other_id = other_db_name + ":" + other_db_id

    ids = ["metanetx:" + line[0], other_id]

    if len(line) > 3:
        return {'ids': ids,
                'EC_num': [line[3]]}
    else:
        return {'ids': ids}


def metanetx_reaction_xref_line_reader(line):
    """
    Reader function for metanetx reaction cross reference file. \n
    :param file_name: name of database file
    :return: dictionary with properties from line
        """
    if (line[0][0] == "#") or (line[0][0:3] == "MNX") or (line[0][0:3] == "mnx") or (line[1] == "EMPTY"):
        return None

    [other_db, other_db_id] = line[0].split(":")
    if "." in other_db:
        other_db_name = other_db.split(".")[0].lower()
        other_id = other_db_name + ":" + other_db_id
    else:
        if line[0].split(":")[0][-1] == "R":
            other_db_name = other_db[:-1].lower()
        else:
            other_db_name = other_db
        other_id = other_db_name + ":" + other_db_id

    ids = ["metanetx:" + line[1], other_id]

    return {'ids': ids}


def bigg_reactions_line_reader(line):
    """
    Line reader for reactions from BiGG database.\n
    :param line: line in file
    :return: dictionary with reaction information from line
    """
    ids = []
    ec_nums = []
    xref_links = []
    names = []

    ids += ['bigg:' + line[0].lower()]
    ids += [('bigg:' + ele) for ele in (line[5].split("; "))]

    names += [line[1]]

    if "http" in line[4]:
        for link in line[4].split('; '):
            link_part = link.split('/')
            other_db = link_part[3]
            other_db_id = link_part[4]

            if other_db.lower() == "ec-code":
                ec_nums += [other_db_id]
            else:
                other_db_name = other_db.split('.')[0].lower()
                if other_db_name in {"seed", "metanetx"}:
                    ids += [other_db_name + ':' + other_db_id]

            link_lower = link.lower()
            if ("seed" in link_lower) or ("kegg" in link_lower) or ("metanetx" in link_lower):
                xref_links += ["http" + link.split("http")[1]]

    prop = {}
    if ids:
        prop["ids"] = ids
    if ec_nums:
        prop["EC_num"] = ec_nums
    if xref_links:
        prop["xref_links"] = xref_links
    if names:
        prop["Name"] = names

    if prop:
        return prop
    else:
        return None


def add_values_to_property_list(property_list, prop_value, for_name=False):
    invalid_values_list = ['\'\'', '\"\"', 'null', '-', '']
    for value in prop_value:
        if not for_name:
            if (value not in invalid_values_list) and (value not in property_list):
                property_list.append(value)
        else:
            property_list_lower = [prop.lower() for prop in property_list]
            if (value not in invalid_values_list) and (value.lower() not in property_list_lower):
                property_list.append(value)


def merge_identifiers(source_id, other_id, for_reac=False):
    """
    Merges the two met IDs into lowest univ ID.
    :param source_id: primary metabolite ID from database being processed
    :param other_id: cross referenced metabolite ID to be mapped to primary met ID
    """
    global met_univ_id_prop_dict, met_univ_id_dict, \
        reac_univ_id_prop_dict, reac_univ_id_dict

    if for_reac:
        source_univ_id = source_id
        other_univ_id = other_id
        id_mapper = reac_univ_id_dict
        prop_mapper = reac_univ_id_prop_dict

    else:
        source_univ_id = met_univ_id_dict.get(source_id)
        other_univ_id = met_univ_id_dict.get(other_id)
        id_mapper = met_univ_id_dict
        prop_mapper = met_univ_id_prop_dict

    source_id_properties = prop_mapper[source_univ_id]
    other_properties = prop_mapper[other_univ_id]

    if other_univ_id != source_univ_id:
        if source_univ_id < other_univ_id:
            for xref_id in other_properties['ids']:
                id_mapper[xref_id] = source_univ_id
            for key, value in other_properties.items():
                if key == 'Name':
                    add_values_to_property_list(source_id_properties[key], value, True)
                else:
                    add_values_to_property_list(source_id_properties[key], value)

            del prop_mapper[other_univ_id]

        else:
            for xref_id in source_id_properties['ids']:
                id_mapper[xref_id] = other_univ_id
            for key, value in source_id_properties.items():
                if key == 'Name':
                    add_values_to_property_list(other_properties[key], value, True)
                else:
                    add_values_to_property_list(other_properties[key], value)

            del prop_mapper[source_univ_id]


def clean_id_mapping_dictionary(id_dictionary, info_dictionary, for_reac=False):
    dict_copy_id_converter = {}
    dict_copy_info = {}
    db_name_dict = {}
    if for_reac:
        db_preference = {'metanetx': 0, 'seed': 1, 'bigg': 2, 'kegg': 3, 'sabiork': 4, 'metacyc': 5}
    else:
        db_preference = {'kegg': 0, 'chebi': 1, 'metanetx': 2, 'bigg': 3, 'seed': 4, 'sabiork': 5}

    for db_id, fl_id in id_dictionary.items():
        [db_name, new_key] = db_id.split(":", 1)
        conflict_fl_id = dict_copy_id_converter.get(new_key)
        if conflict_fl_id:
            if db_preference.get(db_name, maxsize) > db_preference.get(db_name_dict[new_key], maxsize):
                remove_conflicting_id(new_key, fl_id, info_dictionary)
                fl_id = conflict_fl_id
            else:
                remove_conflicting_id(new_key, conflict_fl_id, info_dictionary)

        dict_copy_id_converter[new_key] = fl_id
        if new_key not in db_name_dict:
            db_name_dict[new_key] = db_name

    for univ_id in dict_copy_id_converter.values():
        copied_info = info_dictionary[univ_id].copy()
        copied_info['ids'] = list(set(copied_info['ids']))
        dict_copy_info[univ_id] = copied_info

    log(f"Number of database ids: {len(dict_copy_id_converter)}")
    log(f"Number of univ ids: {len(dict_copy_info)}")
    log("")

    return dict_copy_id_converter, dict_copy_info


def remove_conflicting_id(db_id, fl_id_to_update, info_dict):
    other_db_id = []
    if len(info_dict[fl_id_to_update]['ids']) > 1:
        for other_id in info_dict[fl_id_to_update]['ids']:
            if db_id == other_id.split(":", 1)[1]:
                other_db_id.append(other_id)
        for d_id in other_db_id:
            info_dict[fl_id_to_update]['ids'].remove(d_id)

    if len(info_dict[fl_id_to_update]['ids']) == 0:
        del info_dict[fl_id_to_update]


# Main program
def build_id_mapping(delete_database_files):
    """
    Main function that downloads database files and processes them to merge identifiers into a mapping dictionary.
    Mapping dictionary is serialized and saved.
    """
    global met_univ_id_dict, met_univ_id_prop_dict, \
        reac_univ_id_dict, reac_univ_id_prop_dict

    print("Creating directories")
    create_directories()

    log("Downloading files")
    tic = perf_counter()

    download_database_files()

    toc = perf_counter()
    log("")
    log(f"All files downloaded in {(toc - tic) / 60:0.3f} min")

    log("Processing metabolites")
    tic = perf_counter()

    # Process KeGG database metabolites
    process_kegg_compounds(kegg_metabolites_filename)

    # Process ChEBI database metabolite IDs
    process_met_file(chebi_compound_structure_filename, chebi_compounds_inchi_reader)
    process_met_file(chebi_compounds_filename, chebi_compounds_names)

    # Process MetaNetX database metabolite IDs
    process_met_file(mx_chem_prop_filename, metanetx_chem_prop_line_reader)
    process_cross_ref_info(mx_met_depr_filename, metanetx_chem_depr_reader)
    process_met_file(mx_met_xref_filename, metanetx_chem_xref_line_reader)
    process_cross_ref_info(mx_met_xref_filename, metanetx_chem_xref_reader)

    # Process BiGG database metabolite IDs
    process_met_file(bigg_metabolites_filename, bigg_metabolites_line_reader)
    process_cross_ref_info(bigg_metabolites_filename, bigg_models_xref_reader)

    # Process ModelSEED database metabolite IDs
    process_met_file(ms_met_filename, modelseed_metabolites_line_reader)
    process_cross_ref_info(ms_met_filename, modelseed_metabolites_xref_reader)
    process_cross_ref_info(ms_met_aliases_filename, modelseed_met_aliases_reader)

    # Process reaction IDs from modelSEED, MetaNetX, and BiGG
    process_reac_file(ms_reac_filename, modelSeed_reactions_line_reader)
    process_reac_file(ms_reac_pathways_filename, modelSeed_reaction_pathways_line_reader)
    process_reac_file(mx_reac_prop_filename, metanetx_reaction_prop_line_reader)
    process_reac_file(mx_reac_xref_filename, metanetx_reaction_xref_line_reader)
    process_reac_file(bigg_reactions_filename, bigg_reactions_line_reader)

    log("Cleaning reaction id mapping dictionary")
    reac_univ_id_dict, reac_univ_id_prop_dict = clean_id_mapping_dictionary(reac_univ_id_dict,
                                                                              reac_univ_id_prop_dict,
                                                                              for_reac=True)

    log("Cleaning metabolite id dictionary")
    met_univ_id_dict, met_univ_id_prop_dict = clean_id_mapping_dictionary(met_univ_id_dict,
                                                                            met_univ_id_prop_dict)

    if(delete_database_files):
        log("Deleting downloads")
        shutil.rmtree(files_dir)

    toc = perf_counter()
    log("")
    log(f"All database identifiers processed in {(toc - tic) / 60:0.3f} min")
    log("")
    print(f"New ID mapping tables created.")

    return met_univ_id_dict, met_univ_id_prop_dict, reac_univ_id_dict, reac_univ_id_prop_dict