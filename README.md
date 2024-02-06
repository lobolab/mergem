mergem
======
mergem is a python package and command-line tool for merging, comparing, and translating genome-scale metabolic models.

------


Installation
------
Use pip to install the latest release:

    pip install mergem

------

Usage
------
For detailed usage instructions, please refer to the [documentation](https://mergem.readthedocs.io/en/latest/).

#### Command-line usage
Command-line options can be viewed using "--help" flag, as shown below:

    > mergem --help
    Usage: mergem [OPTIONS] [INPUT_FILENAMES]...

    mergem takes genome-scale metabolic models as input, merges them into a single model and saves the merged model as .xml. Users can optionally select the objective, provide an output filename for the merged model, and translate the models to a different namespace.

    Lobo Lab (https://lobolab.umbc.edu)

    Options:
    -obj TEXT  Set objective: 'merge' all objectives (default) or 1, 2, 3... (objective from one of the input models)  
    -o TEXT    Save model as (filename with format .xml, .sbml, etc.)  
    -v         Print merging statistics
    -up        Update ID mapping table
    -s         Save ID mapping table as CSV
    -e         Uses exact stoichiometry when merging reactions
    -p         Consider protonation when merging reactions
    -a         Extend annotations with mergem database of metabolites and reactions
    -t TEXT    Translate all metabolite and reaction IDs to a target namespace (chebi, metacyc, kegg, reactome, metanetx, hmdb, biocyc, bigg, seed, sabiork, or rhea)
    --version  Show the version and exit.
    --help     Show this message and exit.

 
For merging two models and setting objective of merged model from first model, use:

    mergem -i model1.xml -i model2.xml -obj 1

To print merging statistics, append the "-v" flag:

    mergem -i model1.xml -i model2.xml -obj 1 -v 

#### Python usage

To use mergem  within a python script, simply import the package with:

    import mergem

For merging or processing one, two, or more models, provide a list of models to the merge function:

    results = mergem.merge(input_models, set_objective='merge', exact_sto=False use_prot=False, extend_annot=False, trans_to_db=None)
    merged_model = results['merged_model']
    jacc_matrix = results['jacc_matrix']
    num_met_merged = results['num_met_merged']
    num_reac_merged = results['num_reac_merged']
    met_sources = results['met_sources']
    reac_sources = results['reac_sources']

* `input_models` is a list of one or more COBRApy model objects or strings specifying file names.
* `set_objective` specifies if the objective functions are merged ('merge') or copied from a single model (specifying the index of the model: '1', 2', '3', etc.).
* `exact_sto` use exact stoichiometry when merging reactions.
* `use_prot` consider hydrogen and proton metabolites when merging reactions.
* `add_annot` add additional metabolite and reaction annotations from mergem dictionaries.
* `trans_to_db` translate metabolite and reaction IDs to a target database (chebi, metacyc, kegg, reactome, metanetx, hmdb, biocyc, bigg, seed, sabiork, or rhea)

* `results` a dictionary with all the results, including:
* `merged_model` the merged model.
* `jacc_matrix` metabolite and reaction jaccard distances.
* `num_met_merged` number of metabolites merged.
* `num_reac_merged` number of reactions merged.
* `met_sources` dictionary mapping each metabolite ID in the merged model to the corresponding metabolite IDs from each of the input models.
* `reac_sources` dictionary mapping each reaction ID in the merged model to the corresponding reaction IDs from each of the input models.

The merge function returns a dictionary of results including the merged model,
the metabolite and reaction Jaccard distance matrix between models, and the 
metabolite and reaction model sources. 


The following functions can also be imported from mergem:

    from mergem import translate, load_model, save_model, map_localization, map_metabolite_univ_id, map_reaction_univ_id, get_metabolite_properties, get_reaction_properties, update_id_mapper

* `translate(input_model, trans_to_db)` translates a model to another target database specified in `trans_to_db`.
* `load_model(filename)` loads a model from the given filename/path.
* `save_model(cobra_model, file_name)` takes a cobra model as input and exports it as file `file_name`.
* `map_localization(id_or_model_localization)` converts localization suffixes into common notation.
* `map_metabolite_univ_id(met_id)` maps metabolite id to metabolite universal id.
* `map_reaction_univ_id(reac_id)` maps reaction id to metabolite universal id.
* `get_metabolite_properties(met_univ_id)` retrieves the properties of a metabolite using its universal id
* `get_reaction_properties(reac_univ_id)` retrieves the properties of a reaction using its universal id
* `update_id_mapper(delete_database_files)` updates and build mergem database. It will download the latest source database files, merge the identifiers based on common properties, and save the mapping mapping tables and information internally. This process can take several hours. The parameter specifies if the downloaded intermediate database files are deleted after the update (saves disk space but the next update will take longer; dafault is True).


------
Citation
======
Please cite mergem using: <br>
<br> [mergem: merging, comparing, and translating genome-scale metabolic models using universal identifiers](https://doi.org/10.1093/nargab/lqae010)
<br> A. Hari, A. Zarrabi, D. Lobo
<br> <b>NAR Genomics and Bioinformatics</b>, 6(1), lqae010, 2024

------
Acknowledgements 
======

This package was developed at [The Lobo Lab](https://lobolab.umbc.edu), University of Maryland, Baltimore County.

------

License
======
This package is under GNU GENERAL PUBLIC LICENSE. The package is free for use without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the
use of this software. Permission is granted to anyone to use this software for any purpose, 
subject to the following restrictions:

1. The origin of this software and database must not be misrepresented;
   you must not claim that you wrote the original software.
2. If you use this software and/or database in a work (any production in the scientific, literary, and artistic domain), 
   an acknowledgment and citation (see publication above) in the work is required.
3. This notice may not be removed or altered from any distribution.