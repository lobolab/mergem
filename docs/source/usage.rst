*********************************
Using mergem to merge models
*********************************
mergem can merge two or more genome-scale metabolic models. The command-line execution can take input models in
various COBRApy compatible formats (SBML, JSON, YAML, and MAT).
mergem can be imported into a python script and the merge function can take cobra objects in addition to filenames.
A single objective function from any of the input models can be set as the objective for merged model. Alternatively,
objective functions from all input models can be merged into a single function and set as the objective in the merged
model.


.. _cli:

Command-line
==========================
Once mergem has been installed using pip, the following commands can be run on the command-line.
Printing help text displays all the options.

::

    > mergem --help
    Usage: mergem [INPUT_FILENAMES] [OPTIONS]

    mergem takes genome-scale metabolic models as input, merges them into a
    single model and saves merged model as .xml. Users can optionally select the
    objective and provide an output filename for merged model.

    Options:
    -obj TEXT  Set objective: 'merge' all objectives (default) or 1, 2, 3...
             (objective from one of the input models)
    -o TEXT    Save model as (filename with format .xml, .sbml, etc.)
    -v         Print merging statistics
    -up        Update ID mapping table
    -s         Save ID mapping table as CSV
    -e         Uses exact stoichiometry when merging reactions
    -p         Consider protonation when merging reactions
    -a         Extend annotations with mergem database of metabolites and reactions
    -t         Translate metabolite and reaction IDs to a target namespace (chebi, metacyc, kegg, reactome, metanetx, hmdb, biocyc, bigg, seed, sabiork, or rhea)
    --version  Show the version and exit.
    --help     Show this message and exit.


For merging two models and setting objective of merged model from first model, use:

::

    mergem model1.xml model2.xml


The :code:`-obj` argument can be used to set the objective function of merged model. Allowed values include :code:`merge`
to merge input model objective functions (default) or an integer representing the objective function from the model
in order of input (1, 2, 3, ..):

::

    mergem model1.xml model2.xml -obj 1


To print merge statistics, append the :code:`-v` argument:

::

    mergem model1.xml model2.xml -v


Output model filename can be provided using the :code:`-o` argument followed by desired output filename with file format
specified in the extension (.xml, ...):

::

    mergem model.xml model2.xml -o mergedmodel.xml


Save the ID mapping table as a CSV file by using the :code:`-s` argument:

::

    mergem model1.xml model2.xml -s


By default, reactions are merged when they have both a similar set of reactants and a similar set of products, without comparing their stoichiometry. To merge reactions only when they have the same exact stoichiometry in their reactants and products, use the :code:`-e` argument:

::

    mergem model1.xml model2.xml -e


By default, reactions are compared ignoring the hydrogen and proton metabolites. To consider also the hydrogen and proton metabolites when comparing reactions, use the :code:`-p` argument:

::

    mergem model1.xml model2.xml -p


Metabolite and reaction annotations are merged from all input models. In addition, mergem can extend these annotations using the mergem database. For extending the annotations using mergem dabase, use the :code:`-a` argument:

::

    mergem model1.xml model2.xml -a


Mergem can translate the metabolite and reaction IDs to another database system when using the :code:`-t` argument:

::

    mergem model1.xml -t chebi


.. _python-import:

Python
=======================

Merge models
-----------------

Import the mergem package to use its modules within a python script:

::

    import mergem


Provide the list of models to be merged:

::

    merge_results = mergem.merge(input_models, set_objective='merge', exact_sto=False)

* :code:`input_models` is a list of COBRApy model objects or strings specifying file names.
* :code:`set_objective` specifies if the objective functions are merged ('merge') or copied from a single model (specifying the index of the model: '1', 2', '3', etc.).
* :code:`exact_sto` is set to true for using exact stoichiometry when merging reactions.


Other mergem functions
---------------------------

The following functions can also be imported from mergem:

::

    from mergem import load_model, save_model, map_localization, map_metabolite_univ_id, map_reaction_univ_id,
                        get_metabolite_properties, get_reaction_properties, update_id_mapper


:code:`load_model` loads a model from the given filename/path.

:code:`save_model` takes a cobra model as input and exports as file_name.

:code:`map_localization` converts localization suffixes into common notation.

:code:`map_metabolite_univ_id` maps metabolite id to metabolite universal id.

:code:`map_reaction_univ_id` maps reaction id to metabolite universal id.

:code:`get_metabolite_properties` retrieves the properties of a metabolite using its universal id

:code:`get_reaction_properties` retrieves the properties of a reaction using its universal id

:code:`update_id_mapper` downloads the latest database files, merges the database identifiers based on common properties and saves the mapping tables as pickles.



All the functions can be imported at once with:

::

    from mergem import *


