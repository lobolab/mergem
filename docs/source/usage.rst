
Using mergem to merge models
=================================
mergem can merge two or more genome-scale metabolic models. The command-line execution can take input models in
various COBRApy compatible formats (SBML, JSON, YAML, and MAT).
mergem can be imported into a python script and the merge function can take cobra objects in addition to filenames.
A single objective function from any of the input models can be set as the objective for merged model. Alternatively,
objective functions from all input models can be merged into a single function and set as the objective in the merged
model.


.. _cli:

Command-line execution
--------------------------
Once mergem has been installed using pip, the following commands can be run on the command-line.
Printing help text displays all the options.

.. code-block:: none

    > mergem --help
    Usage: mergem [OPTIONS] [INPUT_FILENAMES]...

    mergem takes genome-scale metabolic models as input, merges them into a
    single model and saves merged model as .xml. Users can optionally select the
    objective and provide an output filename for merged model.

    Options:
    -obj TEXT  Set objective: 'merge' all objectives (default) or 1, 2, 3..
             (objective from one of the input models)
    -o TEXT    Save model as (filename with format .xml, .sbml, etc.)
    -v         Print merging statistics
    -up        Update ID mapping table
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

    mergem model1.xml model2.xml -obj 1 -v


Output model filename can be provided using the :code:`-o` argument followed by desired output filename with file format
specified in the extension (.xml, ...):

::

    mergem model.xml model2.xml -obj 1 -o mergedmodel.xml


.. _python-import:

Python import
---------------------

Import the mergem package to use its modules within a python script:

::

    import mergem


Provide the list of models to be merged:

::

    merge_results = mergem.merge([model1, model2, ..], objective)

where the models can be COBRApy model objects, or filenames and objective can be 'merge' or
model number ('1', 2', '3', etc.).
