**********
mergem
**********

**mergem** is a Python library for merging two or more genome-scale metabolic
models.
The library is publicly available via PyPI at `<https://pypi.org/project/mergem/>`_ and can be pip
installed.
mergem can be used on the command-line and can also be imported within python scripts.
The package can take models in various COBRApy compatible formats such as SBML, JSON, etc. and even
COBRApy model objects, when the package is imported. The results of a single merge include the merged model,
jaccard distances between all pairs of models, number of metabolites and reactions merged, and lists of
models that contain each metabolite and reaction.

For each input model, mergem converts the metabolite IDs into a common namespace using a database ID mapping dictionary.
Reactions are compared using the participating metabolites (after conversion to common namespace). The metabolite ID mapping dictionary contains metabolite identifiers from various databases
such as ModelSEED, KEGG, ChEBI, and MetaNetX that have been unified per metabolite. The dictionary thus allows for model
metabolites to be compared more efficiently. The mapping dictionaries can be updated, during which the latest identifier
information is downloaded from each database and identifiers representing the same metabolite are mapped to one another.

mergem is also available in the user-friendly application `Fluxer <https://fluxer.umbc.edu>`_, which produces tidy flux
graphs that can visually compare the complete metabolic network from multiple models.
Documentation for Fluxer based merging can be found on its `tutorial page <https://fluxer.umbc.edu/tutorial>`_

.. toctree::
   :maxdepth: 15
   :caption: Contents

   installation
   usage
   output
   visualize
   update
   citation
   acknowledgements
   license
   contact






.. note::

   This project is under active development.








