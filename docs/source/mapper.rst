*******************************
Save ID mapping tables
*******************************

Saving mappers from the command-line
========================================

The database ID mapping tables can be downloaded as csv files. Each file contains all the
universal IDs for either metabolites or reactions followed by the corresponding cross-referenced
identifiers in other databases.

The ID mapping tables are saved as :file:`mergem_univ_id_mapper_metabolites.csv` and
:file:`mergem_univ_id_mapper_reactions.csv` in the working directory.

The :code:`-s` argument can be used to save the ID mapping tables on the command line:

::

    mergem -s

Saving mappers using Python script
========================================

When importing mergem into a python script, the :code:`save_mapping_tables()` function can be called to
save the mapping tables as shown below:

::

    import mergem

    mergem.save_mapping_tables()


Custom filenames for saving ID mappers
------------------------------------------

The mapping table filenames can be customized by providing the desired filenames as input to the
function.

::

   mergem.save_mapping_tables('custom_met_mapper_filename.csv', 'custom_reac_mapper_filename.csv')