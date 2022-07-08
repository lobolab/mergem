******************************************
Updating Database ID mapping dictionaries
******************************************

The mergem algorithm uses database ID mapping dictionaries to map model metabolite identifiers into an internal identifier
that allows for comparison and merging of metabolites. Updating the mapping dictionary involves downloading the latest
identifier information from databases such as SEED, MetaNetX, KEGG, etc and linking identifiers across databases if
they represent the same metabolite, after checking for matching properties. The update process can take a few hours and
depends on the internet connection.

The :code:`-up` argument can be used to update the ID mapping pickles on the command line:

::

    mergem -up

When importing mergem into a python script, the :code:`update_id_mapper()` function can be called to update the mappping
dictionaries as shown below:

::

    import mergem

    mergem.update_id_mapper()



