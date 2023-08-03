mergem
======
mergem is a python package for merging genome-scale metabolic models.
The package can be used as a command-line tool or can be imported within a python script.

------


Installation
------
To install the latest release  

    pip install mergem

------

Usage
------
For detailed usage instructions, please refer to the help [documentation](https://mergem.readthedocs.io/en/latest/).

#### Command-line usage
Command-line options can be viewed using "--help" flag, as shown below:

    > mergem --help
    Usage: mergem [OPTIONS]

    Options:
    -i TEXT    Input model filenames
    -obj TEXT  Set objective: 'merge' all objectives (default) or 1, 2, 3..
             (objective from one of the input models)
    -o TEXT    Save model as (filename with format .xml, .sbml, etc.)
    -v         Print merging statistics
    -up        Update metabolite ID mapping table
    -s         Save ID mapping table as CSV
    -e         Uses exact stoichiometry for merging
    --version  Show the version and exit.
    --help     Show this message and exit.
 
For merging two models and setting objective of merged model from first model, use:

    mergem -i model1.xml -i model2.xml -obj 1

To print merging statistics, append the "-v" flag:

    mergem -i model1.xml -i model2.xml -obj 1 -v 

#### Importing mergem

To use mergem modules within a python script, simply import the package within the script:

    import mergem

Provide the list of models to be merged as a list to the merge function:

    merge_results = mergem.merge([model1, model2,..], objective)

where objective can be 'merge' or model index ('1', '2', '3', etc).
The merge function returns a dictionary of results including the merged model,
the metabolite and reaction jaccard distances of each model with respect to first model, and the 
metabolite and reaction model sources. 

------
Citation
======
Please cite mergem using: <br>
<br> [mergem: merging and comparing genome-scale metabolic models using universal identifiers](https://doi.org/10.1101/2022.07.14.499633)
<br> A. Hari, D. Lobo.
<br> <b>bioRxiv</b>, doi:10.1101/2022.07.14.499633, 2022

------

------
Acknowledgements 
======

This package was developed at [The Lobo Lab](https://lobolab.umbc.edu), University of Maryland Baltimore County.

------

References
======
The following publications have contributed towards the development of this package:
* [Fluxer: a web application to compute and visualize genome-scale metabolic flux networks](https://doi.org/10.1093/nar/gkaa409).
* [COBRApy: COnstraints-Based Reconstruction and Analysis for Python](http://dx.doi.org/doi:10.1186/1752-0509-7-74).

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