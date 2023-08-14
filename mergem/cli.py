#!/usr/bin/python
"""
    Uses mergem dictionaries to translate models to common namespace
    and then merges the models to remove metabolite and reaction
    duplicates.
    Resulting model can be FBA simulated.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""
import click
import sys
import os
import mergem
from .__version import _version

_allowed_file_formats = ["sbml", "xml", "mat", "m", "matlab", "json", "yaml"]


@click.command(no_args_is_help=True)
@click.argument('input_filenames', nargs=-1, type=click.Path(exists=True))
@click.option('-obj', nargs=1, default='merge',
              help="Set objective: 'merge' all objectives (default) or 1, 2, 3.. (objective from one of the input models)")
@click.option('-o', nargs=1, help='Save model as (filename with format .xml, .sbml, etc.)')
@click.option('-v', help='Print merging statistics', is_flag=True)
@click.option('-up', help='Update ID mapping table', is_flag=True)
@click.option('-s', help='Save ID mapping table as CSV', is_flag=True)
@click.option('-e', help='Uses exact stoichiometry when merging reactions', is_flag=True)
@click.version_option(_version + "\nLobo Lab (https://lobolab.umbc.edu)")
def main(input_filenames, obj, o=None, v=False, up=False, s=False, e=False):
    """
    mergem takes genome-scale metabolic models as input, merges them into a single model
    and saves the merged model as .xml. Users can optionally select the objective and provide
    an output filename for the merged model.

    Lobo Lab (https://lobolab.umbc.edu)
    """
    model_filenames = input_filenames
    objective = obj
    output_filename = o
    print_stats = v
    input_list_of_models = []

    click.secho(f"mergem, v{_version}")

    if up:
        click.secho('Updating ID mapper. This process may take a few hours.. ')
        mergem.update_id_mapper()
        click.secho('ID mapper updated. ', fg='green')
        if len(model_filenames) == 0:
            sys.exit()

    if s:
        click.secho('Saving ID mapper as CSV file.  ')
        mergem.save_mapping_tables()
        click.secho('ID mapping tables saved. ', fg='green')
        if len(model_filenames) == 0:
            sys.exit()

    if len(model_filenames) < 2:
        click.secho('Error: Enter 2 or more models to merge.', fg='red')
        sys.exit()

    if (objective != 'merge') and (not(int(obj) <= len(model_filenames))):
        click.secho('Error: Invalid objective selected for merged model.', fg='red')
        sys.exit()

    if output_filename is not None:
        file_format = os.path.splitext(output_filename)[1][1:].strip().lower()
        if file_format not in _allowed_file_formats:
            click.secho('Error: Invalid output file format.', fg='red')
            sys.exit()

    for filename in model_filenames:
        try:
            input_model = mergem.load_model(filename)
        except Exception as e:
            click.secho(e, fg='red')
            sys.exit()
        input_list_of_models.append(input_model)

    merge_results = mergem.merge(input_list_of_models, set_objective=objective, exact_sto=e)
    result_merged_model = merge_results['merged_model']

    try:
        if output_filename is None:
            output_filename = result_merged_model.id + ".xml"
        mergem.save_model(result_merged_model, output_filename)

    except Exception as e:
        click.secho(e, fg='red')
        sys.exit()

    click.secho(f"\nMerging models complete. Merged model saved as {output_filename}", fg="green")

    if print_stats:
        click.echo("Jaccard distance matrix: {}".format(merge_results['jacc_matrix']))
        click.echo("Metabolites merged: {}". format(merge_results['num_met_merged']))
        click.echo("Reactions merged: {}".format(merge_results['num_reac_merged']))


if __name__ == "__main__":
    main()



