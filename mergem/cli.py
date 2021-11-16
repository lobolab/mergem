#!/usr/bin/python
"""
    Uses Fluxer dictionaries to translate models to common namespace
    and then merges the models to remove metabolite and reaction
    duplicates.
    Resulting model can be FBA simulated.

    Copyright (c) Lobo Lab (https://lobolab.umbc.edu)
"""
import click
from . import __version
from . import __modelHandling
from .mergeModels import merge, update_id_mapper
import sys
import os

allowed_file_formats = ["sbml", "xml", "mat", "m", "matlab", "json", "yaml"]


@click.command()
@click.argument('input_filenames', nargs=-1, type=click.Path(exists=True))
@click.option('-obj', nargs=1, default='merge',
              help="Set objective: 'merge' all objectives (default) or 1, 2, 3.. (objective from one of the input models)")
@click.option('-o', nargs=1, help='Save model as (filename with format .xml, .sbml, etc.)')
@click.option('-v', help='Print merging statistics', is_flag=True)
@click.option('-up', help='Update ID mapping table', is_flag=True)
@click.version_option(__version.version + "\nLobo Lab (https://lobolab.umbc.edu)")
def main(input_filenames, obj, o=None, v=False, up=False):
    model_filenames = input_filenames
    objective = obj
    output_filename = o
    print_stats = v
    input_list_of_models = []

    click.secho(f"mergem, {__version.version}")
    click.secho("Lobo Lab (https://lobolab.umbc.edu)")

    if up:
        click.secho('Updating Metabolite ID mapper. This process may take a few hours.. ')
        update_id_mapper()
        click.secho('Metabolite ID mapper updated. ', fg='green')
        if len(i) == 0:
            sys.exit()

    if len(model_filenames) < 2:
        click.secho('Error: Enter 2 or more models to merge.', fg='red')
        sys.exit()

    if (objective != 'merge') and (not(int(obj) <= len(model_filenames))):
        click.secho('Error: Invalid objective selected for merged model.', fg='red')
        sys.exit()

    if output_filename is not None:
        file_format = os.path.splitext(output_filename)[1][1:].strip().lower()
        if file_format not in allowed_file_formats:
            click.secho('Error: Invalid output file format.', fg='red')
            sys.exit()

    for filename in model_filenames:
        try:
            input_model = __modelHandling.__load_model(filename)
        except Exception as e:
            click.secho(e, fg='red')
            sys.exit()
        input_list_of_models.append(input_model)

    merge_results = merge(input_list_of_models, objective)
    result_merged_model = merge_results['merged_model']

    try:
        if output_filename is None:
            output_filename = result_merged_model.id + ".xml"
        __modelHandling.__export_merged_model(result_merged_model, output_filename)

    except Exception as e:
        click.secho(e, fg='red')
        sys.exit()

    if print_stats:
        click.echo("Model (Met, Reac) jaccard distances: {}".format(merge_results['jacc_d']))
        click.echo("Mets merged: {}". format(merge_results['Met_merged']))
        click.echo("Reacs merged: {}".format(merge_results['Reac_merged']))

    click.secho(f"\nMerging models complete. Merged model saved as {output_filename}", fg="green")


if __name__ == "__main__":
    main()



