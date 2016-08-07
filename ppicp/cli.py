#!/usr/bin/env python

"""
PPICP -- Protein-Protein Interaction Calculation Pipeline

A scientific software for the calculation of protein-protein interactions (PPIs).
The pipeline runs various third party software to compute PPIs.
It consists of four major steps:
    1. Downloading of input data based on user input (data from PDB, DSSP, and Reduce).
    2. Calculation of PPIs.
    3. Motif detection in PPIs-graphs.
    4. Compiling of statistics based on the previous output.

For more information about the third party software and the pipeline itself please refer to the
documentation.
"""

from __future__ import division

import argparse
import collections
import os
import shutil
import time

import config
import dssp
import hydrogen
import initialize
import output_results
import pdb
import plcc
import statistics
import utilities


def create_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    cli_parser = argparse.ArgumentParser(description="PPICP -- Protein-Protein Interaction "
                                                     "Calculation Pipeline",
                                         epilog="A scientific software for the calculation of "
                                                "protein-protein interactions (PPIs).")

    cli_parser.add_argument('input', type=lambda x: initialize.check_input_path(x, cli_parser),
                            nargs=1, help="Path to a file listing PDB-IDs, to a PDB file, or to a "
                                          "directory containing PDB files.")

    cli_parser.add_argument('output', type=lambda x: initialize.check_output_path(x), nargs=1,
                            help="Path to an output directory.")
    cli_parser.add_argument('--ppi', '-p', choices=['all', 'with-ligands', 'no-pi-effects',
                                                    'ligands-no-pi-effects'])
    cli_parser.add_argument('--pdb', '-pdb', action='store_true')
    cli_parser.add_argument('--dssp', '-d', action='store_true')
    cli_parser.add_argument('--hydrogen', '-r', action='store_true')
    cli_parser.add_argument('--models2chains', '-s', action='store_true')
    cli_parser.add_argument('--motifs', '-m', action='store_true')
    cli_parser.add_argument('--statistics', '-stats', action='store_true')
    cli_parser.add_argument('--planarity', '-pln', action='store_true')
    return cli_parser


def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments.
    """
    args = parser.parse_args()
    return args


def main():
    """
    Run the CLI.
    """

    start_time = time.time()

    # Create the standard config file and remember where it is stored.
    initialize.conf()
    conf_path = initialize.CONF_FILE

    # Remember the path from where the script is called.
    cwd = os.path.abspath(os.getcwd())

    # Parse arguments
    args = parse_args(create_parser())

    # Get the input and output directories. Create the output directory structure.
    in_dir = args.input[0]
    out_dir = args.output[0]
    out_subdirs = initialize.init_output_dir(out_dir)

    # If no arguments are given, run the complete pipeline.
    if not (args.pdb or args.models2chains or args.hydrogen or args.dssp or args.motifs or
            args.statistics or args.planarity or args.ppi is not None):

        # First get the input files (.pdbids). This can be multiple files.
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        # Get the PDB IDs from all .pdbids files.
        for pdb_list_file in input_files:
            if initialize.is_valid_input_file(os.path.join(in_dir, pdb_list_file)):
                pdb_ids += (utilities.get_pdb_ids_from_file(os.path.join(in_dir, pdb_list_file)))

        # Download the PDB files into the pdb subdirectory.
        os.chdir(out_subdirs['pdb'])
        pdb.pdb_download(pdb_ids, config.get_num_pdb_dl_threads(conf_path))

        # Copy the files into the mod_pdb directory. Later, for the PPI calculation, we need all
        # files in the mod_pdb subdir.
        # Since we have to change all PDB files (adding hydrogen atoms), we copy all PDB files over
        # into the mod_pdb subdir and change them inplace there.
        for files in os.listdir(out_subdirs['pdb']):
            shutil.copy(files, out_subdirs['mod_pdb'])

        # Convert models to chains and calculate hydrogen atoms. PDB files will be changed inplace.
        os.chdir(out_subdirs['mod_pdb'])
        for index in pdb_ids:
            plcc.pdb_models_to_chains(index, os.path.join(out_subdirs['mod_pdb'], index + '.pdb'))
            hydrogen.calc_hydrogen(index, os.path.join(out_subdirs['mod_pdb'], index + '.pdb'))

        # Get the list of all modified PDB files.
        pdb_files = []
        for pdb_file in os.listdir(out_subdirs['mod_pdb']):
            if initialize.is_pdb_file(os.path.join(out_subdirs['mod_pdb'], pdb_file)):
                pdb_files.append(pdb_file)

        # Download DSSP files into the dssp subdir.
        dssp.dssp_download(pdb_files, out_subdirs['dssp'],
                           config.get_num_dssp_dl_threads(conf_path))

        # Copy the DSSP files into the mod_pdb subdir. Again, we need those files for the PPI
        # calculations.
        for files in os.listdir(out_subdirs['dssp']):
            shutil.copy(os.path.join(out_subdirs['dssp'], files), out_subdirs['mod_pdb'])

        # Calculate PPIs for all PDB files in the mod_pdb subdir.
        for index in pdb_ids:
            plcc.calculate_ppi(index)

        # Calculate statistics based on the PPI results.
        num_pdb_files = 0
        edges_ppi = 0
        vertices_ppi = 0
        all_contacts = collections.Counter()
        all_aas = 0
        atom_contacts = {}
        for files in os.listdir(out_subdirs['mod_pdb']):
            if files.endswith('.fanmod'):
                edges_ppi += statistics.count_edges_in_ppi_aa_graph(
                    os.path.join(out_subdirs['mod_pdb'], files))
                vertices_ppi += statistics.count_vertices_in_ppi_aa_graph(
                    os.path.join(out_subdirs['mod_pdb'], files))
            elif files.endswith('.stats'):
                all_contacts += collections.Counter(statistics.count_all_contacts(
                    os.path.join(out_subdirs['mod_pdb'], files)))
            elif files.endswith('.csv'):
                atom_contacts.update(statistics.count_atom_num_contacts(
                    os.path.join(out_subdirs['mod_pdb'], files)))
            elif files.endswith('aagraph.gml'):
                all_aas += statistics.count_aas_in_aa_graph(os.path.join(out_subdirs['mod_pdb'],
                                                                         files))
            elif files.endswith('.pdb'):
                num_pdb_files += 1

        output_results.bar_chart_types_of_contacts(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_hydrogen(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_hydrogen_verbose(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_vdw(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_ligands(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_pi_effects(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_pi_effects_verbose(all_contacts, out_subdirs['imgs'])

        end_time = time.time()
        runtime = (end_time - start_time) / 60
        all_atom_atom_contacts = statistics.amount_atom_contacts(all_contacts)
        avg_num_edges_in_graph = edges_ppi / num_pdb_files
        avg_num_edges_on_atom_level = all_atom_atom_contacts / num_pdb_files
        aas_contributing = vertices_ppi
        avg_num_aas_per_graph = all_aas / num_pdb_files
        avg_num_aas_contributing_per_graph = aas_contributing / num_pdb_files
        avg_atom_atom_per_edge = avg_num_edges_on_atom_level / avg_num_edges_in_graph
        num_contacts_per_atom = sorted(atom_contacts.items(), key=lambda x: x[1], reverse=True)

        output_results.html_wrapper(out_subdirs['statistic'],
                                    time.asctime(time.localtime(end_time)),
                                    num_pdb_files,
                                    runtime,
                                    avg_num_edges_in_graph,
                                    avg_num_edges_on_atom_level,
                                    aas_contributing,
                                    avg_num_aas_per_graph,
                                    avg_num_aas_contributing_per_graph,
                                    avg_atom_atom_per_edge,
                                    all_aas,
                                    all_atom_atom_contacts,
                                    edges_ppi,
                                    num_contacts_per_atom)

    # If arguments are given, run those instead.

    # Download PDB files.
    if args.pdb:
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        for pdb_list_file in input_files:
            if initialize.is_valid_input_file(os.path.join(in_dir, pdb_list_file)):
                pdb_ids += (utilities.get_pdb_ids_from_file(os.path.join(in_dir, pdb_list_file)))

        os.chdir(out_subdirs['pdb'])
        pdb.pdb_download(pdb_ids, config.get_num_pdb_dl_threads(conf_path))
        os.chdir(cwd)

    # Convert Models in PDB files to chains.
    if args.models2chains:
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)
        for index in pdb_ids:
            plcc.pdb_models_to_chains(index, os.path.join(out_subdirs['mod_pdb'], index))
        os.chdir(cwd)

    # Add hydrogen atoms to the PDB file.
    if args.hydrogen:
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)
        for index in pdb_ids:
            hydrogen.calc_hydrogen(index, out_subdirs['mod_pdb'] + index)
        os.chdir(cwd)

    # Download DSSP files.
    if args.dssp:
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)
        dssp.dssp_download(pdb_ids, out_subdirs['dssp'], config.get_num_dssp_dl_threads(initialize
                                                                                        .CONF_FILE))
        os.chdir(cwd)

    if args.ppi == 'all':
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)
        os.chdir(in_dir)
        for index in pdb_ids:
            plcc.calculate_ppi(index)
        os.chdir(cwd)
    elif args.ppi == 'with-ligands':
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)
        os.chdir(in_dir)
        for index in pdb_ids:
            plcc.calculate_ppi_incl_ligands(index)
        os.chdir(cwd)
    elif args.ppi == 'no-pi-effects':
        raise NotImplementedError
    elif args.ppi == 'ligands-no-pi-effects':
        raise NotImplementedError

    if args.motifs:
        raise NotImplementedError
    if args.statistics:
        input_files = [f for f in os.listdir(in_dir)]
        num_pdb_files = 0
        edges_ppi = 0
        vertices_ppi = 0
        all_contacts = collections.Counter()
        all_aas = 0
        atom_contacts = {}
        for ppi_files in input_files:
            if ppi_files.endswith('.fanmod'):
                edges_ppi += statistics.count_edges_in_ppi_aa_graph(os.path.join(in_dir, ppi_files))
                vertices_ppi += statistics.count_vertices_in_ppi_aa_graph(os.path.join(in_dir,
                                                                                       ppi_files))
            elif ppi_files.endswith('.stats'):
                all_contacts += collections.Counter(statistics.count_all_contacts(
                    os.path.join(in_dir, ppi_files)))
            elif ppi_files.endswith('.csv'):
                atom_contacts.update(statistics.count_atom_num_contacts(os.path.join(in_dir,
                                                                                     ppi_files)))
            elif ppi_files.endswith('aagraph.gml'):
                all_aas += statistics.count_aas_in_aa_graph(os.path.join(in_dir, ppi_files))
            elif ppi_files.endswith('.pdb'):
                num_pdb_files += 1

        output_results.bar_chart_types_of_contacts(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_hydrogen(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_hydrogen_verbose(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_vdw(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_ligands(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_pi_effects(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_pi_effects_verbose(all_contacts, out_subdirs['imgs'])

        end_time = time.time()
        runtime = (end_time - start_time) / 60
        all_atom_atom_contacts = statistics.amount_atom_contacts(all_contacts)
        avg_num_edges_in_graph = edges_ppi / num_pdb_files
        avg_num_edges_on_atom_level = all_atom_atom_contacts / num_pdb_files
        aas_contributing = vertices_ppi
        avg_num_aas_per_graph = all_aas / num_pdb_files
        avg_num_aas_contributing_per_graph = aas_contributing / num_pdb_files
        avg_atom_atom_per_edge = avg_num_edges_on_atom_level / avg_num_edges_in_graph
        num_contacts_per_atom = sorted(atom_contacts.items(), key=lambda x: x[1], reverse=True)

        output_results.html_wrapper(out_subdirs['statistic'],
                                    time.asctime(time.localtime(end_time)),
                                    num_pdb_files,
                                    runtime,
                                    avg_num_edges_in_graph,
                                    avg_num_edges_on_atom_level,
                                    aas_contributing,
                                    avg_num_aas_per_graph,
                                    avg_num_aas_contributing_per_graph,
                                    avg_atom_atom_per_edge,
                                    all_aas,
                                    all_atom_atom_contacts,
                                    edges_ppi,
                                    num_contacts_per_atom)
    if args.planarity:
        raise NotImplementedError


if __name__ == '__main__':
    main()
