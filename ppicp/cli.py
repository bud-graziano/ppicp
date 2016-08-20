#!/usr/bin/env python

"""
PPICP -- Protein-Protein Interaction Calculation Pipeline
---------------------------------------------------------

A scientific software for the calculation of protein-protein interactions (PPIs).
The pipeline runs various third party software to compute PPIs.
It consists of four major steps:

    1. Downloading of input data based on user input (data from PDB, DSSP, and Reduce).
    2. Calculation of PPIs.
    3. Motif detection in PPIs-graphs.
    4. Compiling of statistics based on the previous output.

For more information about the third party software and the pipeline itself please refer to the
documentation.

ppicp.cli
~~~~~~~~~
"""

from __future__ import division


import argparse
import collections
import os
import shutil
import time
import sys

import clint.textui as clt

PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(PARENT_DIR)

from ppicp import config
from ppicp import dssp
from ppicp import hydrogen
from ppicp import initialize
from ppicp import motifs
from ppicp import output_results
from ppicp import rcsb_pdb
from ppicp import plcc
from ppicp import statistics
from ppicp import utilities


def create_parser():
    """
    Create a CLI parser.

    :return: the parser object.
    """
    cli_parser = argparse.ArgumentParser(description="PPICP -- Protein-Protein Interaction "
                                                     "Calculation Pipeline",
                                         epilog="PPICP is a scientific software to calculate "
                                                "protein-protein interactions from RCSB PDB data. "
                                                "The pipeline includes the execution of the "
                                                "following steps:\n1. Downloading of input data "
                                                "based on a list of PDB-IDs.\n    a. Data "
                                                "from RCSB PDB.\n    b. Data from DSSP.\n2. "
                                                "Converting models in PDB files to separate chains."
                                                "\n3. Removing existing hydrogen atoms and "
                                                "re-adding them with the Reduce software.\n4. "
                                                "Calculation of protein-protein interactions based "
                                                "on graph-theoretical methods.\n5. Detection of "
                                                "motifs in the calculated PPI-graphs with the "
                                                "Fanmod software.\n6. Compiling of statistics "
                                                "based on the PPI calculation output.",
                                         formatter_class=argparse.RawDescriptionHelpFormatter)

    cli_parser.add_argument('input', type=lambda x: initialize.check_input_path(x, cli_parser),
                            nargs=1, help="Path to a directory containing a '.pdbids' file that "
                                          "stores PDB-IDs separated by new lines.")

    cli_parser.add_argument('output', type=lambda x: initialize.check_output_path(x), nargs=1,
                            help="Path to an output directory.")
    cli_parser.add_argument('--ppi', '-p',
                            choices=['standard', 'with-ligands', 'no-pi-effects',
                                     'ligands-no-pi-effects'],
                            help="Calculate protein-protein interactions. Standard runs the PPI "
                                 "calculation including pi-effects but without ligand contacts."
                                 "'with-ligands' includes ligands in the PPI calculation. 'no-pi-"
                                 "effects' excludes pi-effects from the PPI calculation. 'ligands-"
                                 "no-pi-effects' includes ligands and excludes pi-effects. This "
                                 "requires PDB and DSSP files.")
    cli_parser.add_argument('--pdb', '-pdb', action='store_true',
                            help="Download PDB files based on PDB-IDs from a '.pdbids' file.")
    cli_parser.add_argument('--dssp', '-d', action='store_true',
                            help="Download DSSP files based on PDB files.")
    cli_parser.add_argument('--hydrogen', '-r', action='store_true',
                            help="Re-calculate hydrogen atoms for given PDB files.")
    cli_parser.add_argument('--models2chains', '-s', action='store_true',
                            help="Convert models in PDB files to separate chains.")
    cli_parser.add_argument('--motifs', '-m', action='store_true',
                            help="Detect motifs of size 3-8 in the generated PPI-graphs. This "
                                 "requires PPI calculations to be "
                                 "already finished.")
    cli_parser.add_argument('--statistics', '-stats', action='store_true',
                            help="Compile statistics based on PPI and motif calculations.")
    cli_parser.add_argument('--planarity', '-pln', action='store_true',
                            help="Calculate the planarity for aromatic rings in PDB files.")
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

    initialize.conf()

    start_time = time.time()
    print (clt.colored.green('Starting the pipeline. {} \n'.format(time.asctime(
        time.localtime(start_time))), bold=True))

    logger = initialize.init_logger(__name__)
    logger.info('Starting the pipeline.')

    # Create the standard config file and remember where it is stored.
    logger.debug('Initialization.')
    conf_path = initialize.CONF_FILE

    # Remember the path from where the script is called.
    cwd = os.path.abspath(os.getcwd())

    # Parse arguments
    logger.debug('Parsing arguments.')
    args = parse_args(create_parser())
    logger.debug('Supplied arguments are: %s', args)

    # Get the input and output directories.
    in_dir = args.input[0]
    logger.debug('Input directory set as: %s', in_dir)

    out_dir = args.output[0]
    logger.debug('Output directory set as: %s', out_dir)

    # If no arguments are given, run the complete pipeline.
    if not (args.pdb or args.models2chains or args.hydrogen or args.dssp or args.motifs or
            args.statistics or args.planarity or args.ppi is not None):

        logger.info('Starting complete run through the pipeline.')

        # Create the output directory structure.
        out_subdirs = initialize.init_output_dir(out_dir)
        logger.debug('Creating sub-directory structure. %s', out_subdirs)

        # First get the input files (.pdbids). This can be multiple files.
        logger.debug('Collecting PDB-IDs.')
        print (clt.colored.cyan('\nPDB'))

        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []

        # Get the PDB IDs from all .pdbids files.
        print 'Collecting PDB-IDs'
        for pdb_list_file in input_files:
            if initialize.is_valid_input_file(os.path.join(in_dir, pdb_list_file)):
                pdb_ids += (utilities.get_pdb_ids_from_file(os.path.join(in_dir, pdb_list_file)))

        # Download the PDB files into the pdb subdirectory.
        os.chdir(out_subdirs['pdb'])
        logger.info('Start downloading PDB files.')
        print 'Start downloading PDB files'

        rcsb_pdb.pdb_download(pdb_ids, config.get_num_pdb_dl_threads(conf_path))

        # Copy the files into the mod_pdb directory. Later, for the PPI calculation, we need all
        # files in the mod_pdb subdir.
        # Since we have to change all PDB files (adding hydrogen atoms), we copy all PDB files over
        # into the mod_pdb subdir and change them inplace there.
        logger.debug('Copying PDB files into the %s sub-directory.', (out_subdirs['pdb']))

        for files in os.listdir(out_subdirs['pdb']):
            shutil.copy(files, out_subdirs['ppi_results'])

        # Convert models to chains and calculate hydrogen atoms. PDB files will be changed inplace.
        print (clt.colored.cyan('\nModels to Chains and Hydrogen Atoms'))

        os.chdir(out_subdirs['ppi_results'])

        logger.info('Start converting models to chains and re-calculating hydrogen atoms.')
        print 'Start converting models to chains and re-calculating hydrogen atoms.'
        for index in clt.progress.bar(pdb_ids):
            plcc.pdb_models_to_chains(index.lower(), os.path.join(out_subdirs['ppi_results'],
                                                                  index.lower() + '.pdb'))
            hydrogen.calc_hydrogen(index.lower() + '.pdb', os.path.join(out_subdirs['ppi_results'],
                                                                        index.lower() + '.pdb'))

        # Get the list of all PDB files including the modified.
        logger.debug('Get a list of all PDB files in %s', (out_subdirs['ppi_results']))
        pdb_files = []
        for pdb_file in os.listdir(out_subdirs['ppi_results']):
            if initialize.is_pdb_file(os.path.join(out_subdirs['ppi_results'], pdb_file)):
                pdb_files.append(pdb_file)

        # Download DSSP files.
        print (clt.colored.cyan('\nDSSP'))
        logger.info('Start downloading DSSP files.')
        print 'Start downloading DSSP files'
        dssp.dssp_download(pdb_files, out_subdirs['ppi_results'],
                           config.get_num_dssp_dl_threads(conf_path))

        # Calculate PPIs for all PDB files in the mod_pdb subdir.
        print (clt.colored.cyan('\nPPI'))
        logger.info('Start PPI calculations (without ligands).')
        print 'Start PPI calculations (without ligands)'
        for index in clt.progress.bar(pdb_ids):
            plcc.calculate_ppi(index.lower())

        # Combine the calculated single .fanmod files into one network.
        # (out_dir, out_dir) as input here looks strange but it's correct since the files we need
        # to combine are in the output directory already and the file that stores the combined ones
        # should also be in the output directory.
        motifs.combine_fanmod_files(out_subdirs['ppi_results'], out_subdirs['ppi_results'])
        os.chdir(cwd)

        # Detect motifs.
        print (clt.colored.cyan('\nMotifs'))
        logger.info('Start motif detection.')
        print 'Start motif detection'
        for fm_file in clt.progress.bar(os.listdir(out_subdirs['ppi_results'])):
            if fm_file.endswith('.fanmod'):
                for motif_size in xrange(3, 9):
                    logger.debug('Detect motifs of size %d', motif_size)
                    motifs.calculate_motifs(motif_size,
                                            os.path.join(out_subdirs['ppi_results'], fm_file),
                                            os.path.join(out_subdirs['motif'], fm_file[:4]) +
                                            '_{}_fanmod'.format(motif_size))

        # Calculate statistics based on the PPI results.
        print (clt.colored.cyan('\nStatistics'))

        num_pdb_files = 0
        edges_ppi = 0
        vertices_ppi = 0
        all_contacts = collections.Counter()
        all_aas = 0
        atom_contacts = {}
        aa_chem_props_3_ppi = collections.Counter()
        aa_chem_props_5_ppi = collections.Counter()
        aa_chem_props_3_all = collections.Counter()
        aa_chem_props_5_all = collections.Counter()

        logger.info('Start compiling statistics.')
        print 'Start compiling statistics'
        for files in clt.progress.bar(os.listdir(out_subdirs['ppi_results'])):
            if files.endswith('.fanmod'):
                edges_ppi += statistics.count_edges_in_ppi_aa_graph(
                    os.path.join(out_subdirs['ppi_results'], files))
                vertices_ppi += statistics.count_vertices_in_ppi_aa_graph(
                    os.path.join(out_subdirs['ppi_results'], files))
            elif files.endswith('.stats'):
                all_contacts += collections.Counter(statistics.count_all_contacts(
                    os.path.join(out_subdirs['ppi_results'], files)))
            elif files.endswith('.csv'):
                atom_contacts.update(statistics.count_atom_num_contacts(
                    os.path.join(out_subdirs['ppi_results'], files)))
            elif files.endswith('aagraph.gml'):
                all_aas += statistics.count_aas_in_aa_graph(os.path.join(out_subdirs['ppi_results'],
                                                                         files))
                chem_props_3, chem_props_5 = statistics.get_chem_props(
                    os.path.join(out_subdirs['ppi_results'], files))
                aa_chem_props_3_all += collections.Counter(chem_props_3)
                aa_chem_props_5_all += collections.Counter(chem_props_5)
            elif files.endswith('.id'):
                chem_props_3, chem_props_5 = statistics.get_chem_props(
                    os.path.join(out_subdirs['ppi_results'], files))
                aa_chem_props_3_ppi += collections.Counter(chem_props_3)
                aa_chem_props_5_ppi += collections.Counter(chem_props_5)
            elif files.endswith('.pdb'):
                num_pdb_files += 1

        # Create and save the charts.
        logger.debug('Saving charts to %s.', (out_subdirs['imgs']))
        output_results.bar_chart_types_of_contacts(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_hydrogen(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_hydrogen_verbose(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_vdw(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_ligands(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_pi_effects(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_pi_effects_verbose(all_contacts, out_subdirs['imgs'])
        output_results.bar_chart_chem_props_3(aa_chem_props_3_ppi, out_subdirs['imgs'])
        output_results.bar_chart_chem_props_5(aa_chem_props_5_ppi, out_subdirs['imgs'])
        output_results.bar_chart_chem_props_3(aa_chem_props_3_all, out_subdirs['imgs'], False)
        output_results.bar_chart_chem_props_5(aa_chem_props_5_all, out_subdirs['imgs'], False)

        # Calculate statistics.
        logger.debug('Calculating statistics.')
        end_time = time.time()
        runtime = (end_time - start_time) / 60
        all_atom_atom_contacts = statistics.amount_atom_contacts(all_contacts)
        aas_contributing = vertices_ppi

        avg_num_edges_in_graph = -1
        avg_num_edges_on_atom_level = -1
        avg_num_aas_per_graph = -1
        avg_num_aas_contributing_per_graph = -1
        avg_atom_atom_per_edge = -1

        try:
            avg_num_edges_in_graph = edges_ppi / num_pdb_files
            avg_num_edges_on_atom_level = all_atom_atom_contacts / num_pdb_files
            avg_num_aas_per_graph = all_aas / num_pdb_files
            avg_num_aas_contributing_per_graph = aas_contributing / num_pdb_files
            avg_atom_atom_per_edge = avg_num_edges_on_atom_level / avg_num_edges_in_graph
        except ZeroDivisionError as err:
            logger.warning('%s \nEither the number of PDB files was zero or the average number of '
                           'edges in a PPI AA graph was zero.', err)

            if avg_num_edges_in_graph == 0:
                avg_atom_atom_per_edge = -1
            else:
                avg_num_edges_in_graph = -1
                avg_num_edges_on_atom_level = -1
                avg_num_aas_per_graph = -1
                avg_num_aas_contributing_per_graph = -1
                avg_atom_atom_per_edge = -1

        num_contacts_per_atom = sorted(atom_contacts.items(), key=lambda x: x[1], reverse=True)

        # Save everything in the output HTML file.
        logger.info('Save statistics results as HTML in %s', out_dir)
        output_results.html_wrapper(out_subdirs['statistic'],
                                    time.asctime(time.localtime(end_time)),
                                    num_pdb_files,
                                    format(runtime, '.4f'),
                                    format(avg_num_edges_in_graph, '.2f'),
                                    format(avg_num_edges_on_atom_level, '.2f'),
                                    aas_contributing,
                                    format(avg_num_aas_per_graph, '.2f'),
                                    format(avg_num_aas_contributing_per_graph, '.2f'),
                                    format(avg_atom_atom_per_edge, '.2f'),
                                    all_aas,
                                    all_atom_atom_contacts,
                                    edges_ppi,
                                    num_contacts_per_atom)

        print (clt.colored.green('\nClosing the pipeline. {} \n'.format(time.asctime(
            time.localtime(end_time))), bold=True))
        logger.info('Closing the pipeline.')
        sys.exit(0)

    # If arguments are given, run those instead.

    out_dir = initialize.check_output_path(out_dir)
    # Download PDB files.
    if args.pdb:
        logger.debug('Collecting PDB-IDs.')
        print (clt.colored.cyan('\nPDB'))

        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []

        # Get the PDB IDs from all .pdbids files.
        print 'Collecting PDB-IDs'
        for pdb_list_file in input_files:
            if initialize.is_valid_input_file(os.path.join(in_dir, pdb_list_file)):
                pdb_ids += (utilities.get_pdb_ids_from_file(os.path.join(in_dir, pdb_list_file)))

        os.chdir(out_dir)

        # Download the PDB files.
        logger.info('Start downloading PDB files.')
        print 'Start downloading PDB files'

        rcsb_pdb.pdb_download(pdb_ids, config.get_num_pdb_dl_threads(conf_path))
        os.chdir(cwd)

    # Convert Models in PDB files to chains.
    if args.models2chains:
        print (clt.colored.cyan('\nModels to Chains'))

        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []

        logger.debug('Collect input files (PDB) for models to chains conversion.')
        print 'Collecting input files (PDB)'
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)

        logger.info('Start converting models to chains.')
        print 'Start converting models to chains'
        for index in clt.progress.bar(pdb_ids):
            plcc.pdb_models_to_chains(index, os.path.join(out_dir, index))
        os.chdir(cwd)

    # Add hydrogen atoms to the PDB file.
    if args.hydrogen:
        print (clt.colored.cyan('\nHydrogen Atoms'))

        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []

        logger.debug('Collect input files (PDB) for hydrogen calculations.')
        print 'Collecting input files (PDB)'
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)

        logger.info('Start re-calculating hydrogen atoms.')
        print 'Start re-calculating hydrogen atoms'
        for index in clt.progress.bar(pdb_ids):
            hydrogen.calc_hydrogen(index, os.path.join(out_dir, index))
        os.chdir(cwd)

    # Download DSSP files.
    if args.dssp:
        print (clt.colored.cyan('\nDSSP'))

        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []

        logger.debug('Collect input files (PDB) for DSSP calculations.')
        print 'Collecting input files (PDB)'
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)

        logger.info('Start downloading DSSP files.')
        print 'Start downloading DSSP files'
        dssp.dssp_download(pdb_ids, out_dir, config.get_num_dssp_dl_threads(initialize.CONF_FILE))
        os.chdir(cwd)

    if args.ppi == 'standard':
        print (clt.colored.cyan('\nPPI'))

        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []

        logger.debug('Collect input files (PDB, DSSP) for PPI calculations.')
        print 'Collecting input files (PDB, DSSP)'
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)

        logger.info('Start PPI calculations (without ligands).')
        print 'Start PPI calculations (without ligands)'
        for index in clt.progress.bar(pdb_ids):
            plcc.calculate_ppi(index)
        for files in os.listdir(in_dir):
            if files.endswith(('.gml', '.stats', '.fanmod', '.id', '.csv', '.py')) \
                    and in_dir != out_dir:
                shutil.copy(files, out_dir)
                os.remove(files)

        # (out_dir, out_dir) as input here looks strange but it's correct since the files we need
        # to combine are in the output directory already and the file that stores the combined ones
        # should also be in the output directory.
        motifs.combine_fanmod_files(out_dir, out_dir)
        os.chdir(cwd)

    elif args.ppi == 'with-ligands':
        print (clt.colored.cyan('\nPPI'))

        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []

        logger.debug('Collect input files (PDB, DSSP) for PPI calculations (with ligands).')
        print 'Collecting input files (PDB, DSSP)'
        for pdb_id in input_files:
            if initialize.is_pdb_file(os.path.join(in_dir, pdb_id)):
                pdb_ids.append(pdb_id)

        os.chdir(in_dir)

        logger.info('Start PPI calculations (with ligands).')
        print 'Start PPI calculations (with ligands)'
        for index in clt.progress.bar(pdb_ids):
            plcc.calculate_ppi_incl_ligands(index)
        for files in os.listdir(in_dir):
            if files.endswith(('.gml', '.stats', '.fanmod', '.id', '.csv', '.py')) \
                    and in_dir != out_dir:
                shutil.move(files, out_dir)
        os.chdir(cwd)
    elif args.ppi == 'no-pi-effects':
        logger.info('Start PPI calculations (no pi-effects).')
        raise NotImplementedError
    elif args.ppi == 'ligands-no-pi-effects':
        logger.info('Start PPI calculations (with ligands, no pi-effects).')
        raise NotImplementedError

    # Detect motifs.
    if args.motifs:
        print (clt.colored.cyan('\nMotifs'))

        input_files = [f for f in os.listdir(in_dir)]

        logger.info('Start motif detection.')
        print 'Start motif detection'
        for fm_file in clt.progress.bar(input_files):
            if fm_file.endswith('.fanmod'):
                for motif_size in xrange(3, 9):
                    logger.debug('Detect motifs of size %d', motif_size)
                    motifs.calculate_motifs(motif_size,
                                            os.path.join(in_dir, fm_file),
                                            os.path.join(out_dir, fm_file[:4]) +
                                            '_{}_fanmod'.format(motif_size))

    if args.statistics:
        print (clt.colored.cyan('\nStatistics'))

        input_files = [f for f in os.listdir(in_dir)]
        num_pdb_files = 0
        edges_ppi = 0
        vertices_ppi = 0
        all_contacts = collections.Counter()
        all_aas = 0
        atom_contacts = {}
        aa_chem_props_3_ppi = collections.Counter()
        aa_chem_props_5_ppi = collections.Counter()
        aa_chem_props_3_all = collections.Counter()
        aa_chem_props_5_all = collections.Counter()

        logger.info('Start compiling statistics.')
        print 'Start compiling statistics'
        for ppi_files in clt.progress.bar(input_files):
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
                chem_props_3, chem_props_5 = statistics.get_chem_props(os.path.join(in_dir,
                                                                                    ppi_files))
                aa_chem_props_3_all += collections.Counter(chem_props_3)
                aa_chem_props_5_all += collections.Counter(chem_props_5)
            elif ppi_files.endswith('.id'):
                chem_props_3, chem_props_5 = statistics.get_chem_props(os.path.join(in_dir,
                                                                                    ppi_files))
                aa_chem_props_3_ppi += collections.Counter(chem_props_3)
                aa_chem_props_5_ppi += collections.Counter(chem_props_5)
            elif ppi_files.endswith('.pdb'):
                num_pdb_files += 1
        # Create the 'imgs' subdirectory if it is not already there.
        initialize.check_output_path(os.path.join(out_dir, 'imgs'))

        # Create and save the charts.
        logger.debug('Saving charts to %s.', (os.path.join(out_dir, 'imgs')))
        output_results.bar_chart_types_of_contacts(all_contacts, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_hydrogen(all_contacts, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_hydrogen_verbose(all_contacts, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_vdw(all_contacts, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_ligands(all_contacts, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_pi_effects(all_contacts, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_pi_effects_verbose(all_contacts, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_chem_props_3(aa_chem_props_3_ppi, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_chem_props_5(aa_chem_props_5_ppi, os.path.join(out_dir, 'imgs'))
        output_results.bar_chart_chem_props_3(aa_chem_props_3_all, os.path.join(out_dir, 'imgs'),
                                              False)
        output_results.bar_chart_chem_props_5(aa_chem_props_5_all, os.path.join(out_dir, 'imgs'),
                                              False)

        # Calculate statistics.
        logger.debug('Calculating statistics.')
        end_time = time.time()
        runtime = (end_time - start_time) / 60
        all_atom_atom_contacts = statistics.amount_atom_contacts(all_contacts)
        aas_contributing = vertices_ppi

        avg_num_edges_in_graph = -1
        avg_num_edges_on_atom_level = -1
        avg_num_aas_per_graph = -1
        avg_num_aas_contributing_per_graph = -1
        avg_atom_atom_per_edge = -1

        try:
            avg_num_edges_in_graph = edges_ppi / num_pdb_files
            avg_num_edges_on_atom_level = all_atom_atom_contacts / num_pdb_files
            avg_num_aas_per_graph = all_aas / num_pdb_files
            avg_num_aas_contributing_per_graph = aas_contributing / num_pdb_files
            avg_atom_atom_per_edge = avg_num_edges_on_atom_level / avg_num_edges_in_graph
        except ZeroDivisionError as err:
            logger.warning('%s \nEither the number of PDB files was zero or the average number of '
                           'edges in a PPI AA graph was zero.', err)
            if avg_num_edges_in_graph == 0:
                avg_atom_atom_per_edge = -1
            else:
                avg_num_edges_in_graph = -1
                avg_num_edges_on_atom_level = -1
                avg_num_aas_per_graph = -1
                avg_num_aas_contributing_per_graph = -1
                avg_atom_atom_per_edge = -1

        num_contacts_per_atom = sorted(atom_contacts.items(), key=lambda x: x[1], reverse=True)

        # Save everything in the output HTML file.
        logger.info('Save statistics results as HTML in %s.', out_dir)
        output_results.html_wrapper(out_dir,
                                    time.asctime(time.localtime(end_time)),
                                    num_pdb_files,
                                    format(runtime, '.4f'),
                                    format(avg_num_edges_in_graph, '.2f'),
                                    format(avg_num_edges_on_atom_level, '.2f'),
                                    aas_contributing,
                                    format(avg_num_aas_per_graph, '.2f'),
                                    format(avg_num_aas_contributing_per_graph, '.2f'),
                                    format(avg_atom_atom_per_edge, '.2f'),
                                    all_aas,
                                    all_atom_atom_contacts,
                                    edges_ppi,
                                    num_contacts_per_atom)
    if args.planarity:
        logger.info('Start planarity calculations.')
        raise NotImplementedError

    end_time = time.time()
    print (clt.colored.green('\nClosing the pipeline. {} \n'.format(time.asctime(
        time.localtime(end_time))), bold=True))
    logger.info('Closing the pipeline.')
    sys.exit(0)

if __name__ == '__main__':
    main()
