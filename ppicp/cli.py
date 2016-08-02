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

import argparse
import os

import config
import dssp
import initialize
import pdb
import plcc
import hydrogen
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
    :return:
    """
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
        print('Run whole pipeline.')
        raise NotImplementedError

    # If arguments are given, run those instead.

    # Download PDB files.
    if args.pdb:
        input_files = [f for f in os.listdir(in_dir)]
        pdb_ids = []
        for pdb_list_file in input_files:
            if initialize.is_valid_input_file(os.path.join(in_dir, pdb_list_file)):
                pdb_ids += (utilities.get_pdb_ids_from_file(os.path.join(in_dir, pdb_list_file)))
        print(pdb_ids)
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
            plcc.pdb_models_to_chains(index, out_subdirs['mod_pdb'] + index)
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
        pass
    if args.planarity:
        raise NotImplementedError

if __name__ == '__main__':
    main()
