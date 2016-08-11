#!/usr/bin/env python

"""
ppicp.initialize
~~~~~~~~~~~~~~~~

Initialize folder structures and the config file.
Everything that is needed to run the actual application smoothly.
"""

import os
import sys

PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(PARENT_DIR)

from ppicp import config


CONF_FILE = os.path.expanduser('~') + '/.ppicp_settings'
BIN_DIR = os.path.abspath(os.path.join(os.path.join(os.path.dirname(__file__), os.path.pardir),
                                       'bin'))


def conf():
    """
    Write the config file if not already existing.
    """
    if not os.path.isfile(CONF_FILE):
        config.write_config(CONF_FILE)


def init_output_dir(out_dir):
    """
    Initialize the output directory structure.

    :param out_dir: Path to the output directory root.
    :return: a dictionary of the paths to the output subfolders.
    """
    pdb_files = os.path.join(out_dir, 'pdb_original')
    ppi_results = os.path.join(out_dir, 'ppi_results')
    motif_files = os.path.join(out_dir, 'motifs')
    statistic_files = os.path.join(out_dir, 'statistics')
    planarity_files = os.path.join(out_dir, 'planarity')

    check_output_path(pdb_files)
    check_output_path(ppi_results)
    check_output_path(motif_files)
    check_output_path(statistic_files)
    check_output_path(planarity_files)

    img_files = os.path.join(statistic_files, 'imgs')
    check_output_path(img_files)

    return {'pdb': pdb_files, 'ppi_results': ppi_results, 'motif': motif_files,
            'statistic': statistic_files, 'planarity': planarity_files, 'imgs': img_files}


def is_pdb_file(path):
    """
    Checks if the specified path points to a PDB file.

    :param path: to the input file.
    :return: True if the directory exists, false otherwise.
    """
    path = os.path.abspath(path)
    return bool(path.endswith('.pdb') and os.path.isfile(path))


def is_valid_input_file(path):
    """
    Checks if the specified path points to a valid file.

    :param path: to the input file.
    :return: True if the directory exists, false otherwise.
    """
    path = os.path.abspath(path)
    return bool(path.endswith('.pdbids') and os.path.isfile(path))


def is_valid_directory(path):
    """
    Checks if the specified path points to a valid directory.

    :param path: to the input directory.
    :return: True if directory exists, false otherwise.
    """
    path = os.path.abspath(path)
    return bool(os.path.isdir(path))


def check_input_path(path, parser):
    """
    Checks if the specified input is correct.
    That is, the input path points to a a file or a directory.

    :param path: input path.
    :param parser: the parser object.
    :return: the type of the input.
    """
    path = os.path.abspath(path)
    if is_valid_input_file(path):
        return path
    elif is_valid_directory(path):
        return path
    else:
        parser.error('The directory or file \"{}\" does not exist or is not valid!\nYou have to '
                     'either enter the path to a file with the file extension \".pdbids\" that '
                     'contains PDB IDs or you have to enter the path to a directory containing PDB '
                     'files.'.format(path))


def check_output_path(path):
    """
    Checks if the output path already exists, if not a new directory is created at the specified
    path.

    :param path: to the output directory
    :return: type of the output object.
    """
    path = os.path.abspath(path)
    if is_valid_directory(path):
        return path
    else:
        os.mkdir(path)
        print('Created new directory at \"{}\".'.format(path))
        return path
