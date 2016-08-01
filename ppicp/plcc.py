#!/usr/bin/env python

"""
Calculate PPIs (in parallel and without ligands if specified) and convert models in PDB files to
chains.
"""

from __future__ import division

import multiprocessing
import multiprocessing.pool
import platform
import subprocess
import sys

import config
import initialize
import utilities


class PtglWorker(multiprocessing.Process):
    """
    Worker class that makes it possible to spawn multiple processes executing the plcc.jar code.
    """
    def __init__(self, queue):
        super(PtglWorker, self).__init__()
        self.queue = queue

    def run(self):
        pdb_id = self.queue.get()
        calculate_ppi(pdb_id)
        self.queue.task_done()


def calculate_ppi(pdb_path):
    """
    Calculates the protein-protein interaction by executing the plcc.jar code with the necessary
    flags.
    :param pdb_path: PDB ID.
    :return: True if the calculation was successful, False otherwise.
    """
    print("[PPI] Working on: {}".format(pdb_path))
    if config.get_java_exec_path(initialize.CONF_FILE) is None:
        java_exec = utilities.get_jre_path()
    else:
        java_exec = config.get_java_exec_path(initialize.CONF_FILE)
    if platform.system() == 'Windows':
        try:
            print subprocess.check_output([java_exec, '-jar',
                                           initialize.BIN_DIR + '/plcc.jar',
                                           pdb_path.rstrip('.pdb'), '--alt-aa-contacts'])
            return True
        except WindowsError, (errno, strerror):
            print('WindowsError [Error {}]: {}\nPPI calculations failed.'.format(errno, strerror))
            return False
    elif platform.system() == 'Linux':
        try:
            print subprocess.check_output([java_exec, '-jar',
                                           initialize.BIN_DIR + '/plcc.jar',
                                           pdb_path.rstrip('.pdb'), '--alt-aa-contacts'])
            return True
        except OSError, (errno, strerror):
            print('OSError [Error {}]: {}\nPPI calculations failed.'.format(errno, strerror))
            return False
    else:
        print("[ERROR] OS could not be determined.")
        sys.exit(1)


def calculate_ppi_no_ligands(pdb_path):
    """
    Calculates the protein-protein interaction by executing the plcc.jar code with the necessary
    flags.
    :param pdb_path: PDB ID.
    :return: True if the calculation was successful, False otherwise.
    """
    print("[PPI] Working on: {}".format(pdb_path))
    if config.get_java_exec_path(initialize.CONF_FILE) is None:
        java_exec = utilities.get_jre_path()
    else:
        java_exec = config.get_java_exec_path(initialize.CONF_FILE)
    if platform.system() == 'Windows':
        try:
            print subprocess.check_output([java_exec, '-jar',
                                           initialize.BIN_DIR + '/plcc.jar',
                                           pdb_path.rstrip('.pdb'), '--alt-aa-contacts-no-ligands'])
            return True
        except WindowsError, (errno, strerror):
            print('WindowsError [Error {}]: {}\nPPI calculations failed.'.format(errno, strerror))
            return False
    elif platform.system() == 'Linux':
        try:
            print subprocess.check_output([java_exec, '-jar',
                                           initialize.BIN_DIR + '/plcc.jar',
                                           pdb_path.rstrip('.pdb'), '--alt-aa-contacts-no-ligands'])
            return True
        except OSError, (errno, strerror):
            print('OSError [Error {}]: {}\nPPI calculations failed.'.format(errno, strerror))
            return False
    else:
        print("[ERROR] OS could not be determined.")
        sys.exit(1)


def pdb_models_to_chains(pdb_path, out_dir):
    """
    Runs the splitpdb.jar code with necessary flags that takes PDB files containing multiple models
    and generates PDB files that only contain one model. This is needed for later DSSP calculations
    to be correct.
    :param out_dir: Where the converted PDB files are saved.
    :param pdb_path: PDB ID.
    :return: True if the calculation was successful, False otherwise.
    """
    print("[MODELS2CHAINS] Working on: {}".format(pdb_path))
    if config.get_java_exec_path(initialize.CONF_FILE) is None:
        java_exec = utilities.get_jre_path()
    else:
        java_exec = config.get_java_exec_path(initialize.CONF_FILE)
    if platform.system() == 'Windows':
        try:
            print subprocess.check_output([java_exec, '-jar',
                                           initialize.BIN_DIR + '/plcc.jar',
                                           pdb_path.rstrip('.pdb'),
                                           '--convert-models-to-chains',
                                           pdb_path, out_dir + '.split'])
            return True
        except WindowsError, (errno, strerror):
            print('WindowsError [Error {}]: {}\nSplitting calculations failed.'.format(errno,
                                                                                       strerror))
            return False
    elif platform.system() == 'Linux':
        try:
            print subprocess.check_output([java_exec, '-jar',
                                           initialize.BIN_DIR + '/plcc.jar',
                                           pdb_path.rstrip('.pdb'),
                                           '--convert-models-to-chains',
                                           pdb_path, out_dir + '.split'])
            return True
        except OSError, (errno, strerror):
            print('OSError [Error {}]: {}\nSplitting calculations failed.'.format(errno, strerror))
            return False
    else:
        print("[ERROR] OS could not be determined.")
        sys.exit(1)
