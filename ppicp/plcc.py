#!/usr/bin/env python

"""
ppicp.plcc
~~~~~~~~~~

Calculate PPIs (in parallel and without ligands if specified) and convert models in PDB files to
chains.
"""

from __future__ import division

import multiprocessing
import multiprocessing.pool
import os
import platform
import subprocess
import sys

PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(PARENT_DIR)

from ppicp import config
from ppicp import initialize
from ppicp import utilities


LOGGER = initialize.init_logger(__name__)
LOGGER_PLCC = initialize.init_plcc_logger(__name__)


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
    Calculates the protein-protein interactions by executing the plcc.jar code with the necessary
    flags.

    :param pdb_path: PDB ID.
    :return: True if the calculation was successful, False otherwise.
    """
    LOGGER.info("[PPI] Working on: %s", pdb_path)
    if config.get_java_exec_path(initialize.CONF_FILE) is None:
        java_exec = utilities.get_jre_path()
    else:
        java_exec = config.get_java_exec_path(initialize.CONF_FILE)
    if platform.system() == 'Windows':
        try:
            command = [java_exec, '-jar', initialize.BIN_DIR + '/plcc.jar', pdb_path[:4],
                       '--alt-aa-contacts']
            subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = subp.communicate()
            LOGGER.debug(stdout)

            if stderr != '':
                LOGGER.debug(stderr)
                LOGGER_PLCC.error('{%s}\n %s', pdb_path, stderr)

            retcode = subp.poll()
            if retcode:
                cmd = command
                raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

            return True
        except (WindowsError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nPPI calculations failed.', err)
            return False
    elif platform.system() == 'Linux':
        try:
            command = [java_exec.rstrip('\n'), '-jar', os.path.join(initialize.BIN_DIR, 'plcc.jar'),
                       pdb_path[:4], '--alt-aa-contacts']
            subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = subp.communicate()
            LOGGER.debug(stdout)
            if stderr != '':
                LOGGER.debug(stderr)
                LOGGER_PLCC.error('{%s}\n %s', pdb_path, stderr)

            retcode = subp.poll()
            if retcode:
                cmd = command
                raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

            return True
        except (OSError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nPPI calculations failed.', err)
            return False
    else:
        LOGGER.critical("[ERROR] OS could not be determined.")
        sys.exit(1)


def calculate_ppi_incl_ligands(pdb_path):
    """
    Calculates the protein-protein interactions including ligands by executing the ``plcc.jar`` code
    with the necessary flags.

    :param pdb_path: PDB ID.
    :return: True if the calculation was successful, False otherwise.
    """
    LOGGER.info("[PPI] Working on: %s", pdb_path)
    if config.get_java_exec_path(initialize.CONF_FILE) is None:
        java_exec = utilities.get_jre_path()
    else:
        java_exec = config.get_java_exec_path(initialize.CONF_FILE)
    if platform.system() == 'Windows':
        try:
            command = [java_exec, '-jar', initialize.BIN_DIR + '/plcc.jar', pdb_path[:4],
                       '--alt-aa-contacts-ligands']
            subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = subp.communicate()
            LOGGER.debug(stdout)
            if stderr != '':
                LOGGER.debug(stderr)
                LOGGER_PLCC.error('{%s}\n %s', pdb_path, stderr)

            retcode = subp.poll()
            if retcode:
                cmd = command
                raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

            return True
        except (WindowsError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nPPI calculations failed.', err)
            return False
    elif platform.system() == 'Linux':
        try:
            command = [java_exec.rstrip('\n'), '-jar', os.path.join(initialize.BIN_DIR, 'plcc.jar'),
                       pdb_path[:4], '--alt-aa-contacts-ligands']

            subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = subp.communicate()
            LOGGER.debug(stdout)
            if stderr != '':
                LOGGER.debug(stderr)
                LOGGER_PLCC.error('{%s}\n %s', pdb_path, stderr)

            retcode = subp.poll()
            if retcode:
                cmd = command
                raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

            return True
        except (OSError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nPPI calculations failed.', err)
            return False
    else:
        LOGGER.critical("[ERROR] OS could not be determined.")
        sys.exit(1)


def pdb_models_to_chains(pdb_path, out_dir):
    """
    Runs the ``splitpdb.jar`` code with necessary flags that takes PDB files containing multiple
    modelsand generates PDB files that only contain one model. This is needed for later DSSP
    calculations to be correct.

    :param out_dir: Where the converted PDB files are saved.
    :param pdb_path: PDB ID.
    :return: True if the calculation was successful, False otherwise.
    """
    LOGGER.info("[MODELS2CHAINS] Working on: %s", pdb_path)
    if config.get_java_exec_path(initialize.CONF_FILE) is None:
        java_exec = utilities.get_jre_path()
    else:
        java_exec = config.get_java_exec_path(initialize.CONF_FILE)
    if platform.system() == 'Windows':
        try:
            command = [java_exec, '-jar', initialize.BIN_DIR + '/plcc.jar', pdb_path[:4],
                       '--convert-models-to-chains', out_dir]
            subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = subp.communicate()
            LOGGER.debug(stdout)
            if stderr != '':
                LOGGER.debug(stderr)
                LOGGER_PLCC.error('{%s}\n%s', pdb_path, stderr)

            retcode = subp.poll()
            if retcode:
                cmd = command
                raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

            return True
        except (WindowsError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nSplitting calculations failed.', err)
            return False
    elif platform.system() == 'Linux':
        try:
            command = [java_exec.rstrip('\n'), '-jar', os.path.join(initialize.BIN_DIR, 'plcc.jar'),
                       pdb_path[:4], '--convert-models-to-chains', out_dir]
            subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = subp.communicate()
            LOGGER.debug(stdout)
            if stderr != '':
                LOGGER.debug(stderr)
                LOGGER_PLCC.error('{%s}\n %s', pdb_path, stderr)

            retcode = subp.poll()
            if retcode:
                cmd = command
                raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

            return True
        except (OSError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nSplitting calculations failed.', err)
            return False
    else:
        LOGGER.critical("[ERROR] OS could not be determined.")
        sys.exit(1)
