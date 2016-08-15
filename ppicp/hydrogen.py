#!/usr/bin/env python

"""
ppicp.hydrogen
~~~~~~~~~~~~~~

Functions that handle the stripping and re-adding/calculating of hydrogen atoms to PDB files.
"""

import os
import platform
import subprocess
import sys

PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(PARENT_DIR)

from ppicp import config
from ppicp import initialize


LOGGER = initialize.init_logger(__name__)
LOGGER_HYD = initialize.init_hyd_logger(__name__)


def calc_hydrogen(pdb_path, out_dir):
    """
    Calculates hydrogen atoms for a given PDB file and saves the added hydrogen atoms to a new PDB
    file. This is done by invoking the ``reduce`` software.

    :param pdb_path: The PDB file.
    :param out_dir: Where the resulting files are saved.
    :return: True if successful, False otherwise.
    """
    LOGGER.info("[HYDROGEN] Working on: %s", pdb_path)
    if platform.system() == 'Windows':
        raise NotImplementedError
    elif platform.system() == 'Linux':
        if config.get_hydrogen_app(initialize.CONF_FILE) == 'reduce':
            LOGGER.debug('Using Reduce to re-calculate hydrogen atoms.')
            try:
                LOGGER.debug('Stripping hydrogen atoms for %s', pdb_path)
                command = [os.path.join(initialize.BIN_DIR, 'reduce'), '-Trim', '-Quiet', pdb_path]
                subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                stdout, stderr = subp.communicate()

                with open(out_dir + 'stripped', 'w') as out:
                    out.write(stdout)

                if stderr != '':
                    LOGGER.debug(stderr)
                    LOGGER_HYD.error('{%s}\n %s', pdb_path, stderr)

                retcode = subp.poll()
                if retcode:
                    cmd = command
                    raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

                LOGGER.debug('Re-adding hydrogen atoms for %s', pdb_path)
                command = [os.path.join(initialize.BIN_DIR, 'reduce'), '-build', '-DB',
                           os.path.join(initialize.BIN_DIR, 'reduce_wwPDB_het_dict.txt'),
                           '-Quiet', out_dir + 'stripped']

                subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                stdout, stderr = subp.communicate()

                with open(out_dir, 'w') as out:
                    out.write(stdout)

                if stderr != '':
                    LOGGER.debug(stderr)
                    LOGGER_HYD.error('{%s}\n %s', pdb_path, stderr)

                # Get rid of the overhead.
                os.remove(out_dir + 'stripped')

                return True
            except (OSError, subprocess.CalledProcessError) as err:
                LOGGER.error('%s \nHydrogen calculations failed.', err)
                return False
        else:
            LOGGER.debug('Using user-specified application.')
            hydrogen_app = config.get_hydrogen_app(initialize.CONF_FILE)
            try:
                LOGGER.debug(subprocess.check_output([hydrogen_app]))
                return True
            except (OSError, subprocess.CalledProcessError) as err:
                LOGGER.error('%s \nHydrogen calculations failed.', err)
                return False
    else:
        LOGGER.critical("[ERROR] OS could not be determined.")
        sys.exit(1)
