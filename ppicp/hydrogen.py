#!/usr/bin/env python

"""
Functions that handle the stripping and re-adding/calculating of hydrogen atoms to PDB files.
"""

import os
import platform
import subprocess
import sys

import config
import initialize


def calc_hydrogen(pdb_path, out_dir):
    """
    Calculates hydrogen atoms for a given PDB file and saves the added hydrogen atoms to a new PDB
    file.
    :param pdb_path: The PDB file.
    :param out_dir: Where the resulting files are saved.
    :return: True if successful, False otherwise.
    """
    if platform.system() == 'Windows':
        raise NotImplementedError
    elif platform.system() == 'Linux':
        if config.get_hydrogen_app(initialize.CONF_FILE) is 'reduce':
            try:
                print subprocess.check_output([os.path.join(initialize.BIN_DIR, 'reduce'),
                                               '-Trim', '-Quiet',
                                               pdb_path + ' > ' + out_dir + '.stripped_h'])
                print subprocess.check_output([os.path.join(initialize.BIN_DIR, 'reduce'),
                                               '-build', '-DB',
                                               os.path.join(initialize.BIN_DIR,
                                                            'reduce_wwPDB_het_dict.txt'), '-Quiet',
                                               out_dir + '.stripped_h' + ' > ' + out_dir])
                return True
            except (OSError, subprocess.CalledProcessError) as e:
                print('[Error {}]: {}\nHydrogen calculations failed.'.format(e.returncode,
                                                                             e.output))
                return False
        else:
            hydrogen_app = config.get_hydrogen_app(initialize.CONF_FILE)
            try:
                print subprocess.check_output([hydrogen_app])
                return True
            except (OSError, subprocess.CalledProcessError) as e:
                print('OSError [Error {}]: {}\nHydrogen calculations failed.'.format(e.returncode,
                                                                                     e.output))
                return False
    else:
        print("[ERROR] OS could not be determined.")
        sys.exit(1)
