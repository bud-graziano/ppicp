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

import config
import initialize


def calc_hydrogen(pdb_path, out_dir):
    """
    Calculates hydrogen atoms for a given PDB file and saves the added hydrogen atoms to a new PDB
    file. This is done by invoking the ``reduce`` software.

    :param pdb_path: The PDB file.
    :param out_dir: Where the resulting files are saved.
    :return: True if successful, False otherwise.
    """
    if platform.system() == 'Windows':
        raise NotImplementedError
    elif platform.system() == 'Linux':
        if config.get_hydrogen_app(initialize.CONF_FILE) == 'reduce':
            try:
                with open(out_dir + 'stripped', 'w') as out:
                    out.write(subprocess.check_output([os.path.join(initialize.BIN_DIR, 'reduce'),
                                                       '-Trim', '-Quiet', pdb_path]))

                with open(out_dir, 'w') as out:
                    subprocess.call([os.path.join(initialize.BIN_DIR, 'reduce'), '-build', '-DB',
                                     os.path.join(initialize.BIN_DIR, 'reduce_wwPDB_het_dict.txt'),
                                     '-Quiet', out_dir + 'stripped'], stdout=out)

                # Get rid of the overhead.
                os.remove(out_dir + 'stripped')
                return True
            except (OSError, subprocess.CalledProcessError) as err:
                print('{}\nHydrogen calculations failed.'.format(err))
                return False
        else:
            hydrogen_app = config.get_hydrogen_app(initialize.CONF_FILE)
            try:
                print subprocess.check_output([hydrogen_app])
                return True
            except (OSError, subprocess.CalledProcessError) as err:
                print('{}\nHydrogen calculations failed.'.format(err))
                return False
    else:
        print("[ERROR] OS could not be determined.")
        sys.exit(1)
