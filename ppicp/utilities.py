#!/usr/bin/env python

"""
Toolkit of various functions.
"""

import errno
import os
import platform
import subprocess
import string
if platform.system() == 'Windows':
    from ctypes import windll


def get_pdb_ids_from_file(path_to_file):
    """
    Takes a path to a file that stores line-separated PDB-IDs and returns a list of all found
    PDB-IDs.
    :param path_to_file: path to a file containing line separated PDB-IDs.
    :return: List of all found PDB-IDs.
    """
    with open(path_to_file, 'r') as f:
        content = f.read().splitlines()
        pdb_ids = []
        for line in content:
            pdb_ids.append(line)
    return pdb_ids


def get_drives():
    """
    Get all the available drives on the computer.
    :return: the drives.
    """
    drives = []
    bitmask = windll.kernel32.GetLogicalDrives()
    for letter in string.uppercase:
        if bitmask & 1:
            drives.append(letter)
        bitmask >>= 1
    return drives


def get_jre_path():
    """
    Get the path to the java runtime environment.
    :return: the path to the jre.
    """
    if platform.system() == 'Windows':
        # jre_regex = re.compile('[C-Z]:\\\\Program Files\\\\Java\\\\jdk1\.7\.(?!.*?jre)(?=.*?bin)',
        #  re.IGNORECASE)
        drives = get_drives()
        for drive in drives:
            paths = os.walk(drive + ':\\')
            for path in paths:
                if 'java.exe' in path[2]:
                    return os.path.join(path[0], 'java.exe')

    elif platform.system() == 'Linux':
        jre = subprocess.check_output(['which', 'java'])
        if jre is not None:
            return jre
    raise OSError(errno.ENOENT, os.strerror(errno.ENOENT, 'java.exe'))
