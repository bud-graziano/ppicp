#!/usr/bin/env python

"""
Handles the config of the PPI Calculation Pipeline.
"""

import ConfigParser
import datetime


def write_config(out_path):
    """
    Write the config file with default parameters to the user's home directory if it not already
    exists there.
    :param out_path: Where the config file should be stored.
    """
    config = ConfigParser.SafeConfigParser()

    config.add_section('java')
    config.set('java', '#executable_path', 'set/abspath/to/java/executable')
    config.add_section('io')
    config.set('io', 'num_pdb_download_threads', '2')
    config.set('io', 'num_dssp_download_threads', '5')
    config.add_section('external')
    config.set('external', 'calc_hydrogens', 'reduce')
    config.set('external', 'calc_motifs', 'fanmod')

    with open(out_path, 'w') as f:
        f.write('#These are the settings for PPICP\n#For further details on the values please '
                'refer to the documentation.\n#Lines starting with # are comments and will '
                'therefore be ignored.\n#{}\n\n'.format(datetime.datetime.now()
                                                        .strftime('%a %b %Y %H:%M:%S')))
        config.write(f)


def get_num_pdb_dl_threads(conf_path):
    """
    Get the number of threads from the config file that should be used to download PDB files.
    :param conf_path: Path to the config file.
    :return: number of threads.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.getint('io', 'num_pdb_download_threads')


def get_num_dssp_dl_threads(conf_path):
    """
    Get the number of threads from the config file that should be used to download DSSP files.
    :param conf_path: Path to the config file.
    :return: number of threads.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.getint('io', 'num_dssp_download_threads')


def get_java_exec_path(conf_path):
    """
    Get the path to the java executable. If no path is found, application will later try to find
    the path on its own.
    :param conf_path: Path to the config file.
    :return: path to the java executable.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    try:
        return config.get('java', 'executable_path')
    except ConfigParser.NoOptionError:
        pass


def get_hydrogen_app(conf_path):
    """
    Get the command-line string that will be run in order to calculate the hydrogen atoms. If the
    option is set to "reduce", the Reduce software will be run to do this task.
    :param conf_path: Path to the config file.
    :return: the command-line string.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('external', 'calc_hydrogens')


def get_motifs_app(conf_path):
    """
    Get the command-line string that will be run in order to calculate the motifs. If the option is
    set to "fanmod", the Fanmod software will be used to do this task.
    :param conf_path: Path to the config file.
    :return: the command-line string.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('external', 'calc_motifs')
