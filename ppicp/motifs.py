#!/usr/bin/env python

"""
Functions that handle the detection of motifs in networks.
"""

import os
import platform
import subprocess
import sys

import config
import initialize


def calculate_motifs(motif_size, in_file, out_file):
    """
    Detect motifs/subgraphs of a given motif size by invoking the Fanmod software. Only motifs of
    size 3 to 8 are supported.
    The standard parameters are taken from the config file.

    ---------------------------------------------------------------------------
    FANMOD command-line version

    input parameters (numbered in order):
    algorithm options
        1- subgraph (motif) size [default(as appears in original FANMOD GUI): 3]
        2- # of samples used to determine approx. # of subgraphs [100000]
        3- full enumeration? 1(yes)/0(no) [1]
        sampling probabilities will be input at end of command line, as their number may vary!

    input file
        4- infile name
        5- directed? 1(yes)/0(no) [1]
        6- colored vertices? 1(yes)/0(no) [0]
        7- colored edges? 1(yes)/0(no) [0]

    random networks
        8- random type: 0(no regard)/1(global const)/2(local const) [2]
        9- regard vertex colors? 1(yes)/0(no) [0]
        10- regard edge colors? 1(yes)/0(no) [0]
        11- re-estimate subgraph number? 1(yes)/0(no) [0]
        12- # of random networks [1000]
        13- # of exchanges per edge [3]
        14- # of exchange attempts per edge [3]

    output file
        15- outfile name
        16- outfile format 1(ASCII - human readable)/0(CSV - for easy parsing) [1]
        17- create dumpfile? 1(yes)/0(no) [0]

    sampling probabilities (if NOT full enumeration)
        18->25 (depending on subgraph size) [.5]

    Each run also creates a log file in addition to the regular OUT file created by FANMOD.
    This log file (<outputfilename>.log) includes run times for each analyzed network (original and
    randomized networks) or run-preventing errors.
    A 'dump' file (<outputfilename>.dump) that contains a list of all subnetworks found in the input
    network as a list (including their adjacency matrix and the participating vertices) can also be
    created.
    ---------------------------------------------------------------------------
    :param motif_size: Size of the motif to detect.
    :param in_file: Path to the input file (must be correctly formatted).
    :param out_file: Path to the output files.
    :return: True if successful, False otherwise.
    """

    if not 2 < motif_size < 9:
        return False

    config_file = initialize.CONF_FILE

    if platform.system() == 'Windows':
        raise NotImplementedError
    elif platform.system() == 'Linux':
        try:
            print subprocess.check_output([os.path.join(initialize.BIN_DIR, 'fanmod_linux'),
                                           str(motif_size),
                                           config.get_fanmod_num_samples(config_file),
                                           config.get_fanmod_enum(config_file),
                                           in_file,
                                           config.get_fanmod_directed(config_file),
                                           config.get_fanmod_colored_vertices(config_file),
                                           config.get_fanmod_colored_edges(config_file),
                                           config.get_fanmod_random_type(config_file),
                                           config.get_fanmod_regard_vert_color(config_file),
                                           config.get_fanmod_regard_edg_color(config_file),
                                           config.get_fanmod_reestimate(config_file),
                                           config.get_fanmod_num_rand_networks(config_file),
                                           config.get_fanmod_num_exchanges_per_edge(config_file),
                                           config.get_fanmod_num_exchange_attempts(config_file),
                                           out_file,
                                           config.get_fanmod_outfile_format(config_file),
                                           config.get_fanmod_dumpfile(config_file)])
            return True
        except (OSError, subprocess.CalledProcessError) as err:
            print('{}\nMotif calculation failed.'.format(err))
            return False
    else:
        print('[ERROR] OS could not be determinded.')
        sys.exit(1)
