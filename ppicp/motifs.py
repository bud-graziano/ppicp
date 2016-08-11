#!/usr/bin/env python

"""
ppicp.motifs
~~~~~~~~~~~~

Functions that handle the detection of motifs in networks.
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


def calculate_motifs(motif_size, in_file, out_file):
    """
    Detect motifs/subgraphs of a given motif size by invoking the ``fanmod_linux`` software. Only
    motifs of size 3 to 8 are supported.
    The standard parameters are taken from the config file.

    FANMOD command-line version

    input parameters (numbered in order):

    algorithm options
        1. subgraph (motif) size [default(as appears in original FANMOD GUI): 3]
        2. # of samples used to determine approx. # of subgraphs [100000]
        3. full enumeration? 1(yes)/0(no) [1]

    sampling probabilities will be input at end of command line, as their number may vary!

    input file
        4. infile name
        5. directed? 1(yes)/0(no) [1]
        6. colored vertices? 1(yes)/0(no) [0]
        7. colored edges? 1(yes)/0(no) [0]

    random networks
        8. random type: 0(no regard)/1(global const)/2(local const) [2]
        9. regard vertex colors? 1(yes)/0(no) [0]
        10. regard edge colors? 1(yes)/0(no) [0]
        11. re-estimate subgraph number? 1(yes)/0(no) [0]
        12. # of random networks [1000]
        13. # of exchanges per edge [3]
        14. # of exchange attempts per edge [3]

    output file
        15. outfile name
        16. outfile format 1(ASCII - human readable)/0(CSV - for easy parsing) [1]
        17. create dumpfile? 1(yes)/0(no) [0]

    sampling probabilities (if NOT full enumeration)
        18.->25 (depending on subgraph size) [.5]

    Each run also creates a log file in addition to the regular OUT file created by FANMOD.
    This log file (<outputfilename>.log) includes run times for each analyzed network (original and
    randomized networks) or run-preventing errors.
    A 'dump' file (<outputfilename>.dump) that contains a list of all subnetworks found in the input
    network as a list (including their adjacency matrix and the participating vertices) can also be
    created.

    :param motif_size: Size of the motif to detect.
    :param in_file: Path to the input file (must be correctly formatted).
    :param out_file: Path to the output files.
    :return: True if successful, False otherwise.
    """
    LOGGER.info('[MOTIF] Working on: %s and size %d}', in_file, motif_size)
    if not 2 < motif_size < 9:
        LOGGER.warning('Motif size not between 3-8.')
        return False

    config_file = initialize.CONF_FILE

    if platform.system() == 'Windows':
        raise NotImplementedError
    elif platform.system() == 'Linux':
        try:
            LOGGER.debug(subprocess.check_output([os.path.join(initialize.BIN_DIR, 'fanmod_linux'),
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
                                                  config.get_fanmod_num_exchanges_per_edge(
                                                      config_file),
                                                  config.get_fanmod_num_exchange_attempts(
                                                      config_file),
                                                  out_file,
                                                  config.get_fanmod_outfile_format(config_file),
                                                  config.get_fanmod_dumpfile(config_file)]))
            return True
        except (OSError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nMotif calculation failed.', err)
            return False
    else:
        LOGGER.critical('[ERROR] OS could not be determinded.')
        sys.exit(1)
