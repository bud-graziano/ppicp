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
from ppicp import utilities

LOGGER = initialize.init_logger(__name__)
LOGGER_MOTIF = initialize.init_motif_logger(__name__)


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

    if utilities.is_empty_file(in_file):
        LOGGER.warning('%s is empty.', in_file)
        return False

    config_file = initialize.CONF_FILE

    if platform.system() == 'Windows':
        raise NotImplementedError
    elif platform.system() == 'Linux':
        try:
            command = [os.path.join(initialize.BIN_DIR, 'fanmod_linux'),
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
                       config.get_fanmod_dumpfile(config_file)]
            subp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stdout, stderr = subp.communicate()

            if stdout != '':
                LOGGER.debug(stdout)
            if stderr != '':
                LOGGER.debug(stderr)
                LOGGER_MOTIF.error('{%s|%d}\n %s', in_file, motif_size, stderr)

            retcode = subp.poll()
            if retcode:
                cmd = command
                raise subprocess.CalledProcessError(retcode, cmd, output=stdout)

            return True
        except (OSError, subprocess.CalledProcessError) as err:
            LOGGER.error('%s \nMotif calculation failed.', err)
            return False
    else:
        LOGGER.critical('[ERROR] OS could not be determinded.')
        sys.exit(1)


def combine_fanmod_files(in_dir, out_dir):
    """
    Combines several ``.fanmod`` files from a given input directory into one file.
    This means that all single graphs are combined into one network consisting of those single
    graphs.

    :param in_dir: Path to the input directory.
    :param out_dir: Path to the output directory.
    """
    # Get the path to all files ending with .fanmod in the specified directory.
    input_files = [os.path.join(in_dir, fanmod_file) for fanmod_file in os.listdir(in_dir) if
                   fanmod_file.endswith('.fanmod')]

    # Separate index to keep track of which ID in the new file belongs to which ID in what old file.
    separate_index = {}

    # The output that gets saved.
    out_combined = []

    # Will store the maximal ID that we need to add to a current ID to get a continuous list of IDs.
    max_id = 0

    for enum, fanmod_file in enumerate(input_files):
        # Get the ID of the input file
        tmp = os.path.split(fanmod_file)
        pdb_id = tmp[len(tmp) - 1].split('_')[0]

        # It is possible that there are empty .fanmod files, we don't want those.
        if not utilities.is_empty_file(fanmod_file):
            indices = []
            if enum > 0:
                # Read the input file.
                with open(fanmod_file, 'r') as inp:
                    lines = inp.read().splitlines()
                for line in lines:

                    # While going through all lines of the input file, store the IDs in a list.
                    # The file format is id1 id2 where id1 and id2 are vertices and two ids on the
                    # same line indicate an edge.
                    indices.append(int(line.split(" ")[0]))
                    indices.append(int(line.split(" ")[1]))

                    # Since we want to combine all input .fanmod files into one big file, we have to
                    # add the edges from every file to the big file. We also need a continuous ID
                    # range, i.e. we add the currently highest ID in the big file to the new IDs to
                    # be added to the big file, to ensure this behaviour.
                    edge_1 = (int(line.split(" ")[0])) + max_id
                    edge_2 = (int(line.split(" ")[1])) + max_id

                    # Save the final edge in a list that can later be saved to a file.
                    out_combined.append([edge_1, edge_2])

                    # We also keep track of which new ID corresponds to which old ID. Therefor we
                    # save the PDB ID and the old ID. This can be useful if we later want to convert
                    # the big .fanmod file back to the original files.
                    separate_index[edge_1] = [pdb_id, int(line.split(" ")[0])]
                    separate_index[edge_2] = [pdb_id, int(line.split(" ")[1])]

                # We stored all the IDs of the current .fanmod file in a list. Now take the maximal
                # ID, add 1 and add it to the previous ID. This way, we get the new max ID for the
                # next iteration.
                max_id += max(indices) + 1
            # If we are on the first iteration, we simply need to copy the input IDs to the output.
            # There is no need to add any value to the IDs. Still, we need to calculate the new
            # max ID for the next iteration.
            else:
                with open(fanmod_file, 'r') as inp:
                    lines = inp.read().splitlines()
                for line in lines:
                    # Same as above, this time we just don't add a max id yet.
                    indices.append(int(line.split(" ")[0]))
                    indices.append(int(line.split(" ")[1]))
                    edge_1 = (int(line.split(" ")[0]))
                    edge_2 = (int(line.split(" ")[1]))
                    out_combined.append([edge_1, edge_2])
                    separate_index[edge_1] = [pdb_id, int(line.split(" ")[0])]
                    separate_index[edge_2] = [pdb_id, int(line.split(" ")[1])]

                max_id = max(indices) + 1
        else:
            LOGGER.debug('Empty .fanmod file for PDB ID %s', pdb_id)

    # Save the new big .fanmod file.
    with open(os.path.join(out_dir, 'all_ppi_networks.comb_fanmod'), 'w') as output:
        for edge in out_combined:
            output.write(str(edge[0]) + ' ' + str(edge[1]) + '\n')

    # Save additional information regarding the vertices origin.
    with open(os.path.join(out_dir, 'suppl_all_ppi_networks.suppl_fanmod'), 'w') as suppl_output:
        for key, value in separate_index.iteritems():
            suppl_output.write(str(key) + ' ' + str(value) + '\n')
