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
    config.add_section('fanmod')
    config.set('fanmod', '# = algorithm options', '')
    config.set('fanmod', '# number of samples used to determine approx. # of subgraphs', '[100000]')
    config.set('fanmod', 'num_samples_for_approx_num_subgraphs', '100000')
    config.set('fanmod', '# full enumeration? 1(yes)/0(no)', '[1]')
    config.set('fanmod', 'full_enumeration', '1')
    config.set('fanmod', '# = input file', '')
    config.set('fanmod', '# directed? 1(yes)/0(no)', '[1]')
    config.set('fanmod', 'directed', '0')
    config.set('fanmod', '# colored vertices? 1(yes)/0(no)', '[0]')
    config.set('fanmod', 'colored_vertices', '0')
    config.set('fanmod', '# colored edges? 1(yes)/0(no)', '[0]')
    config.set('fanmod', 'colored_edges', '0')
    config.set('fanmod', '# = random networks', '')
    config.set('fanmod', '# random type: 0(no regard)/1(global const)/2(local const)', '[2]')
    config.set('fanmod', 'random_type', '2')
    config.set('fanmod', '# regard vertex colors? 1(yes)/0(no)', '[0]')
    config.set('fanmod', 'regard_vertex_colors', '0')
    config.set('fanmod', '# regard edge colors? 1(yes)/0(no)', '[0]')
    config.set('fanmod', 'regard_edge_colors', '0')
    config.set('fanmod', '# reestimate subgraph number? 1(yes)/0(no)', '[0]')
    config.set('fanmod', 'reestimate_subgraph_num', '0')
    config.set('fanmod', '# number of random networks', '[1000]')
    config.set('fanmod', 'num_rand_networks', '1000')
    config.set('fanmod', '# number of exchanges per edge', '[3]')
    config.set('fanmod', 'num_exchanges_per_edge', '3')
    config.set('fanmod', '# number of exchange attempts per edge', '[3]')
    config.set('fanmod', 'num_exchange_attempts', '3')
    config.set('fanmod', '# = output file', '')
    config.set('fanmod', '# outfile format 1(ASCII - human readable)/0(CSV - for easy parsing)',
               '[1]')
    config.set('fanmod', 'outfile_format', '0')
    config.set('fanmod', '# create dumpfile? 1(yes)/0(no)', '[0]')
    config.set('fanmod', 'create_dumpfile', '1')
    config.add_section('external')
    config.set('external', 'calc_hydrogens', 'reduce')

    with open(out_path, 'w') as f:
        f.write('# These are the settings for the Protein-Protein Interaction Calculation Pipeline.'
                '\n# For further details on the values please refer to the documentation.\n# Lines '
                'starting with # are comments and will therefore be ignored.\n# {}\n\n'
                .format(datetime.datetime.now().strftime('%a %b %Y %H:%M:%S')))
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


def get_fanmod_num_samples(conf_path):
    """
    Get the number of samples that will be used to determine the approximate number of subgraphs.
    :param conf_path: Path to the config file.
    :return: the number of samples.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'num_samples_for_approx_num_subgraphs')


def get_fanmod_enum(conf_path):
    """
    Check whether Fanmod should run a full enumeration to calculate the subgraphs or sampling.
    :param conf_path: Path to the config file.
    :return: 1 if full enumeration, 0 if sampling.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'full_enumeration')


def get_fanmod_directed(conf_path):
    """
    Check whether the edges are treated as directed or undirected edges.
    :param conf_path: Path to the config file.
    :return: 1 if directed, 0 if undirected.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'directed')


def get_fanmod_colored_vertices(conf_path):
    """
    Check whether the vertices are treated as colored vertices or not.
    :param conf_path: Path to the config file.
    :return: 1 if colored, 0 if not colored.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'colored_vertices')


def get_fanmod_colored_edges(conf_path):
    """
    Check whether the edges are treated as colored edges or not.
    :param conf_path: Path to config file.
    :return: 1 if colored, 0 if not colored.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'colored_edges')


def get_fanmod_random_type(conf_path):
    """
    Check which type of random networks should be generated.
    :param conf_path: Path to the config file.
    :return: 0 if no regard, 1 if global constant, or 2 if local constant.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'random_type')


def get_fanmod_regard_vert_color(conf_path):
    """
    Check if the vertex colors should be considered during the generation of random networks.
    :param conf_path: Path to the config file.
    :return: 1 if vertex colors are considered, 0 if not.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'regard_vertex_colors')


def get_fanmod_regard_edg_color(conf_path):
    """
    Check if the edge colors should be considered during the generation of random networks.
    :param conf_path: Path to the config file.
    :return: 1 if edge colors are considered, 0 if not.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'regard_edge_colors')


def get_fanmod_reestimate(conf_path):
    """
    Check if the number of subgraphs should be re-estimated.
    :param conf_path: Path to the config file.
    :return: 1 if subgraphs should be re-estimated, 0 if not.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'reestimate_subgraph_num')


def get_fanmod_num_rand_networks(conf_path):
    """
    Get the number of random networks that will be generated.
    :param conf_path: Path to the config file.
    :return: the number of random networks.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'num_rand_networks')


def get_fanmod_num_exchanges_per_edge(conf_path):
    """
    Get the number of exchanges per edge for the random networks.
    :param conf_path: Path to the config file.
    :return: the number of exchanges per edge.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'num_exchanges_per_edge')


def get_fanmod_num_exchange_attempts(conf_path):
    """
    Get the number of exchange attempts per edge for the random networks.
    :param conf_path: Path to the config file.
    :return: the number of exchange attempts per edge.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'num_exchange_attempts')


def get_fanmod_outfile_format(conf_path):
    """
    Get the output file format. The options are ASCII format or CSV format.
    :param conf_path: Path to the config file.
    :return: 1 if ASCII output file format, 0 if CSV.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'outfile_format')


def get_fanmod_dumpfile(conf_path):
    """
    Check if a dumpfile should be generated. The dumpfile will contain a list of all subnetworks
    found in the input network including their adjacency matrix and the participating vertices.
    :param conf_path: Path to the config file.
    :return: 1 if dumpfile is created, 0 if not.
    """
    config = ConfigParser.SafeConfigParser()
    config.read(conf_path)
    return config.get('fanmod', 'create_dumpfile')
