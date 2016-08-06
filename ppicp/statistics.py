#!/usr/bin/env python

"""
A collection of functions that retrieve information from the files generated by the PPI calculations
which can be used for statistics reports.
"""

import re
import csv

def count_edges_in_ppi_aa_graph(path):
    """
    Takes the path to a file that stores edges for one protein complex (.fanmod file) and counts the
     edges.
    :param path: Path to file that contains the edges.
    :return: the number of edges in the PPI AA graph.
    """
    with open(path, 'r') as f:
        simple_file = f.readlines()
    return len(simple_file)


def count_vertices_in_ppi_aa_graph(path):
    """
    Takes the path to a file that stores edges for one protein complex (.fanmod file) and counts the
     vertices that are involved in edges.
    :param path: Path to file that contains the vertices.
    :return: the number of vertices in the PPI AA graph.
    """
    with open(path, 'r') as f:
        simple_file = f.read()

    vertices = re.split('\n| ', simple_file)
    vertices_id = []
    for v in vertices:
        if not v == '':
            vertices_id.append(int(v))
    try:
        return max(vertices_id) + 1  # +1 because .fanmod starts counting its ids at 0
    except ValueError as err:
        print "{}\nFile was most likely empty.".format(err)
        return None


def count_aas_in_aa_graph(aagraph_path):
    """
    Takes the path to a aagraph.gml file and counts how many vertices/amino acids are part of that
    protein in total.
    :param aagraph_path: Path to file that contains the vertices/amino acids.
    :return: the number of vertices/amino acids in that AA graph.
    """
    with open(aagraph_path, 'r') as f:
        aagraph = f.readlines()

    vertices_ids = []
    for v in aagraph:
        if v.startswith("    id"):  # indention is important -> gets all the vertex ids
            vertices_ids.append(int(v.split("id ")[1].strip()))
    return max(vertices_ids)


def count_all_contacts(stats_file):
    """
    Takes the path to a .stats file and count the different types of detected contacts on atom
    level. These contacts are not the same as edges in the PPI AA graph, since we are now counting
    contacts on an atom-atom level while the PPI AA graph shows contacts on a residue leve.
    :param stats_file: Path to the .stats file.
    :return: dictionary of all contact types and their abundances.
    """
    all_contact_types = {'BBHB': 0, 'BBBH': 0, 'IVDW': 0, 'ISS': 0, 'BCHB': 0, 'BCBH': 0, 'CBHB': 0,
                         'CBBH': 0, 'CCHB': 0, 'CCBH': 0, 'BB': 0, 'CB': 0, 'BC': 0, 'CC': 0,
                         'BL': 0, 'LB': 0, 'CL': 0, 'LC': 0, 'LL': 0, 'NHPI': 0, 'PINH': 0,
                         'CAHPI': 0, 'PICAH': 0, 'CNHPI': 0, 'PICNH': 0, 'SHPI': 0, 'PISH': 0,
                         'XOHPI': 0, 'PIXOH': 0, 'PROCDHPI': 0, 'PIPROCDH': 0, 'CCACOH': 0,
                         'CCOCAH': 0, 'BCACOH': 0, 'BCOCAH': 0}

    with open(stats_file, 'r') as f:
        text = f.readlines()

    for line in text:
        if ':' in line:     # to make sure this a valid line and no empty line at file end.
            c_type = line.split(": ")[0]
            c_num = int(line.split(": ")[1].rstrip('\n'))
            all_contact_types[c_type] += c_num

    return all_contact_types


def count_atom_num_contacts(csv_file):
    """
    Takes the path to a *atom_atom_contacts.csv file and counts how many contacts one single atom
    is part of.
    :param csv_file: Path to the *atom_atom_contacts.csv file.
    :return: a dictionary with all atoms and the amount of contacts they part of.
    """
    pdb_id = csv_file[len(csv_file)-27:len(csv_file)-23]
    all_atom_contacts = {}
    with open(csv_file, 'rb') as csvfile:
        text = csv.reader(csvfile, delimiter=';')

        for line in text:
            # Find chain_id and atom_id for residue A.
            chain_id = re.compile(r'\b({0}.?)\b'.format('Chain='),
                                  flags=re.IGNORECASE).search(line[1]).group().split(r'=')[1]
            atom_id = re.compile(r'\b({0}.*?\d)\b'.format('Atom'),
                                 flags=re.IGNORECASE).search(line[1]).group().split(r'#')[1]

            # An atom id equal to 0 usually denotes H atoms added by the Reduce software which we
            # want to exclude here atm.
            if not atom_id == '0':
                key = (pdb_id, chain_id, atom_id)
                if key not in all_atom_contacts:
                    all_atom_contacts[key] = 1
                else:
                    all_atom_contacts[key] += 1

            # Find chain_id and atom_id for residue B.
            chain_id = re.compile(r'\b({0}.?)\b'.format('Chain='), flags=re.IGNORECASE).search(
                line[2]).group().split(r'=')[1]
            atom_id = re.compile(r'\b({0}.*?\d)\b'.format('Atom'), flags=re.IGNORECASE).search(
                line[2]).group().split(r'#')[1]

            # An atom id equal to 0 usually denotes H atoms added by the Reduce software which we
            # want to exclude here atm.
            if not atom_id == '0':
                key = (pdb_id, chain_id, atom_id)
                if not all_atom_contacts.has_key(key):
                    all_atom_contacts[key] = 1
                else:
                    all_atom_contacts[key] += 1

    return all_atom_contacts


def count_pis(pdb_file):
    """
    Takes the path to a PDB file and counts how many aromatic pi-ring systems are present.
    :param pdb_file: Path to a PDB file.
    :return: the number of pi-ring systems.
    """
    aas = []

    # Go through all SEQRES entries and save all amino acid names.
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('SEQRES'):
                aas.extend(line[19:71].split(" "))

    trp = aas.count("TRP")
    tyr = aas.count("TYR")
    phe = aas.count("PHE")

    return 2 * trp + tyr + phe  # 2 times Tryptophane since it has a aromatic 5 and 6 ring
    # print('TRP: {} | TYR: {} | PHE: {} | TOTAL: {}'.format(trp, tyr, phe, trp + tyr + phe))
