#!/usr/bin/env python

"""
ppicp.output_results
~~~~~~~~~~~~~~~~~~~~

Collection of functions to generate graphs/charts that visualize the results of the PPI
calculations.
"""

import os
import pygal


def bar_chart_types_of_contacts(all_contacts_dict, out_dir):
    """
    Generate and save a bar chart to give an overview of all detected contact types.

    :param all_contacts_dict: Dictionary containing all contact types.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = 'Types of Contacts'

    total_hb = all_contacts_dict['BBHB'] + all_contacts_dict['BBBH'] + all_contacts_dict['BCHB'] \
               + all_contacts_dict['BCBH'] + all_contacts_dict['CBHB'] + all_contacts_dict['CBBH'] \
               + all_contacts_dict['CCHB'] + all_contacts_dict['CCBH']

    total_lig = all_contacts_dict['BL'] + all_contacts_dict['LB'] + all_contacts_dict['CL'] \
                + all_contacts_dict['LC'] + all_contacts_dict['LL']

    total_pi = all_contacts_dict['NHPI'] + all_contacts_dict['PINH'] + all_contacts_dict['CAHPI'] \
               + all_contacts_dict['PICAH'] + all_contacts_dict['CNHPI'] \
               + all_contacts_dict['PICNH'] + all_contacts_dict['SHPI'] + all_contacts_dict['PISH']\
               + all_contacts_dict['XOHPI'] + all_contacts_dict['PIXOH'] \
               + all_contacts_dict['PROCDHPI'] + all_contacts_dict['PIPROCDH'] \
               + all_contacts_dict['CCACOH'] + all_contacts_dict['CCOCAH'] \
               + all_contacts_dict['BCACOH'] + all_contacts_dict['BCOCAH']

    total_vdw = all_contacts_dict['IVDW']

    total_disulphide = all_contacts_dict['ISS']

    bar_chart.add('Hydrogen Bonds', total_hb)
    bar_chart.add('Van der Waals Interactions', total_vdw)
    bar_chart.add(u'\u03C0-Effects', total_pi)
    bar_chart.add('Disulphide Bridges', total_disulphide)
    bar_chart.add('Ligands', total_lig)
    bar_chart.render_to_file(os.path.join(out_dir, 'all_contacts.svg'))


def bar_chart_hydrogen_verbose(all_contacts_dict, out_dir):
    """
    Generate and save a bar chart showing detailed information about detected hydrogen bonds.

    :param all_contacts_dict: Dictionary containing all contact types.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = 'Hydrogen Bonds (verbose)'

    bar_chart.add('BBHB', all_contacts_dict['BBHB'])
    bar_chart.add('BBBH', all_contacts_dict['BBBH'])
    bar_chart.add('BCHB', all_contacts_dict['BCHB'])
    bar_chart.add('BCBH', all_contacts_dict['BCBH'])
    bar_chart.add('CBHB', all_contacts_dict['CBHB'])
    bar_chart.add('CBBH', all_contacts_dict['CBBH'])
    bar_chart.add('CCHB', all_contacts_dict['CCHB'])
    bar_chart.add('CCBH', all_contacts_dict['CCBH'])

    bar_chart.render_to_file(os.path.join(out_dir, 'hb_verbose.svg'))


def bar_chart_hydrogen(all_contacts_dict, out_dir):
    """
    Generate and save a bar chart showing summarized information about detected hydrogen bonds.

    :param all_contacts_dict: Dictionary containing all contact types.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = 'Hydrogen Bonds'

    bb = all_contacts_dict['BBHB'] + all_contacts_dict['BBBH']
    cc = + all_contacts_dict['CCHB'] + all_contacts_dict['CCBH']
    bccb = all_contacts_dict['BCHB'] + all_contacts_dict['BCBH'] + all_contacts_dict['CBHB'] \
           + all_contacts_dict['CBBH']

    bar_chart.add('BB', bb)
    bar_chart.add('BC/CB', bccb)
    bar_chart.add('CC', cc)

    bar_chart.render_to_file(os.path.join(out_dir, 'hb.svg'))


def bar_chart_vdw(all_contacts_dict, out_dir):
    """
    Generate and save a bar chart showing information about detected van der Waals interactions.

    :param all_contacts_dict: Dictionary containing all contact types.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = 'Van der Waals Interactions'

    bar_chart.add('BB', all_contacts_dict['BB'])
    bar_chart.add('BC/CB', all_contacts_dict['BC'] + all_contacts_dict['CB'])
    bar_chart.add('CC', all_contacts_dict['CC'])

    bar_chart.render_to_file(os.path.join(out_dir, 'vdw.svg'))


def bar_chart_ligands(all_contacts_dir, out_dir):
    """
    Generate and save a bar chart showing information about detected ligand contacts.

    :param all_contacts_dir: Dictionary containing all contact types.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = 'Ligand Contacts'

    bar_chart.add('LB/BL', all_contacts_dir['LB'] + all_contacts_dir['BL'])
    bar_chart.add('LC/CL', all_contacts_dir['LC'] + all_contacts_dir['CL'])
    bar_chart.add('LL', all_contacts_dir['LL'])

    bar_chart.render_to_file(os.path.join(out_dir, 'ligand.svg'))


def bar_chart_pi_effects_verbose(all_contacts_dict, out_dir):
    """
    Generate and save a bar chart showing detailed information about detected pi-effects.

    :param all_contacts_dict: Dictionary containing all contact types.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = u'\u03C0-Effects (verbose)'

    bar_chart.add('NHPI', all_contacts_dict['NHPI'])
    bar_chart.add('PINH', all_contacts_dict['PINH'])
    bar_chart.add('CAHPI', all_contacts_dict['CAHPI'])
    bar_chart.add('PICAH', all_contacts_dict['PICAH'])
    bar_chart.add('CNHPI', all_contacts_dict['CNHPI'])
    bar_chart.add('PICNH', all_contacts_dict['PICNH'])
    bar_chart.add('SHPI', all_contacts_dict['SHPI'])
    bar_chart.add('PISH', all_contacts_dict['PISH'])
    bar_chart.add('XOHPI', all_contacts_dict['XOHPI'])
    bar_chart.add('PIXOH', all_contacts_dict['PIXOH'])
    bar_chart.add('PROCDHPI', all_contacts_dict['PROCDHPI'])
    bar_chart.add('PIPROCDH', all_contacts_dict['PIPROCDH'])
    bar_chart.add('CCACOH', all_contacts_dict['CCACOH'])
    bar_chart.add('CCOCAH', all_contacts_dict['CCOCAH'])
    bar_chart.add('BCACOH', all_contacts_dict['BCACOH'])
    bar_chart.add('BCOCAH', all_contacts_dict['BCOCAH'])

    bar_chart.render_to_file(os.path.join(out_dir, 'pi_verbose.svg'))


def bar_chart_pi_effects(all_contacts_dict, out_dir):
    """
    Generate and save a bar chart showing summarized information about detected pi-effects.

    :param all_contacts_dict: Dictionary containing all contact types.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = u'\u03C0-Effects'

    bar_chart.add('NHPI', all_contacts_dict['NHPI'] + all_contacts_dict['PINH'])
    bar_chart.add('CAHPI', all_contacts_dict['CAHPI'] + all_contacts_dict['PICAH'])
    bar_chart.add('ArgLysNHPI', all_contacts_dict['CNHPI'] + all_contacts_dict['PICNH'])
    bar_chart.add('CYSSHPI', all_contacts_dict['SHPI'] + all_contacts_dict['PISH'])
    bar_chart.add('XOHPI', all_contacts_dict['XOHPI'] + all_contacts_dict['PIXOH'])
    bar_chart.add('CAHOC', all_contacts_dict['CCACOH'] + all_contacts_dict['CCOCAH']
                  + all_contacts_dict['BCACOH'] + all_contacts_dict['BCOCAH'])
    bar_chart.add('PROCDHPI', all_contacts_dict['PROCDHPI'] + all_contacts_dict['PIPROCDH'])

    bar_chart.render_to_file(os.path.join(out_dir, 'pi.svg'))


def bar_chart_chem_props_3(chem_props_dict, out_dir, ppi=True):
    """
    Generate and save a bar chart showing information about the chemical properties of amino acids
    according to the three class system.

    :param chem_props_dict: Dictionary containing all chemical properties according to the three
     class system and how often they occur.
    :param ppi: Specify whether the input is from all amino acids or amino acids contributing to
     PPIs. That is a ``aagraph.gml`` or a ``aagraph_simple.id`` file, respectively.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = 'Amino Acid Chemical Properties (3 class)'

    bar_chart.add('Polar Uncharged', chem_props_dict['POLAR_UNCHARGED'])
    bar_chart.add('Polar Charged', chem_props_dict['POLAR_CHARGED'])
    bar_chart.add('Hydrophobic', chem_props_dict['HYDROPHOBIC'])
    bar_chart.add('Unknown', chem_props_dict['UNKNOWN'])

    if ppi:
        bar_chart.render_to_file(os.path.join(out_dir, 'chem_props_3_ppi.svg'))
    else:
        bar_chart.render_to_file(os.path.join(out_dir, 'chem_props_3_all.svg'))


def bar_chart_chem_props_5(chem_props_dict, out_dir, ppi=True):
    """
    Generate and save a bar chart showing information about the chemical properties of amino acids
    according to the five class system.

    :param chem_props_dict: Dictionary containing all chemical properties according to the five
     class system and how often they occur.
    :param ppi: Specify whether the input is from all amino acids or amino acids contributing to
     PPIs. That is a ``aagraph.gml`` or a ``aagraph_simple.id`` file, respectively.
    :param out_dir: Path to the output directory.
    """
    bar_chart = pygal.Bar()
    bar_chart.title = 'Amino Acid Chemical Properties (5 class)'

    bar_chart.add('Small Apolar', chem_props_dict['SMALL_APOLAR'])
    bar_chart.add('Hydrophobic', chem_props_dict['HYDROPHOBIC'])
    bar_chart.add('Polar', chem_props_dict['POLAR'])
    bar_chart.add('Negative Charge', chem_props_dict['NEGATIVE_CHARGE'])
    bar_chart.add('Positive Charge', chem_props_dict['POSITIVE_CHARGE'])
    bar_chart.add('Unknown', chem_props_dict['UNKNOWN'])

    if ppi:
        bar_chart.render_to_file(os.path.join(out_dir, 'chem_props_5_ppi.svg'))
    else:
        bar_chart.render_to_file(os.path.join(out_dir, 'chem_props_5_all.svg'))


def html_wrapper(out_dir,
                 date_time,
                 processed_files,
                 runtime, edges_ppi,
                 edges_atom,
                 aas_contributing,
                 aas_per_graph,
                 aas_contributing_per_graph,
                 atom_per_edge,
                 all_aas,
                 all_atom_atom_contacts,
                 all_resdiue_residue_contacts,
                 num_contacts_per_atom):
    """
    Generate and save a html file which shows all statistics after the PPI calculations are
    finished. This includes general information about the runtime, etc., as well as various graphs
    and other statistics.

    :param out_dir: Path to the output directory.
    :param date_time: The date and time the results are generated.
    :param processed_files: Number of processed files.
    :param runtime: Overall runtime.
    :param edges_ppi: Avg. number of edges in the PPI AA graph.
    :param edges_atom: Avg. number of edges on atom-atom level.
    :param aas_contributing: Number of amino acids contributing to PPIs.
    :param aas_per_graph: Avg. number of amino acids per AA graph.
    :param aas_contributing_per_graph: Avg. number of amino acids per PPI AA graph that contribute
     to PPIs.
    :param atom_per_edge: Avg. number of contacts on atom-atom level per edge in a PPI AA graph.
    :param all_aas: Number of all amino acids.
    :param all_atom_atom_contacts: Number of all contacts on atom level.
    :param all_resdiue_residue_contacts: Number of all contacts on residue level.
    :param num_contacts_per_atom: List of atoms and the number of contacts they are part of.
    """
    top_one = 'N.A.'
    top_two = 'N.A.'
    top_three = 'N.A.'
    top_four = 'N.A.'
    top_five = 'N.A.'
    try:
        top_one = num_contacts_per_atom[0]
        top_two = num_contacts_per_atom[1]
        top_three = num_contacts_per_atom[2]
        top_four = num_contacts_per_atom[3]
        top_five = num_contacts_per_atom[4]
    except IndexError:
        pass

    html_framework = r"""<!DOCTYPE html>
<html>
    <head>
        <meta charset="utf-8">
        <meta name="description" content="PPICP - Protein-Protein Interaction Calculation Pipeline">
        <meta name="author" content="The MolBI group">
        <link href='https://fonts.googleapis.com/css?family=Open+Sans:300italic,300,400italic,400,600italic,600,700italic,700,800italic,800' rel='stylesheet' type='text/css'>
        <title>PPICP - Protein-Protein Interaction Calculation Pipeline</title>

        <style type="text/css">
        h1 {{
            font-family: 'Open Sans', sans-serif;
            font-weight: 800;
        }}

        h2 {{
            font-family: 'Open Sans', sans-serif;
            font-weight: 700;
        }}

        h3 {{
            font-family: 'Open Sans', sans-serif;
            font-weight: 600;
        }}

        body {{
            font-family: 'Open Sans', sans-serif;
            font-weight: 400;
        }}

        div.container {{
            display:inline-block;
            width: 40vw;
            padding-left: 4%;
            padding-right: 1%;
          }}

        div.centering {{
            margin:0 auto;
            width: 50%;
        }}

        figure figcaption {{
            font-family: 'Open Sans', sans-serif;
            font-weight: 300;
            text-align: center;
        }}

        #inline {{
            display: inline;
            padding-right: 20px;
        }}

            #tiny {{
            font-weight: 300;
        }}
        </style>
    </head>

    <body>
        <h1 id="inline">Results</h1><span id="tiny">{}</span>

        <div>
            <h2>General Info</h2>
            {} processed PDB files<br>
            Runtime: {} min
        </div>

        <h2>Contact Statistics</h2>
        {} contacts on atom level.<br>
        {} contacts on residue level.
            <div class="centering">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/all_contacts.svg" width=100%/>
                    <figcaption>Fig 1. - Types of contacts.</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/hb_verbose.svg" />
                    <figcaption>Fig 2. - Hydrogen bonds (verbose).</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/hb.svg" />
                    <figcaption>Fig 3. - Hydrogen bonds.</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/vdw.svg" />
                    <figcaption>Fig 4. - Van der Waals interactions.</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/ligand.svg" />
                    <figcaption>Fig 5. - Ligand contacts.</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/pi_verbose.svg" />
                    <figcaption>Fig 6. - &pi;-effects (verbose).</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/pi.svg" />
                    <figcaption>Fig 7. - &pi;-effects.</figcaption>
                </figure>
            </div>


        <div>
        <h2>Motifs</h2>
        </div>


        <div>
        <h2>Additional Info</h2>
        <h3>Contact Statistics</h3>
            <p>&#216; number of contacts on residue level (of edges in PPI graph): {}<br>
            &#216; number of contacts on atom-atom level: {}<br>
            &#216; number of atom-atom contacts per edge in PPI graph: {}</p>
        <h3>Amino Acid Statistics</h3>
            <p>Number of all amino acids: {}<br>
            Number of amino acids contributing to PPIs: {}<br>
            &#216; number of amino acids per graph: {}<br>
            &#216; number of amino acids contributing to PPIs per graph: {}<br>
            <br>
            <h4>Chemical properties of all amino acids:</h4>
            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/chem_props_3_all.svg" />
                    <figcaption>Fig 8. - Amino acid chemical properties (3 class).</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/chem_props_5_all.svg" />
                    <figcaption>Fig 9. - Amino acid chemical properties (5 class).</figcaption>
                </figure>
            </div>
            <h4>Chemical properties of amino acids that contribute to PPIs:</h4>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/chem_props_3_ppi.svg" />
                    <figcaption>Fig 10. - Amino acid chemical properties (3 class).</figcaption>
                </figure>
            </div>

            <div class="container">
                <figure>
                    <embed type="image/svg+xml" src="./imgs/chem_props_5_ppi.svg" />
                    <figcaption>Fig 11. - Amino acid chemical properties (5 class).</figcaption>
                </figure>
            </div>



            <h4>Top 5 atoms with most contacts:</h4>
            Shown as: (PDB ID, Chain ID, Atom ID), &#35; contacts
            <ol>
                <li>{}</li>
                <li>{}</li>
                <li>{}</li>
                <li>{}</li>
                <li>{}</li>
            </ol>
        </div>
    </body>
</html>""".format(date_time, processed_files, runtime, all_atom_atom_contacts,
                  all_resdiue_residue_contacts, edges_ppi, edges_atom, atom_per_edge, all_aas,
                  aas_contributing, aas_per_graph, aas_contributing_per_graph, top_one, top_two,
                  top_three, top_four, top_five)

    with open(os.path.join(out_dir, 'results.html'), 'w') as f:
        f.write(html_framework)
