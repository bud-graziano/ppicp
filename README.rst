.. image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: http://protein-protein-interaction-calculation-pipeline.readthedocs.io/en/latest/
   :alt: PPICP Documentation

PPICP -- Protein-Protein Interaction Calculation Pipeline
=========================================================

Overview
--------
PPICP is a scientific software to calculate protein-protein interactions from RCSB PDB data. The pipeline includes the execution of the following steps:

1. Downloading of input data based on a list of PDB-IDs.

   a. Data from RCSB PDB.
   b. Data from DSSP.

2. Converting models in PDB files to separate chains.
3. Removing existing hydrogen atoms and re-adding them with the Reduce software.
4. Calculation of protein-protein interactions based on graph-theoretical methods.
5. Detection of motifs in the calculated PPI-graphs with the Fanmod software.
6. Compiling of statistics based on the PPI calculation output.

Note
~~~~
The project is still WIP.

Requirements
------------
- Python 2.7
- Works on Linux
- Works on Windows if only those steps are executed that do not depend on Reduce or Fanmod

Install
-------
First, download the source code to your machine either by

.. code-block:: bash

   $ git clone https://github.com/bud-graziano/ppicp.git

or you simply download the ZIP archive.

On your machine you can then run

.. code-block:: python

   python setup.py install

to install the pipeline.

The pipeline can then be executed with:

.. code-block:: bash

   $ ppicp --help


Documentation
-------------
The documentation is available online at `Read The Docs <http://protein-protein-interaction-calculation-pipeline.readthedocs.io/en/latest/>`_ and in the ``docs`` directory.

Related Projects
----------------
https://github.com/dfsp-spirit/vplg
