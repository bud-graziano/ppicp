#!/usr/bin/env python

"""
Handles all the interaction with the PDB servers, i.e. downloading PDB files.
"""

import gzip
import os
import Queue
import threading
import urllib


PDB_URL_FTP = r'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/divided/'
PDB_URL_FTP_ALT = r'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/divided/'
ALT_PDB_URL = r'http://www.rcsb.org/pdb/files/'


class PdbDownloader(threading.Thread):
    """
    Worker class that allows to download multiple PDB files from PDB's servers in parallel.
    """
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.queue = queue

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            pdb_id = self.queue.get()
            success = retrieve_pdb_bio_assembly_file(pdb_id, PDB_URL_FTP)
            if not success:     # retry to download the file
                print('Re-trying to download file {}.'.format(pdb_id))
                retrieve_pdb_bio_assembly_file(pdb_id, ALT_PDB_URL)
            self.queue.task_done()


def retrieve_pdb_bio_assembly_file(pdb_id, url):
    """
    Download the biological assembly (biological unit) PDB file from the PDB's ftp server.
    :param pdb_id: PDB ID for which the PDB file should be downloaded.
    :param url: from where the PDB files are downloaded.
    :return: True if download was successful, False otherwise.
    """
    pdb_code = pdb_id.lower()
    final_file = pdb_code + '.pdb'

    # Filename on the server
    archive = '{}.pdb1.gz'.format(pdb_code)

    # Check if the url is ftp or http since the url has to be encoded differently in each case
    if url.startswith('ftp'):
        url = r'{}{}/{}'.format(url, pdb_code[1:3], archive)
    else:
        url = r'{}{}'.format(url, archive)

    # Retrieve the file
    try:
        urllib.urlretrieve(url, archive)

        # Show status message
        print("Retrieving PDB file: {}".format(pdb_code))

        # Uncompress the archive and delete it when done
        zip_file = gzip.open(archive, 'rb')
        with open(final_file, 'wb') as out:
            out.writelines(zip_file)
        zip_file.close()
        os.remove(archive)
        return True
    except (IOError, urllib.ContentTooShortError) as err:
        print('{}\nUnable to download file {}.'.format(err, pdb_id))
        return False


def pdb_download(pdb_ids, num_threads):
    """
    Main function to be run when this file is executed. Downloads PDB files.
    Takes a list of PDB IDs from user input and downloads the corresponding PDB files.
    :param pdb_ids: List of PDB-IDs.
    :param num_threads: determines how many files should be downloaded in parallel.
    :return: True if finished.
    """
    # Store PDB IDs in a queue
    pdb_queue = Queue.Queue()
    for pdb_id in pdb_ids:
        pdb_queue.put(pdb_id)
    for unused in xrange(num_threads):  # ugly why to spawn n threads
        PdbDownloader(pdb_queue).start()
    pdb_queue.join()
    return True
