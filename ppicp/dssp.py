#!/usr/bin/env python

"""
ppicp.dssp
~~~~~~~~~~

Handles all the interaction with the DSSP server, i.e. calculating and downloading DSSP files.
"""

import json
import os
import Queue
import threading
import time
import sys

import requests

PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(PARENT_DIR)

from ppicp import initialize


DSSP_API_URL = r'http://www.cmbi.ru.nl/xssp/'

LOGGER = initialize.init_logger(__name__)


class DsspDownloader(threading.Thread):
    """
    Worker class that allows to calculate and download multiple DSSP files from the DSSP server in
    parallel.
    """
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.queue = queue

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            pdb_id, url, out_dir = self.queue.get()
            pdb_to_dssp(pdb_id, url, out_dir)
            self.queue.task_done()


def pdb_to_dssp(pdb_id, rest_url, out_dir):
    """
    Uploads a PDB file to the DSSP server, which calculates the corresponding DSSP file, and then
    downloads the generated DSSP file.

    :param out_dir: Where the DSSP files are saved.
    :param pdb_id: PDB file for which the DSSP file should be generated.
    :param rest_url: Url of the rest API of the DSSP server.
    :return: True if calculation and download was successful, False otherwise.
    """
    # Remember the PDB ID
    # tmp = pdb_file_path.split('/')
    # pdb_id = tmp[len(tmp) - 1].split('.')[0].lower()

    # Read the pdb file data into a variable
    try:
        files = {'file_': open(pdb_id, 'rb')}
    # Send a request to the server to create dssp data from the pdb file data.
    # If an error occurs, an exception is raised and the program exits. If the
    # request is successful, the id of the job running on the server is
    # returned.
        url_create = '{}api/create/pdb_file/dssp/'.format(rest_url)
        req = requests.post(url_create, files=files)
        req.raise_for_status()

        job_id = json.loads(req.text)['id']
        LOGGER.debug("Job %s submitted successfully. Id is: '%s'", pdb_id, job_id)

        # Loop until the job running on the server has finished, either successfully
        # or due to an error.
        ready = False
        while not ready:
            # Check the status of the running job. If an error occurs an exception
            # is raised and the program exits. If the request is successful, the
            # status is returned.
            url_status = '{}api/status/pdb_file/dssp/{}/'.format(rest_url, job_id)
            req = requests.get(url_status)
            req.raise_for_status()

            status = json.loads(req.text)['status']
            #print "Job status is: '{}'".format(status)

            # If the status equals SUCCESS, exit out of the loop by changing the
            # condition ready. This causes the code to drop into the `else` block
            # below.
            #
            # If the status equals either FAILURE or REVOKED, an exception is raised
            # containing the error message. The program exits.
            #
            # Otherwise, wait for five seconds and start at the beginning of the
            # loop again.
            if status == 'SUCCESS':
                LOGGER.info("Job %s status is: %s", pdb_id, status)
                ready = True
            elif status in ['FAILURE', 'REVOKED']:
                LOGGER.warning(json.loads(req.text)['message'])
                raise Exception(json.loads(req.text)['message'])
            else:
                time.sleep(5)
        else:
            # Requests the result of the job. If an error occurs an exception is
            # raised and the program exits. If the request is successful, the result
            # is returned.
            url_result = '{}api/result/pdb_file/dssp/{}/'.format(rest_url, job_id)
            req = requests.get(url_result)
            req.raise_for_status()
            result = json.loads(req.text)['result']

            # Save the result as a text file with file extension .dssp
            with open(out_dir + '/' + pdb_id[:4] + '.' + 'dssp', 'w') as res:
                res.write(result)
            # Return the result to the caller, which prints it to the screen.
            return ready
    except OSError as err:
        LOGGER.error('%s \nMissing PDB file for %s', err, pdb_id)


def dssp_download(pdb_ids, out_dir, num_threads):
    """
    Main function to be run when this file is executed.
    Takes the path to a directory with PDB files from user input and downloads the corresponding
    DSSP file.

    :param out_dir: Where the DSSP files are saved.
    :param pdb_ids: Path to a directory containing PDB files.
    :param num_threads: determines how many files should be downloaded in parallel.
    :return: True if finished.
    """
    # Threaded execution of dssp file downloading to gain some speedup because it's very I/O heavy
    dssp_queue = Queue.Queue()
    for pdb_id in pdb_ids:
        dssp_queue.put((pdb_id, DSSP_API_URL, out_dir))
    for unused in xrange(num_threads):
        dssp_worker = DsspDownloader(dssp_queue)
        dssp_worker.start()
    dssp_queue.join()
    return True
