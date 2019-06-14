import os
import re
import logging
import constants
import hpc_connection
from pssh.exceptions import SessionError


class RavenHPCProcess(object):

    def __init__(self, process_name, connection_cfg_dict=None):
        """
        :param process_name: 'raven' or 'ostrich'
        :param connection_cfg_dict: should contain at least two entries:
          - src_data_path: directory containing the input data
          - ssh_key_filename: filename of the private ssh key belonging to HPC user, e.g. crim01
        """
        self.logger = logging.getLogger(__name__)
        self.logger.debug("Creating connection")
        self.hpc_connection = hpc_connection.HPCConnection(connection_cfg_dict)
        self.process_name = process_name
        self.live_job_id = None
        self.template_path = constants.template_path
        self.shub_hostname = constants.shub_server
        self.last_progress = 0

    def check_connection(self):
        """
        Verifies that 1) hpc account is accessible and 2) the singularity registry is up
        """
        status, msg = self.hpc_connection.check_connection()
        if not status:
            return status, msg

        response = os.system("ping -c 1 -q " + self.shub_hostname + " > /dev/null 2>&1")
        if response == 0:
            return True, None
        else:
            return False, "shub server down"

    def submit(self, dataset, est_duration="01:00:00"):
        """
        Submits a job to slurm.
        1) copies the input data files
        2) creates the batch script needed by slurm
        3) queues the job
        :param dataset: name of the dataset used in the experiment, e.g. 'mohyse-salmon', 'hmets-salmon', etc.
        :param est_duration: estimated job duration, in format hh:mm:ss
        """
        self.logger.debug("Submitting a job (dataset {}, duration {}".format(dataset, est_duration))
        self.logger.debug("Copy data to hpc")
        self.live_job_id = 0
        try:
            self.hpc_connection.copy_data_to_remote(dataset)
            self.logger.debug("Copy batch script (exec {} selected)".format(self.process_name))
            remote_abs_script_fname = self.hpc_connection.copy_batchscript(self.process_name, est_duration, dataset,
                                                                           "batch_template.txt", self.shub_hostname)
            if self.process_name == "ostrich":
                # In addition, copy  raven script

                srcfilee = os.path.join(self.template_path, "Ost-RAVEN.sh")
                self.hpc_connection.copy_singlefile_to_remote(srcfilee, is_executable=True)
                self.hpc_connection.create_remote_subdir("model/output")
            self.logger.debug("Submit the job")
            jobid = self.hpc_connection.submit_job(remote_abs_script_fname)
            self.logger.debug("Job id = {}".format(jobid))
            self.live_job_id = jobid

        except Exception as e:
            if re.search("tar",repr(e)) is not None:
                raise Exception("Unable to submit job: bad input folder?")
            raise Exception("Unable to submit job: {}".format(e))

        self.last_progress = 0

    def retrieve(self, output_folder):
        """
        Copies the files generated during job execution to the output folder
        :param output_folder: folder where the files should be copied.
        :return:
        """
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        self.hpc_connection.copy_data_from_remote(self.live_job_id, output_folder)
    """
    def job_ended_normally(self):

        output_ok = True

        # Check for string such as DUE TO TIME LIMIT in slurm output file
        if self.hpc_connection.check_slurmoutput_for(r'DUE TO TIME LIMIT', self.live_job_id, ):
            output_ok = False
        print("todo: provide slurm output")
        return output_ok
    """
    def cancel(self):

        self.hpc_connection.cancel_job(self.live_job_id)

    def monitor(self):

        # job_status, progressfilecontent
        progressfile = None
        if self.process_name == 'raven':
            progressfile = 'out/Raven_progress.txt'
        if self.process_name == 'ostrich':
            progressfile = 'OstProgress0.txt'
        s = None
        try:

            s = self.hpc_connection.get_status(self.live_job_id)
            if s == "RUNNING":

                if progressfile is not None:
                    progressfile_content = self.hpc_connection.read_from_remote(progressfile)

                    for line in progressfile_content:

                        match_obj = re.search(r'progress\": (\d*)', line, re.M | re.I)
                        if match_obj:
                            progress = match_obj.group(1)
                            self.last_progress = progress

        except SessionError:
            self.logger.debug("Lost connection, reconnecting")
            self.hpc_connection.reconnect()

        return s, self.last_progress

    def cleanup(self):
        self.logger.debug("Cleaning up...")
        # remote
        self.hpc_connection.cleanup(self.live_job_id)

