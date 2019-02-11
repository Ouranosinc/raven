import os
import logging
import constants
import hpc_connection
from pssh.exceptions import SessionError


class RavenHPCProcess(object):

    def __init__(self, process_name, connection_cfg_dict=None):
        """

        :param process_name: 'raven' or 'ostrich'
        :param connection_cfg_dict:
        """
        logger = logging.getLogger(__name__)


        logger.debug("Creating connection")
        self.hpc_connection = hpc_connection.HPCConnection(connection_cfg_dict)
        self.process_name = process_name
        self.live_job_id = None
        self.shub_hostname = connection_cfg_dict.get("SHubHostname", constants.shub_server)

    def check_connection(self):

        status, msg = self.hpc_connection.check_connection()
        if not status:
            return status, msg

        response = os.system("ping -c 1 " + self.shub_hostname)
        if response == 0:
            return True, None
        else:
            return False, "shub server down"

    def submit(self, dataset):

        self.hpc_connection.copy_data_to_remote(dataset)

        remote_abs_script_fname = self.hpc_connection.copy_batchscript(self.process_name, dataset,
                                                                       "batch_template.txt", self.shub_hostname)
        print("Running " + remote_abs_script_fname)
        jobid = self.hpc_connection.submit_job(remote_abs_script_fname)
        print("job id = " + jobid)
        self.live_job_id = jobid

    def retrieve(self, output_folder):

        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        self.hpc_connection.copy_data_from_remote(self.live_job_id, output_folder)

    def job_ended_normally(self):

        output_ok = True

        # Check for string such as DUE TO TIME LIMIT in slurm output file
        if self.hpc_connection.check_slurmoutput_for(r'DUE TO TIME LIMIT', self.live_job_id, ):
            output_ok = False
        print("todo: provide slurm output")
        return output_ok

    def monitor(self):

        # job_status, progressfilecontent
        progressfile = None
        if self.process_name == 'raven':
            progressfile = 'out/RavenProgress.txt'
        if self.process_name == 'ostrich':
            progressfile = 'OstProgress0.txt'
        s = progress = None
        try:

            s, progress = self.hpc_connection.get_status(self.live_job_id, progressfile)
        except SessionError:
            print("reconnecting")
            self.hpc_connection.reconnect()

        return s, progress

    def cleanup(self):

        # remote
        # self.hpc_connection.cleanup(self.live_job_id)
        pass
