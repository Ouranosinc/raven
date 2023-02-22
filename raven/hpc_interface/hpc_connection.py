import logging
import os
import re
import subprocess
from random import sample
from string import ascii_lowercase, ascii_uppercase, digits

import constants
from gevent import joinall
from pssh.clients import ParallelSSHClient
from pssh.exceptions import (
    AuthenticationException,
    ConnectionErrorException,
    SCPError,
    UnknownHostException,
)


def rand_fname(length=8):
    chars = ascii_lowercase + ascii_uppercase + digits

    fname = "tmp-" + "".join(sample(chars, length))

    return fname


class HPCConnection:
    def __init__(self, external_init_dict=None):
        self.logger = logging.getLogger(constants.logging_name)
        init_dict = {}
        clsname = self.__class__.__name__
        if external_init_dict is not None:
            self.logger.debug(f"{clsname}: initializing from external dict")
            init_dict = external_init_dict
        else:
            self.logger.debug(f"{clsname}: initializing with default values")

        self.hostname = constants.hpc_hostname
        self.user = constants.user
        self.home_dir = os.path.join(constants.cc_working_dir, self.user)
        self.src_data_path = init_dict.get("src_data_path", "./data")
        self.template_path = constants.template_path
        self.logger.debug(
            f"Host being used is {self.hostname}, under username {self.user}"
        )
        self.keypath = init_dict.get("ssh_key_filename", constants.ssh_key_filename)
        self.client = ParallelSSHClient(
            [self.hostname], pkey=self.keypath, user=self.user, keepalive_seconds=300
        )
        self.remote_abs_working_folder = None
        self.remote_working_folder = None
        self.active_dataset_name = None
        self.live_job_id = None

    def check_connection(self):
        status = True
        msg = None
        self.logger.debug("Testing connection...")
        try:
            self.client.run_command("ls")
            self.logger.debug("... ok")
        except (
            AuthenticationException,
            UnknownHostException,
            ConnectionErrorException,
        ) as e:
            status = False
            msg = str(e)
            self.logger.debug(f"... failed ({msg})")

        return status, msg

    def copy_data_to_remote(self, dataset_, remote_temp_folder=None):
        """
        Copies data contained in a local directory over to a remote destination
        """

        self.logger.debug(
            "Copying data to remote location (from {} to {})".format(
                self.src_data_path, self.home_dir
            )
        )
        remote_base_path = self.home_dir
        local_datapath = self.src_data_path
        if remote_temp_folder is None:
            remote_temp_folder = rand_fname()

        full_remote_path = os.path.join(remote_base_path, remote_temp_folder)
        remote_tar = os.path.join(full_remote_path, "data.tar")
        self.remote_abs_working_folder = full_remote_path
        self.remote_working_folder = remote_temp_folder
        self.active_dataset_name = dataset_
        self.logger.debug(f"Creating remote folder {full_remote_path}")
        self.client.run_command("mkdir " + full_remote_path)

        #    data_path_content = os.listdir(path=src_data_path)
        #    assert(len(data_path_content) == 1)
        #    df_basename = data_path_content[0]
        df_basename = dataset_

        # self.logger.debug("system cmd: " + "tar cvf /tmp/" + remote_temp_folder + ".tar -C "
        #                   + os.path.join(local_datapath, df_basename) + " .")
        self.logger.debug(
            "system cmd: tar cvf /tmp/{}.tar -C {} .".format(
                remote_temp_folder, os.path.join(local_datapath, df_basename)
            )
        )

        os.system(
            "tar cf /tmp/"
            + remote_temp_folder
            + ".tar -C "
            + os.path.join(local_datapath, df_basename)
            + " ."
        )
        try:
            self.logger.debug("Copying data tar file")
            g = self.client.scp_send("/tmp/" + remote_temp_folder + ".tar", remote_tar)
            joinall(g, raise_error=True)
        except SCPError as e:
            self.logger.error(f"Copy failed (scp error {e})")
        except Exception as e:
            self.logger.error(f"Copy failed: {e}")
            raise Exception("scp_send failed")

        s = "tar xvf " + remote_tar + " -C " + full_remote_path
        self.logger.debug(f"Untarring remote data: {s}")

        output = self.client.run_command(s)
        self.client.join(output)

        errmsg = next(output[self.hostname]["stderr"], None)
        if errmsg is not None:
            self.logger.error("Error: " + errmsg)
            raise Exception("Error untarring data file: " + errmsg)

        errmsg = next(output[self.hostname]["stdout"], None)
        if errmsg is not None:
            self.logger.debug("stdout: " + errmsg)

        self.logger.debug(
            "Remove remote temp tar file " + "/tmp/" + remote_temp_folder + ".tar"
        )
        os.remove("/tmp/" + remote_temp_folder + ".tar")

    # output files in base_dir/jobname/out
    def copy_data_from_remote(self, jobid, absolute_local_out_dir, cleanup_temp=True):
        self.logger.debug("Copying data from remote")
        absolute_tar_fname = os.path.join(
            self.remote_abs_working_folder, self.remote_working_folder + "_out.tar"
        )

        absolute_output_data_path = os.path.join(self.remote_abs_working_folder, "out")
        stdout_file = os.path.join(self.home_dir, "slurm-" + jobid + ".out")
        self.logger.debug(f"  Remote data is located in {absolute_output_data_path}")
        self.logger.debug(f"  Slurm output file is {stdout_file}")

        try:
            self.logger.debug(f"  Copying slurm file to {absolute_output_data_path}")
            output = self.client.run_command(
                "cp " + stdout_file + " " + absolute_output_data_path
            )
            self.client.join(output)
            self.logger.debug(output)
            self.logger.debug("  Tarring remote folder")
            output = self.client.run_command(
                "tar cf "
                + absolute_tar_fname
                + " -C "
                + absolute_output_data_path
                + " ."
            )
            self.client.join(output)
            self.logger.debug(output)
            # time.sleep(30)  # patch since run_command sems non-blocking
            self.logger.debug("Picking up tar file size")
            output = self.client.run_command("du -sb " + absolute_tar_fname)
            self.client.join(output)
            self.logger.debug(output)
            line = ""
            for char in output[self.hostname].stdout:
                line += char
            # print(line)
            tar_size = int(re.match("[0-9]*", line).group(0))

            self.logger.info(f"{tar_size} bytes to copy from remote")
            local_tar = "/tmp/" + self.remote_working_folder + "_out.tar"
            # g = self.client.scp_recv(absolute_tar_fname, local_tar)
            self.logger.debug(f"Remote tar file is {absolute_tar_fname}")

            tries = 0
            while tries < 3:
                self.logger.debug("Copying tar file to /tmp")
                g = self.client.copy_remote_file(
                    absolute_tar_fname, local_tar
                )  # scp_recv
                joinall(g, raise_error=True)

                output = subprocess.check_output(
                    "du -sb " + local_tar + "_" + self.hostname, shell=True
                )
                recv_tar_size = int(re.match("[0-9]*", output.decode("utf-8")).group(0))

                self.logger.debug(f"Received: {recv_tar_size} bytes")
                if recv_tar_size == tar_size:
                    break
                tries += 1

            if tries == 3:
                raise Exception("Unable to copy tar file from remote end")

            if not os.path.exists(absolute_local_out_dir):
                # shutil.rmtree(absolute_local_out_dir)
                self.logger.debug(
                    "Local destination folder {} does not exist, creating".format(
                        absolute_local_out_dir
                    )
                )
                os.mkdir(absolute_local_out_dir)

            # os.mkdir(path.join(absolute_local_out_dir,jobname)
            self.logger.debug(f"Untarring received file to {absolute_local_out_dir}")
            os.system(
                "tar xf "
                + local_tar
                + "_"
                + self.hostname
                + " -C "
                + absolute_local_out_dir
            )
            if cleanup_temp:
                # print("todo: cleanup tmp file")
                os.remove(local_tar + "_" + self.hostname)

        except Exception as e:
            self.logger.error(f"Exception during file transfer from remote: {e}")

    def copy_singlefile_to_remote(
        self, local_filename, remote_path=".", is_executable=False
    ):
        r = os.path.join(
            self.remote_abs_working_folder,
            remote_path,
            os.path.basename(local_filename),
        )
        g = self.client.copy_file(local_filename, r)
        joinall(g, raise_error=True)
        if is_executable:
            self.client.run_command("chmod ugo+x " + r)

    def create_remote_subdir(self, remote_subdir):
        self.client.run_command(
            "mkdir -p " + os.path.join(self.remote_abs_working_folder, remote_subdir)
        )
        self.client.run_command(
            "chmod 777 " + os.path.join(self.remote_abs_working_folder, remote_subdir)
        )

    # executable_ is either raven or ostrich

    def copy_batchscript(
        self,
        executable_,
        guessed_duration,
        datafile_basename,
        batch_tmplt_fname,
        shub_hostname,
    ):
        template_file = open(os.path.join(self.template_path, batch_tmplt_fname))
        abs_remote_output_dir = os.path.join(self.remote_abs_working_folder, "out")
        tmplt = template_file.read()
        tmplt = tmplt.replace("ACCOUNT", constants.cc_account_info)
        tmplt = tmplt.replace("DURATION", guessed_duration)
        tmplt = tmplt.replace("TEMP_PATH", self.remote_abs_working_folder)
        tmplt = tmplt.replace("INPUT_PATH", self.remote_abs_working_folder)
        tmplt = tmplt.replace("OUTPUT_PATH", abs_remote_output_dir)
        tmplt = tmplt.replace("DATAFILE_BASENAME", datafile_basename)
        tmplt = tmplt.replace("SHUB_HOSTNAME", shub_hostname)
        tmplt = tmplt.replace("EXEC", executable_)

        # subst_template_file, subst_fname = tempfile.mkstemp(suffix=".sh")
        subst_fname = self.remote_working_folder + ".sh"
        file = open("/tmp/" + subst_fname, "w")
        file.write(tmplt)
        file.close()

        self.client.run_command("mkdir " + abs_remote_output_dir)
        self.client.run_command("chmod 777 " + self.remote_abs_working_folder)
        self.client.run_command("chmod 777 " + abs_remote_output_dir)
        g = self.client.copy_file(
            "/tmp/" + subst_fname,
            os.path.join(self.remote_abs_working_folder, subst_fname),
        )
        joinall(g, raise_error=True)
        self.client.run_command(
            "chmod ugo+x " + os.path.join(self.remote_abs_working_folder, subst_fname)
        )
        os.remove("/tmp/" + subst_fname)

        return os.path.join(self.remote_abs_working_folder, subst_fname)

    def submit_job(self, script_fname):
        self.logger.debug(f"Submitting job {script_fname}")
        # output = self.client.run_command("cd {}; ".format(self.home_dir) + constants.sbatch_cmd +
        #                                  " --parsable " + script_fname)
        output = self.client.run_command(
            "cd {}; {} --parsable {}".format(
                self.home_dir, constants.sbatch_cmd, script_fname
            )
        )
        self.client.join(output)
        errmsg = next(output[self.hostname]["stderr"], None)

        if errmsg is not None:
            for e in output[self.hostname]["stderr"]:
                errmsg += e + "\n"

            self.logger.error(f"  Error: {errmsg}")
            raise Exception("Error: " + errmsg)

        self.live_job_id = next(output[self.hostname]["stdout"])
        self.logger.debug(f"  Job id {self.live_job_id}")

        return self.live_job_id

    def read_from_remote(self, remote_filename):
        filecontent = []
        self.logger.debug("read_from_remote")
        retry = True
        # maybe remote file is being overwritten, try again if remote copy fails
        while True:
            try:
                local_filename = os.path.join(
                    "/tmp", self.remote_working_folder + "_progress.json"
                )
                g = self.client.copy_remote_file(
                    os.path.join(self.remote_abs_working_folder, remote_filename),
                    local_filename,
                )
                joinall(g, raise_error=True)
                suffixed_local_filename = local_filename + "_" + self.hostname
                self.logger.debug("  Opening copied file")
                with open(suffixed_local_filename) as f:
                    for line in f:
                        self.logger.debug(line)
                        filecontent.append(line)
                break
            #        except SFTPIOError:
            #            print("SFTPIOError")
            #            return False
            except Exception as e:
                if retry:
                    self.logger.debug(
                        f"exception {e}, retrying"
                    )  # pass # e.g. missing progress file as execution starts
                    retry = False
                else:
                    break

        self.logger.debug("End read_from_remote")
        return filecontent

    def get_status(self, jobid):
        """
        :param jobid:
        :return:
        """
        self.logger.debug("Inside get_status: executing sacct")
        cmd = constants.squeue_cmd + f" -j {jobid} -n -p -b"

        output = self.client.run_command(cmd)
        self.client.join(output)
        status_output = None  # 1 line expected

        errmsg = next(output[self.hostname]["stderr"], None)
        if errmsg is not None:
            for e in output[self.hostname]["stderr"]:
                errmsg += e + "\n"
            self.logger.debug(f"  stderr: {errmsg}")

            raise Exception("Error: " + errmsg)

        stdout_str = ""
        for line in output[self.hostname]["stdout"]:  # errmsg is None
            stdout_str += line + "\n"
            fields = line.split("|")
            if len(fields) >= 2:
                if fields[0] == jobid:
                    status_output = fields[1].split()[0]

        if status_output is None:
            raise Exception(f"Error parsing sacct output: {stdout_str}")

        if status_output not in [
            "COMPLETED",
            "PENDING",
            "RUNNING",
            "TIMEOUT",
            "CANCELLED",
        ]:
            raise Exception(f"Status error: state {status_output} unknown")

        return status_output

    def cancel_job(self, jobid):
        """
        :param jobid:
        :return:
        """
        cmd = constants.scancel_cmd + f" {jobid}"

        output = self.client.run_command(cmd)
        self.client.join(output)
        errmsg = next(output[self.hostname]["stderr"], None)
        if errmsg is not None:
            for e in output[self.hostname]["stderr"]:
                errmsg += e + "\n"
            self.logger.debug(f"  stderr: {errmsg}")

            raise Exception("Cancel error: " + errmsg)

        stdout_str = ""
        for line in output[self.hostname]["stdout"]:  # errmsg is None
            stdout_str += line + "\n"
        if len(stdout_str) > 0:
            raise Exception("Cancel error: " + stdout_str)

    def reconnect(self):
        self.client = ParallelSSHClient(
            [self.hostname], pkey=self.keypath, user=self.user, keepalive_seconds=300
        )

    """
    def check_slurmoutput_for(self, substr, jobid):

            slurmfname = "slurm-" + jobid + ".out"
            local_slurmfname = os.path.join("/tmp", slurmfname)
            stdout_file = os.path.join(self.home_dir, slurmfname)
            found = False
            try:
                g = self.client.copy_remote_file(stdout_file, local_slurmfname)
                joinall(g, raise_error=True)
                # scan file for substr
                with open(local_slurmfname + "_" + self.hostname) as f:
                    for line in f:
                        print("comparing {} with {}".format(substr,line))
                        match_obj = re.search(substr, line)
                        print(match_obj)
                        if match_obj:
                            found = True
                            print("found")

                os.remove(local_slurmfname + "_" + self.hostname)

            except Exception as e:
                print("Exception inside check_slurmoutput_for")
                print(e)
                pass

            return found
    """

    def cleanup(self, jobid):
        try:
            self.logger.debug("Deleting the remote folder")
            output1 = self.client.run_command(
                "rm -rf {}".format(
                    os.path.join(self.home_dir, self.remote_abs_working_folder)
                )
            )
            self.logger.debug("Deleting the slurm log file")
            logfilepath = os.path.join(self.home_dir, f"slurm-{jobid}.out")
            output2 = self.client.run_command(f"rm {logfilepath}")
            self.logger.debug("Deleting the local progress file")
            local_filename = os.path.join(
                "/tmp", self.remote_working_folder + "_progress.json"
            )
            suffixed_local_filename = local_filename + "_" + self.hostname
            os.remove(suffixed_local_filename)

            self.logger.debug(next(output1[self.hostname]["stdout"]))
            self.logger.debug(next(output2[self.hostname]["stdout"]))
            self.logger.debug(next(output1[self.hostname]["stderr"]))
            self.logger.debug(next(output2[self.hostname]["stderr"]))

        except Exception as e:
            self.logger.debug(f"Hmm file cleanup failed: {e}")
