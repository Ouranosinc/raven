from random import sample
from string import digits, ascii_uppercase, ascii_lowercase

from pssh.clients import ParallelSSHClient
from pssh.exceptions import AuthenticationException, UnknownHostException, ConnectionErrorException, SessionError, SFTPIOError

import pprint
import sys
import getopt
import time
import os
import shutil
import tempfile
import re
import json
import subprocess
#src_data_path = "./data/"
#do_cleanup = True
#template_dir = None
#dataset = None
#batch_cmd_template = "batch_template.txt"
#output_folder=""
from subprocess import run
from gevent import joinall

def rand_fname(length=8):
    chars = ascii_lowercase + ascii_uppercase + digits

    fname = 'tmp-' + ''.join(sample(chars, length))

    return fname

class HPCConnection(object):

    def __init__(self, external_init_dict=None):

        init_dict = {}
        if external_init_dict is not None:
            init_dict = external_init_dict

        self.hostname = init_dict.get("hostname", "cedar.computecanada.ca")
        self.user = init_dict.get("user", "lalondem")
        self.home_dir = init_dict.get("home_dir", os.path.join("/home",self.user))
        self.src_data_path = init_dict.get("src_data_path", "./data")
        self.template_path = init_dict.get("template_path", "./")
        self.batchscript_tmplt_filename = init_dict.get("batchscript_fname", "./batch_template.txt")

        self.client = ParallelSSHClient([self.hostname], user=self.user)
        self.remote_abs_working_folder = None
        self.remote_working_folder = None
        self.active_dataset_name = None
        self.live_job_id = None


    def check_connection(self):

        status = True
        msg = None

        try:
            self.client.run_command("ls")
        except (AuthenticationException, UnknownHostException, ConnectionErrorException) as e:
            status = False
            msg = str(e)

        return status, msg

    def copy_data_to_remote(self, dataset_, remote_temp_folder=None):
        """
        Copies data contained in a local directory over to a remote destination
        """

        remote_base_path = self.home_dir
        local_datapath = self.src_data_path
        if remote_temp_folder is None:
            remote_temp_folder = rand_fname()

        full_remote_path = os.path.join(remote_base_path, remote_temp_folder)
        remote_tar = os.path.join(full_remote_path, "data.tar")
        self.remote_abs_working_folder = full_remote_path
        self.remote_working_folder = remote_temp_folder
        self.active_dataset_name = dataset_
        self.client.run_command("mkdir "+full_remote_path)

    #    data_path_content = os.listdir(path=src_data_path)
    #    assert(len(data_path_content) == 1)
    #    df_basename = data_path_content[0]
        df_basename = dataset_

        print("system cmd: "+"tar cvf /tmp/"+remote_temp_folder+".tar -C "+ os.path.join(local_datapath,df_basename)+" .")
        os.system("tar cvf /tmp/"+remote_temp_folder+".tar -C "+ os.path.join(local_datapath,df_basename)+" .")
        try:
            g = self.client.scp_send("/tmp/"+remote_temp_folder+".tar", remote_tar)
            joinall(g, raise_error=True)
        except:
            raise Exception("scp_send failed")

        output = self.client.run_command("tar xvf "+remote_tar+" -C "+full_remote_path)
        errmsg = next(output[self.hostname]["stderr"], None)
        if errmsg is not None:
            print("Error: " + errmsg)
            raise Exception("Error untarring data file: "+errmsg)

        errmsg = next(output[self.hostname]["stdout"], None)
        if errmsg is not None:
            print("stdout: " + errmsg)

        os.remove("/tmp/"+remote_temp_folder+".tar")


# output files in base_dir/jobname/out
    def copy_data_from_remote(self, jobid, absolute_local_out_dir, cleanup_temp=True):

        absolute_tar_fname = os.path.join(self.remote_abs_working_folder, self.remote_working_folder+"_out.tar")
        absolute_output_data_path = os.path.join(self.remote_abs_working_folder, "out")
        stdout_file = os.path.join(self.home_dir,"slurm-"+jobid+".out")

        try:
            print("copy_data_from_remote: copying stdout: cp "+stdout_file+" "+ absolute_output_data_path)

            output = self.client.run_command("cp "+stdout_file+" "+ absolute_output_data_path)
            print("copy_data_from_remote: running tar")
            output = self.client.run_command("tar cf "+absolute_tar_fname+" -C "+absolute_output_data_path+" .")
            time.sleep(30)  # patch since run_command sems non-blocking
            output = self.client.run_command("du -sb " + absolute_tar_fname)
            line = ""
            for l in output[self.hostname].stdout:
                line += l
            print(line)
            tar_size = int(re.match("[0-9]*",line).group(0))


            print("copy_data_from_remote: {} bytes".format(tar_size))
            local_tar = "/tmp/"+self.remote_working_folder + "_out.tar"
           # g = self.client.scp_recv(absolute_tar_fname, local_tar)
            print("absolute tar fname={}".format(absolute_tar_fname))
            tries = 0
            while tries < 3:
                g = self.client.copy_remote_file(absolute_tar_fname, local_tar)  #scp_recv
                joinall(g, raise_error=True)
                print("end copy")
                output = subprocess.check_output("du -sb "+local_tar+"_"+self.hostname, shell=True)
                recv_tar_size = int(re.match("[0-9]*",output.decode('utf-8')).group(0))

                print("received: {} bytes".format(recv_tar_size))
                if recv_tar_size == tar_size:
                    break
                tries += 1

            if tries==3:
                raise Exception("Unable to copy tar file from remote end")

            if not os.path.exists(absolute_local_out_dir):
                #shutil.rmtree(absolute_local_out_dir)
                os.mkdir(absolute_local_out_dir)

            #os.mkdir(path.join(absolute_local_out_dir,jobname)
            print("tar xf "+local_tar+"_"+self.hostname+" -C "+absolute_local_out_dir)
            os.system("tar xf "+local_tar+"_"+self.hostname+" -C "+absolute_local_out_dir)
            if cleanup_temp:
                #os.remove(local_tar+"_" + self.hostname)
                pass

        except Exception as e:
            print("No output from job: {}".format(e))
        print("received")


    # executable_ is either raven or ostrich
    def copy_batchscript(self,executable_, datafile_basename, batch_tmplt_fname, shub_hostname):


        template_file = open(os.path.join(self.template_path, batch_tmplt_fname), "r")
        abs_remote_output_dir = os.path.join(self.remote_abs_working_folder, "out")
        tmplt = template_file.read()
        tmplt = tmplt.replace("ACCOUNT","def-fouchers")
        tmplt = tmplt.replace("DURATION","00:25:00")
        tmplt = tmplt.replace("TEMP_PATH", self.remote_abs_working_folder)
        tmplt = tmplt.replace("INPUT_PATH", self.remote_abs_working_folder)
        tmplt = tmplt.replace("OUTPUT_PATH", abs_remote_output_dir)
        tmplt = tmplt.replace("DATAFILE_BASENAME", datafile_basename)
        tmplt = tmplt.replace("SHUB_HOSTNAME", shub_hostname)
        tmplt = tmplt.replace("EXEC", executable_)

        #    subst_template_file, subst_fname = tempfile.mkstemp(suffix=".sh")
        subst_fname = self.remote_working_folder+".sh"
        file = open("/tmp/"+subst_fname, 'w')
        file.write(tmplt)
        file.close()

        self.client.run_command("mkdir " + abs_remote_output_dir)
        self.client.run_command("chmod 777 " + self.remote_abs_working_folder)
        self.client.run_command("chmod 777 " + abs_remote_output_dir)
        g = self.client.copy_file("/tmp/"+subst_fname, os.path.join(self.remote_abs_working_folder,subst_fname))
        joinall(g, raise_error=True)
        self.client.run_command("chmod ugo+x " + os.path.join(self.remote_abs_working_folder,subst_fname))
        os.remove("/tmp/"+subst_fname)

        if executable_ == "ostrich":
            # In addition, copy  raven script
            r = os.path.join(self.remote_abs_working_folder,"Ost-RAVEN.sh")
            srcfilee = os.path.join(self.template_path,"Ost-RAVEN.sh")
            g=self.client.copy_file(srcfilee, r)
            joinall(g, raise_error=True)
            self.client.run_command("chmod ugo+x " + r)
            self.client.run_command("mkdir " + self.remote_abs_working_folder + "/model/output")

            self.client.run_command("chmod 777 " + self.remote_abs_working_folder + "/model/output")


        return os.path.join(self.remote_abs_working_folder,subst_fname)

    def submit_job(self, script_fname):

        output = self.client.run_command("sbatch --parsable " + script_fname)
        errmsg = next(output[self.hostname]["stderr"], None)
        if errmsg is not None:
            raise Exception("Error: "+errmsg)

        self.live_job_id = next(output[self.hostname]["stdout"])

        return self.live_job_id

    def get_status(self, jobid, progressfile=None):

        cmd = "squeue -o '%T' -j {} --noheader".format(jobid)

        output = self.client.run_command(cmd)
        pp = pprint.PrettyPrinter(indent=4)
#        pp.pprint(output)
        status_output = []  # 1 line expected

        errmsg = next(output[self.hostname]["stderr"], None)
        if errmsg is not None:
            if errmsg.find("Invalid job id specified"):
                return "ERR", None
            print("stderr: "+errmsg)
            raise Exception("Error: "+errmsg)

        for line in output[self.hostname]["stdout"]:
            #print("stdout: " + line)
            status_output.append(line)

        if len(status_output) == 0:
            return "DONE", None
        elif len(status_output) >1:
            print("too many outputs:")
            pp.pprint(status_output)
            raise Exception("Error: too many outputs")
        s = status_output[0]

        progress = None
        if progressfile is not None:

            got_json = False
            try:
                    local_progress_file = os.path.join("/tmp", self.remote_working_folder + "_progress.json")
                    print("Trying to read progress file ->"+local_progress_file)
                    g = self.client.copy_remote_file(os.path.join(self.remote_abs_working_folder, progressfile), local_progress_file)  # scp_recv
                    joinall(g, raise_error=True)

                #load json

                    with open(local_progress_file+"_"+self.hostname) as f:
                        for line in f:
                            matchObj = re.search(r'progress\": (\d*)', line, re.M | re.I)
                            if matchObj:
                                progress = matchObj.group(1)
                        #progress_data = json.load(f)

            except SFTPIOError:
                print("SFTPIOError")
                pass
            except SessionError as e:
                    print(e)

        return s, progress

    def reconnect(self):

        self.client = ParallelSSHClient([self.hostname], user=self.user)

    def cleanup(self, jobid):

        output1 = self.client.run_command("rm -rf {}".format(self.remote_abs_working_folder))
        output2 = self.client.run_command("rm slurm-{}.out".format(jobid))

        try:
            print(next(output1[self.hostname_]["stdout"]))
            print(next(output2[self.hostname_]["stdout"]))
            print(next(output1[self.hostname_]["stderr"]))
            print(next(output2[self.hostname_]["stderr"]))
        except:
            pass


class RavenHPCProcess(object):

    def __init__(self, process_name, connection_cfg_dict=None):
        """

        :param process_name: 'raven' or 'ostrich'
        :param connection_cfg_dict:
        """

        self.hpc_connection = HPCConnection(connection_cfg_dict)
        self. process_name = process_name
        self.live_job_id = None
        self.shub_hostname = connection_cfg_dict.get("SHubHostname","132.217.141.54")

    def check_connection(self):

        status, msg = self.hpc_connection.check_connection()
        if status == False:
            return status, msg
        response = os.system("ping -c 1 " + self.shub_hostname)
        if response == 0:
            return True, None
        else:
            return False, "shub server down"



    def submit(self, dataset):

        self.hpc_connection.copy_data_to_remote(dataset)

        remote_abs_script_fname = self.hpc_connection.copy_batchscript(self.process_name, dataset, "batch_template.txt", self.shub_hostname)
        print("Running " + remote_abs_script_fname)
        jobid = self.hpc_connection.submit_job(remote_abs_script_fname)
        print("job id = " + jobid)
        self.live_job_id = jobid

    def retrieve(self, output_folder):

        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        self.hpc_connection.copy_data_from_remote(self.live_job_id, output_folder)

    def monitor(self):

        #job_status, progressfilecontent
        progressfile = None
        if self.process_name == 'raven':
            progressfile = 'out/RavenProgress.txt'
        if self.process_name == 'ostrich':
            progressfile = 'OstProgress0.txt'
        s = p = progress = None
        try:

            s, progress = self.hpc_connection.get_status(self.live_job_id, progressfile)
        except SessionError:
            print("reconnecting")
            self.hpc_connection.reconnect()

        return s, progress

    def cleanup(self):

        # remote
        #self.hpc_connection.cleanup(self.live_job_id)
        pass


def process_cmd(executable_, client, hostname, option, jobname="", jobid=""):

    if option == 'Retrieve':
#        output_folder_temp = path.join(output_folder, remote_temp_folder)
        if not os.path.exists(output_folder):
            #shutil.rmtree(output_folder_temp)
            os.mkdir(output_folder)

        copy_data_from_remote(client, base_dir=home_dir, jobname=jobname, jobid=jobid,
                              absolute_local_out_dir=output_folder)

    if option == 'Monitor':
        return get_status(client, hostname, jobid)

    if option == "Submit":
        remote_temp_folder = rand_fname()
        remote_path = "/home/" + user + "/"
        print("remote folder is " + remote_temp_folder)
        datafile_basename = copy_data_to_remote(client, src_data_path, dataset, remote_path, remote_temp_folder)

        script_fname = copy_batchscript(executable_, client, remote_path, remote_temp_folder, dataset, batch_cmd_template)
        print("Running " + script_fname)
        j = submit_job(client, script_fname)
        print("job id = " + j)
        return (j, remote_temp_folder)


def cleanup(client, hostname_, jobname, jobid):

    # remote
    print("cleanup:")
    output1 = client.run_command("rm -rf {}".format(path.join(home_dir,jobname)))
    output2 = client.run_command("rm slurm-{}.out".format(jobid))

    try:
        print(next(output1[hostname_]["stdout"]))
        print(next(output2[hostname_]["stdout"]))
        print(next(output1[hostname_]["stderr"]))
        print(next(output2[hostname_]["stderr"]))
    except:
        pass

    #locally
    os.remove(path.join("/tmp",jobname + "_out.tar_" + hostname_))


"""
Expected data structure (raven)
src_data_dir/datasetname/datasetname.rv?
"""






def newmainfct(argv):
    #        print("main.py command (Submit|Monitor|Retrieve)")

    try:
        opts, args = getopt.getopt(argv, "e:d:i:o:t:n")
    except getopt.GetoptError:

        print('main.py -e executable (raven|ostrich) -d workingdir -i dataset -o outputdir -t template_dir -n (no cleanup)')
        sys.exit(2)

    src_data_dir = "./"
    out_dir = "./"
    executable = None
    for opt, arg in opts:
        if opt == '-d':
            src_data_dir = arg
        elif opt == '-e':
            if arg == "raven":
                executable = "raven"
            elif arg == "ostrich":
                executable = "ostrich"
        elif opt == "-i":
            dataset = arg
        elif opt == "-o":
            out_dir = arg
        elif opt == "-t":
            template_dir = arg
        elif opt == "-n":
            do_cleanup = False

    if dataset is None or template_dir is None:
        print("Missing arg!")
        exit(2)
    if executable is None:
        print("Missing executable!")
        exit(2)

#    src_data_path = src_data_dir

    raven_process = RavenHPCProcess(executable,{"src_data_path":src_data_dir, "template_path":template_dir})

    status, msg = raven_process.check_connection()
    if status == True:
        print("Network connections are ok")
    else:
        print("Network error: "+msg)

#    jobinfo = process_cmd(executable, client,hostname,"Submit")
    raven_process.submit(dataset)

    job_finished = False
    while(not job_finished):

        time.sleep(60)
        try:

            out, p = raven_process.monitor()
            if out == "PENDING":
                print("Still pending")
            if out == "RUNNING":
                if p is not None:
                    print("Running ({}%)".format(p))
                else:
                    print("Running (?%)")
            if out == "DONE":
                job_finished = True
            if out == "ERR":
                print("error")
                job_finished = True
            if out is None:
                print("Temp error")

        except Exception as e:
            print("Exception @monitor")
            print(e)
            job_finished = True

    raven_process.retrieve(out_dir)
    raven_process.cleanup()

    print("done.")




def mainfct(argv):
    #        print("main.py command (Submit|Monitor|Retrieve)")

    try:
        opts, args = getopt.getopt(argv, "e:d:i:o:t:n")
    except getopt.GetoptError:

        print('main.py -e executable (raven|ostrich) -d workingdir -i dataset -o outputdir -t template_dir -n (no cleanup)')
        sys.exit(2)

    src_data_dir = "./"
    out_dir = "./"
    executable = None
    global dataset
    global template_dir
    global do_cleanup
    for opt, arg in opts:
        if opt == '-d':
            src_data_dir = arg
        elif opt == '-e':
            if arg == "raven":
                executable = "raven"
            elif arg == "ostrich":
                executable = "ostrich"
        elif opt == "-i":
            dataset = arg
        elif opt == "-o":
            out_dir = arg
        elif opt == "-t":
            template_dir = arg
        elif opt == "-n":
            do_cleanup = False

    if dataset is None or template_dir is None:
        print("Missing arg!")
        exit(2)
    if executable is None:
        print("Missing executable!")
        exit(2)

    global output_folder
    output_folder = out_dir
    global src_data_path
    src_data_path = src_data_dir


    hosts = [hostname]


    #global output_folder
    #output_folder = sys.argv[1]

    if True:

        jobinfo = process_cmd(executable, client,hostname,"Submit")

        job_finished = False
        while(not job_finished):

            time.sleep(60)
            try:
                out = process_cmd(executable, client, hostname, "Monitor", jobid=jobinfo[0] )
                if out == "PENDING":
                    print("Still pending")
                if out == "RUNNING":
                    print("Running")
                if out == "DONE":
                    job_finished = True
            except Exception as e:
                print("Exception @monitor")
                print(e)
                job_finished = True

        process_cmd(executable,client, hostname, "Retrieve", jobname=jobinfo[1], jobid=jobinfo[0])
        if do_cleanup:
            cleanup(client, hostname, jobname=jobinfo[1], jobid=jobinfo[0])


    print("done.")


# -d workingdir -i rootdatafile -o outdir
if __name__ == '__main__':


    newmainfct(sys.argv[1:])
