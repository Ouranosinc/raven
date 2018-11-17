from random import sample
from string import digits, ascii_uppercase, ascii_lowercase
from os import path
#from fabric import Connection
from pssh.clients import ParallelSSHClient
import pprint
import sys
import time
import os
import shutil
import tempfile

hostname = "cedar.computecanada.ca"
user = "lalondem"
home_dir = path.join("/home",user)
src_data_path = "./data/"
batch_cmd_template = "batch_template.txt"
output_folder=""
from subprocess import run
from gevent import joinall

def rand_fname(length=8):
    chars = ascii_lowercase + ascii_uppercase + digits

    fname = 'tmp-' + ''.join(sample(chars, length))

    return fname

def copy_data_to_remote(c, datapath, rootdir, rootfname):

    c.run_command("mkdir "+rootdir+rootfname)
#    for root, dirs, files in os.walk(datapath):
#       for filename in files:
#            print("Transferring "+filename)
#            c.put(path.join(root,filename), remote=tempdest)  # in remote home directory for now
#       for dirname in dirs:
#            print("Transferring "+dirname)
#           c.put(path.join(root,dirname), remote=tempdest)  # in remote home directory for now

    data_path_content = os.listdir(path=src_data_path)
    assert(len(data_path_content) == 1)
    df_basename = data_path_content[0]

    os.system("tar cvf "+rootfname+".tar -C "+ path.join(datapath,df_basename)+" .")
    g=c.scp_send(rootfname+".tar", rootdir+rootfname+"/data.tar")
    joinall(g, raise_error=True)
    destpath = path.join(rootdir+rootfname,"data.tar")
    print(destpath)
    output = c.run_command("tar xvf "+destpath+" -C "+rootdir+rootfname)
    errmsg = next(output[hostname]["stderr"], None)
    if errmsg is not None:
        print("Error: " + errmsg)
    errmsg = next(output[hostname]["stdout"], None)
    if errmsg is not None:
        print("stdout: " + errmsg)

    print("tar xvf "+destpath+" -C "+rootdir+rootfname)
    os.remove(rootfname+".tar")
    #print(next(output[hostname]["stdout"]))
    return df_basename

# output files in base_dir/jobname/out
def copy_data_from_remote(c, base_dir, jobname, jobid, absolute_local_out_dir):

    absolute_tar_fname = path.join(base_dir, jobname, jobname+"_out.tar")
    absolute_output_data_path = path.join(base_dir, jobname,"out")
    stdout_file = path.join(home_dir,"slurm-"+jobid+".out")

    try:
        print("copy_data_from_remote: copying stdout: cp "+stdout_file+" "+ absolute_output_data_path)

        output = c.run_command("cp "+stdout_file+" "+ absolute_output_data_path)
        print("copy_data_from_remote: running tar")
        output = c.run_command("tar cf "+absolute_tar_fname+" -C "+absolute_output_data_path+" .")
        print("copy_data_from_remote: recv")
        #g = c.scp_recv(absolute_tar_fname, path.join(absolute_local_out_dir,jobname+"_out.tar"))
        g = c.copy_remote_file(absolute_tar_fname, jobname + "_out.tar")  #scp_recv
        joinall(g, raise_error=True)
        if os.path.exists(absolute_local_out_dir):
            shutil.rmtree(absolute_local_out_dir)
        os.mkdir(absolute_local_out_dir)
        print("tar xf /media/sf_develop/pavics-hydro/2cc/"+jobname+"_out.tar_"+hostname+" -C "+absolute_local_out_dir)
        os.system("tar xf /media/sf_develop/pavics-hydro/2cc/"+jobname+"_out.tar_"+hostname+" -C "+absolute_local_out_dir)

    except:
        print("No output from job")
    print("received")



def copy_batchscript(c, rootdir, rootfname, datafile_basename, batch_tmplt_fname):

    template_file = open(batch_tmplt_fname, "r")
    tmplt = template_file.read()
    tmplt = tmplt.replace("ACCOUNT","def-fouchers")
    tmplt = tmplt.replace("DURATION","00:20:00")
    tmplt = tmplt.replace("TEMP_PATH", rootdir+rootfname)
    tmplt = tmplt.replace("INPUT_PATH", rootdir + rootfname)
    tmplt = tmplt.replace("OUTPUT_PATH", rootdir + rootfname+"/out")
    tmplt = tmplt.replace("DATAFILE_BASENAME", datafile_basename)

    #    subst_template_file, subst_fname = tempfile.mkstemp(suffix=".sh")
    subst_fname = rootfname+".sh"
    file = open(subst_fname, 'w')
    file.write(tmplt)
    file.close()

    c.run_command("mkdir " + rootdir + rootfname + "/out")
    c.run_command("chmod 777 " + rootdir + rootfname)
    c.run_command("chmod 777 " + rootdir + rootfname + "/out")
    g = c.copy_file(subst_fname,path.join(rootdir+rootfname,subst_fname))
    joinall(g, raise_error=True)
    c.run_command("chmod ugo+x "+path.join(rootdir+rootfname,subst_fname))
    os.remove(subst_fname)
    return path.join(rootdir+rootfname,subst_fname)

def submit_job(client, script_fname):

    output = client.run_command("sbatch --parsable " + script_fname)
    errmsg = next(output[hostname]["stderr"], None)
    if errmsg is not None:
        raise Exception("Error: "+errmsg)

    jobid = next(output[hostname]["stdout"])

    return jobid

def get_status(c, hostname, jobid):

    cmd = "squeue -o '%T' -j {} --noheader".format(jobid)

    output = c.run_command(cmd)
    pp = pprint.PrettyPrinter(indent=4)
#    pp.pprint(output)
    status_output = []  # 1 line expected

    errmsg = next(output[hostname]["stderr"], None)
    if errmsg is not None:
        if errmsg.find("Invalid job id specified"):
            return "DONE"
        print("stderr: "+errmsg)
        raise Exception("Error: "+errmsg)

    for line in output[hostname]["stdout"]:
        #print("stdout: " + line)
        status_output.append(line)

    if len(status_output) == 0:
        return "DONE"
    elif len(status_output) >1:
        print("too many outputs:")
        pp.pprint(status_output)
        raise Exception("Error: too many outputs")

    return status_output[0]


def process_cmd(client, hostname, option, jobname="", jobid=""):
    if option == 'Retrieve':
        if os.path.exists(output_folder):
            shutil.rmtree(output_folder)
        os.mkdir(output_folder)

        copy_data_from_remote(client, base_dir=home_dir, jobname=jobname, jobid=jobid,
                              absolute_local_out_dir="/media/sf_develop/pavics-hydro/2cc/"+output_folder)

    if option == 'Monitor':
        return get_status(client, hostname, jobid)

    if option == "Submit":
        rootfname = rand_fname()
        rootdir = "/home/" + user + "/"
        print("remote folder is " + rootfname)
        datafile_basename = copy_data_to_remote(client, src_data_path, rootdir, rootfname)

        script_fname = copy_batchscript(client, rootdir, rootfname, datafile_basename, batch_cmd_template)
        print("Running " + script_fname)
        j = submit_job(client, script_fname)
        print("job id = " + j)
        return (j, rootfname)


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
    os.remove(jobname + "_out.tar_" + hostname_)



def mainfct():
    #        print("main.py command (Submit|Monitor|Retrieve)")

    hosts = [hostname]
    client = ParallelSSHClient(hosts, user=user)

    #copy_data_from_remote(client, "/home/" + user, jobname="tmp-tA4UxjKu", jobid="13436409",
    #                      absolute_local_out_dir="/media/sf_develop/pavics-hydro/2cc/toto/")
    #exit(1)

#    if len(sys.argv) > 1:
#       option = sys.argv[1]
#        if len(sys.argv==3):
#            remote_datadir = sys.argv[2]
#        else:
#            remote_datadir = None
#
#        process_cmd(client,hostname,option,jobname=remote_datadir)
#
#    else:
    global output_folder
    output_folder = sys.argv[1]

    if True:

        jobinfo = process_cmd(client,hostname,"Submit")

        job_finished = False
        while(not job_finished):

            time.sleep(60)
            try:
                out = process_cmd(client, hostname, "Monitor", jobid=jobinfo[0] )
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

        process_cmd(client, hostname, "Retrieve", jobname=jobinfo[1], jobid=jobinfo[0])
        cleanup(client, hostname, jobname=jobinfo[1], jobid=jobinfo[0])


    print("done.")


if __name__ == '__main__':
    mainfct()