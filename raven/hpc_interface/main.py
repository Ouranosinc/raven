import sys
import getopt
import time
import raven_process
import logging

"""
Expected data structure (raven)
src_data_dir/datasetname/datasetname.rv?
"""


def newmainfct(argv):
    # print("main.py command (Submit|Monitor|Retrieve)")

    logger = logging.getLogger()
    h = logging.FileHandler("hpclog.txt", mode='w')

    try:
        opts, args = getopt.getopt(argv, "e:d:i:o:t:n")
    except getopt.GetoptError:

        print('main.py -e executable (raven|ostrich) -d workingdir -i dataset -o outdir -t template_dir -n (nocleanup)')
        sys.exit(2)

    src_data_dir = "./"
    out_dir = "./"
    dataset = executable = template_dir = None
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

    # src_data_path = src_data_dir

    raven_proc = raven_process.RavenHPCProcess(executable, {"src_data_path": src_data_dir,
                                                            "template_path": template_dir})

    status, msg = raven_proc.check_connection()
    if status:
        print("Network connections are ok")
    else:
        print("Network error: "+msg)

    # jobinfo = process_cmd(executable, client,hostname,"Submit")
    raven_proc.submit(dataset)

    job_finished = False
    error_found = False
    while not job_finished:

        time.sleep(60)
        try:

            out, p = raven_proc.monitor()
            print(out)
            if out == "RUNNING":
                if p is not None:
                    print("Running ({}%)".format(p))
                else:
                    print("Running (?%)")
            if out == "COMPLETED":
                job_finished = True
            if out == "TIMEOUT" or out == "CANCELLED":
                print("error job cancelled/timeout")
                error_found = True
                job_finished = True
            if out is None:
                print("Temp error")

        except Exception as e:
            print("Exception @monitor")
            print(e)
            job_finished = True

    # Check if job ended  normally
    if error_found == False: #raven_proc.job_ended_normally():
        raven_proc.retrieve(out_dir)
    else:
        print("Job ended abnormally")

    raven_proc.cleanup()

    print("done.")


# -d workingdir -i rootdatafile -o outdir
if __name__ == '__main__':

    newmainfct(sys.argv[1:])
