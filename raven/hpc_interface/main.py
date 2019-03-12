import sys
import getopt
import time
import raven_process
import logging
from select import select

"""
Expected data structure (raven)
src_data_dir/datasetname/datasetname.rv?
"""


def newmainfct(argv):

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG, filename="hpclog.txt", filemode='w')

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
    print("(Hit <Return> to kill job)")
    raven_proc = raven_process.RavenHPCProcess(executable, {"src_data_path": src_data_dir,
                                                            "template_path": template_dir})

    status, msg = raven_proc.check_connection()
    if status:
        print("Network connections are ok")
    else:
        print("Network error: "+msg)

    # jobinfo = process_cmd(executable, client,hostname,"Submit")
    print("Submitting job...")
    raven_proc.submit(dataset, "00:30:00")
    print(raven_proc.live_job_id)

    job_finished = False
    abnormal_ending = False
    while not job_finished:

        time.sleep(60)
        try:

            out, p = raven_proc.monitor()
            print(out)
            if out == "RUNNING":
                print("{}%".format(p))
            if out == "COMPLETED":
                job_finished = True
            if out == "TIMEOUT" or out == "CANCELLED":
                print("Uhoh: job "+out)
                abnormal_ending = True
                job_finished = True
            if out is None:
                print("Temp error")
            #dr, dw, de = select([sys.stdin], [], [], 0)
            if sys.stdin in select([sys.stdin], [], [], 0)[0]:
                raven_proc.cancel()

        except Exception as e:
            print("Exception @monitor")
            print(e)
            #job_finished = True

    # Check if job ended  normally
    if abnormal_ending:
        print("Job ended abnormally")
    else:
        raven_proc.retrieve(out_dir)

    raven_proc.cleanup()

    print("done.")


# -d workingdir -i rootdatafile -o outdir
if __name__ == '__main__':

    newmainfct(sys.argv[1:])
