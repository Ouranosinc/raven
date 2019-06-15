#!/bin/bash

# submit with:
#       sbatch submit-seq-to-graham.sh     

#SBATCH --account=def-btolson
#SBATCH --mem-per-cpu=1024M                        # memory; default unit is megabytes
# #SBATCH --mail-user=juliane.mai@uwaterloo.ca     # email address for notifications
# #SBATCH --mail-type=FAIL                         # email send only in case of failure
#SBATCH --time=0-00:15                             # time (DD-HH:MM);  budget = 200
# #SBATCH --time=0-02:30                           # time (DD-HH:MM);  budget = 2000
# #SBATCH --time=1-01:00                           # time (DD-HH:MM);  budget = 20000
#SBATCH --job-name=ostrich-LOWRL                   # name of job in queque
./Ostrich_v20171219_GCC_Linux.exe                  # job
