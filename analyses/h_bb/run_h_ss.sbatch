#!/bin/bash
#SBATCH -J h_ss_job                           # job name
#SBATCH -o h_ss_output.txt              # output written to a text file
#SBATCH -e h_ss_error.txt                # errors written to a text file
#SBATCH -t 05:00:00                           # time limit for the job
#SBATCH -n 1                                  # number of tasks
#SBATCH --cpus-per-task=8                     # number of CPU cores per task
#SBATCH --mem=20GB                            # memory per node
#SBATCH -p submit                             # partition name
#SBATCH --mail-type=BEGIN,END                 # send an email when starting and when done
#SBATCH --mail-user=isanford@terpmail.umd.edu # email address to send confirmation to

# import the sourcing for using fccanalysis and location of the python folder
source /work/submit/jaeyserm/software/FCCAnalyses/setup.sh
export PYTHONPATH=$PYTHONPATH:/home/submit/isanford/FCCAnalyzer/

# running the script
fccanalysis run FCCAnalyzer/analyses/h_bb/h_ss.py