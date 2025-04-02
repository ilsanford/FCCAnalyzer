#!/bin/bash
#SBATCH -J h_ss_job                           # job name
#SBATCH -o h_ss_output_%A_%a.txt              # output written to a text file with job number
#SBATCH -e h_ss_error_%A_%a.txt               # errors written to a text file with job number
#SBATCH -t 05:00:00                           # time limit for the job (HH:MM:SS)
#SBATCH -n 1                                  # number of tasks
#SBATCH --cpus-per-task=16                    # number of CPU cores per task
#SBATCH --mem=50GB                            # memory allocated
#SBATCH -p submit-gpu                         # partition name (pick between submit and submit-gpu)
#SBATCH --mail-type=BEGIN,END                 # send an email when done
#SBATCH --mail-user=isabellalynn622@gmail.com # email address to send confirmation to

# import the sourcing for using fccanalysis and location of the python folder
source /work/submit/jaeyserm/software/FCCAnalyses/setup.sh
export PYTHONPATH=$PYTHONPATH:/home/submit/isanford/FCCAnalyzer/

# running the script
fccanalysis run FCCAnalyzer/analyses/h_bb/h_ss.py --threads=16
