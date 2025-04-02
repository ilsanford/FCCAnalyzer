#!/bin/bash
#SBATCH --job-name=h_ss_job      # job name
#SBATCH --output=h_ss_output.txt # output written to a text file
#SBATCH --error=h_ss_error.txt   # errors written to a text file
#SBATCH --time=05:00:00          # time limit for the job
#SBATCH -  n 1               # number of tasks
#SBATCH --cpus-per-task=4        # number of CPU cores per task
#SBATCH --mem=20G                # memory per node
#SBATCH --partition=submit       # partition name

# import the sourcing for using fccanalysis and location of the python folder
source /work/submit/jaeyserm/software/FCCAnalyses/setup.sh
export PYTHONPATH=$PYTHONPATH:/home/submit/isanford/FCCAnalyzer/

# running the script
fccanalysis run FCCAnalyzer/analyses/h_bb/h_ss.py

