#!/bin/bash

#SBATCH --job-name=hbb_scan
#SBATCH --output=res_%j.txt
#SBATCH --error=err_%j.txt

#SBATCH -N 1
#SBATCH --time=3:00:00
#SBATCH --ntasks=48
#SBATCH --exclusive

source setup.sh
python analyses/h_bb/h_bb.py --maxFiles 4 --nThreads 48
