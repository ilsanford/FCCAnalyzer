#!/bin/bash

#SBATCH --job-name=FCC-Hbb
#SBATCH --output=res_%j.txt
#SBATCH --error=err_%j.txt

#SBATCH -N 1
#SBATCH --time=12:00:00
#SBATCH --ntasks=64

source setup.sh
python analyses/h_bb/h_bb.py --nThreads 64
