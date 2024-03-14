#!/bin/bash

#SBATCH --job-name=test_hbb
#SBATCH --output=res_%j.txt
#SBATCH --error=err_%j.txt

#SBATCH -N 1
#SBATCH --time=12:00:00
#SBATCH --ntasks=64

source setup.sh
python analyses/h_bb/h_bb.py --maxFiles 16 --nThreads 64
