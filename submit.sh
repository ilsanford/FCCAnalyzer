#!/bin/bash

#SBATCH --job-name=test_hbb
#SBATCH --output=res_%j.txt
#SBATCH --error=err_%j.txt

#SBATCH -N 1
#SBATCH --time=10:00
#SBATCH --ntasks=32

source setup.sh
python analyses/h_bb/h_bb.py --maxFiles 2 --nThreads 32
