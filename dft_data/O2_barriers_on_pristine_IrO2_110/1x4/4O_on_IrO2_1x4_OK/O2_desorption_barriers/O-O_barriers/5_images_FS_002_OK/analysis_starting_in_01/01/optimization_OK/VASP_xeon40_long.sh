#!/bin/bash
#SBATCH -J job
#SBATCH -p xeon40	### xeon8; xeon16; xeon24; xeon40
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 50:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

srun /home/cat/zwang/bin/cattheory/vasp/5.4.4/vasp_std
