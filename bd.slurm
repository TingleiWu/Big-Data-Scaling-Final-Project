#!/bin/bash
#SBATCH --mail-user=jingting.xu@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --cpus-per-task=20

module load Anaconda3/2022.05
python /home/xuj31/bd.py

