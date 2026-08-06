#!/usr/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem 1200G

source ~/.bashrc

source /opt/software/anaconda3/2023.03-py3.10/etc/profile.d/conda.sh
conda activate env_R4.4.0
Rscript "${1}"
