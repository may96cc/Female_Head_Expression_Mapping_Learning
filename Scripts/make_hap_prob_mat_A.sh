#!/bin/bash

#SBATCH -A bac                 # Account name
#SBATCH -p bac,general                  # Partition name
#SBATCH --cpus-per-task=1              # Number of CPUs per task
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --mem=60G                       # Memory allocation
#SBATCH -t 12:00:00                     # Time limit
#SBATCH -J mk_hap_prb_mat_A              # Job name
#SBATCH -e mk_hap_prb_mat_A_%j.err       # Error log file
#SBATCH --mail-user=myork@missouri.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -o mk_hap_prb_mat_A.out          # Standard output file

# Activate conda environment
source /cluster/software/SPACK/SPACK_v0.20_dev_a2/spack/opt/spack/linux-almalinux8-x86_64/gcc-12.3.0/miniconda3-4.10.3-c6moxpqnii2vbelazwz5onnnnsh3cbzm/bin/activate /home/may96c/data/conda_envs/female_head_exp_project

# Run the R script for DSPR B panel
Rscript make_hap_prob_mat_A.R

