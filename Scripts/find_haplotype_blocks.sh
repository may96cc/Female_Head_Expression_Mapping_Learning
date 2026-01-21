#!/bin/bash
#SBATCH -A bac                # Account name
#SBATCH -p general,bac                # Partition name
#SBATCH --cpus-per-task=8             # Number of CPUs (increase for faster Polars IO)
#SBATCH --nodes=1                      # Single node
#SBATCH --mem=60G                      # Memory allocation
#SBATCH -t 12:00:00                    # Time limit
#SBATCH -J find_haplotype_blocks             # Job name
#SBATCH -e find_haplotype_blocks_%j.err      # Error log file
#SBATCH -o find_haplotype_blocks_%j.out      # Standard output
#SBATCH --mail-user=myork@missouri.edu
#SBATCH --mail-type=FAIL,END

# ----------------------------
# Activate conda environment
# ----------------------------
source /cluster/software/SPACK/SPACK_v0.20_dev_a2/spack/opt/spack/linux-almalinux8-x86_64/gcc-12.3.0/miniconda3-4.10.3-c6moxpqnii2vbelazwz5onnnnsh3cbzm/bin/activate \
       /home/may96c/data/conda_envs/female_head_exp_project

# ----------------------------
# Optional: set Polars threads
# ----------------------------
export POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK

# ----------------------------
# Run Python script
# ----------------------------
/home/may96c/data/conda_envs/female_head_exp_project/bin/python find_haplotype_blocks.py


