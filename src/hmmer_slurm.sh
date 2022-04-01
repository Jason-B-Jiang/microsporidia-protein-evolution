#!/bin/bash

#SBATCH --output=hmmer_slurm_report
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=23
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-user=jasonbitan.jiang@mail.utoronto.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load StdEnv/2020
module load nixpkgs/16.09
module load gcc/7.3.0
module load hmmer/3.2.1
module load parallel/20160722

# Give permission for all hmmscan commands to be executed
chmod u+x ../data/hmmscan_cmds/*_S_cere

# Execute all commands in parallel
parallel ::: ../data/hmmscan_cmds/*_S_cere
