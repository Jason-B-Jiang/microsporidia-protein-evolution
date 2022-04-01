#SBATCH --output=report
        #SBATCH --cpus-per-task=4
        #SBATCH --time=04:00:00
        #SBATCH --mem-per-cpu=2GB
        #SBATCH --mail-user=jasonbitan.jiang@mail.utoronto.ca
        #SBATCH --mail-type=BEGIN
        #SBATCH --mail-type=END
        #SBATCH --mail-type=FAIL
        
        module load StdEnv/2018.3
        module load nixpkgs/16.09
        module load gcc/7.3.0
        module load blast+/2.7.1
        module load parallel/20160722
        
        parallel < ../data/orthofinder_blast/blast_commands
        