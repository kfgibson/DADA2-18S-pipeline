#!/bin/bash

#SBATCH --account=def-genome
#SBATCH --output=%x-%j.out           # Output log to <job_name>-<job_id>.out
#SBATCH --mail-user=<kfgibson@sfu.ca> # Your username
#SBATCH --mail-type=BEGIN,END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --job-name=DADA2
#SBATCH --time=168:00:00              # Time limit hrs:min:sec
#SBATCH --nodes=1                    # Run on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --mem=384000mb                # Total memory limit

min_len=$1
trunc_len_R1=$2
trunc_len_R2=$3
overlap=$4

cd $SLURM_SUBMIT_DIR
echo "Working directory = `pwd`"
echo "Hostname          = `hostname`"
echo "Start time        = `date`"
echo "Minimum read length: $min_len"
echo "R1 truncation length: $trunc_len_R1"
echo "R2 truncation length: $trunc_len_R2"
echo "Minimum overlap: $overlap"

#run stuff
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0

Rscript /home/kfgibson/scripts/DADA2_workflow_part2.R reads.fofn $min_len $trunc_len_R1 $trunc_len_R2 $overlap
Rscript /home/kfgibson/scripts/DADA2_workflow_part3.R

cut -f3  ASVs_species_pr2.tsv | paste ASVs_taxonomy_pr2.tsv - > ASVs_taxonomy_species_pr2.tsv
cut -f3  ASVs_species_silva.tsv | paste ASVs_taxonomy_silva.tsv - > ASVs_taxonomy_species_silva.tsv


RC=$?
echo "Job finished with exit code ${RC} at: `date`"
exit ${RC}
