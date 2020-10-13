#!/bin/bash

#SBATCH --account=def-genome
#SBATCH --output=%x-%j.out           # Output log to <job_name>-<job_id>.out
#SBATCH --mail-user=<kfgibson@sfu.ca> # Your username
#SBATCH --mail-type=BEGIN,END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --job-name=DADA2
#SBATCH --time=1:00:00              # Time limit hrs:min:sec
#SBATCH --nodes=1                    # Run on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --mem=48000mb                # Total memory limit

fprimer=$1
rprimer=$2
cutadapt_min=$3
cutadapt_max=$4

cd $SLURM_SUBMIT_DIR
echo "Working directory = `pwd`"
echo "Hostname          = `hostname`"
echo "Start time        = `date`"
echo "Fwd primer: $fprimer"
echo "Rev primer: $rprimer"
echo "Cutadapt min read length: $cutadapt_min"
echo "Cutadapt max read length: $cutadapt_max" 

#run stuff
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0
module load python/3.7.4

source ~/virtual_env/bin/activate
Rscript /home/kfgibson/scripts/DADA2_workflow_part1.R reads.fofn $fprimer $rprimer $cutadapt_min $cutadapt_max

RC=$?
echo "Job finished with exit code ${RC} at: `date`"
exit ${RC}
