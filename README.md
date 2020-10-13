This pipeline was desgined to run DADA2 on fastq files from 18S amplicon sequencing. 

TO RUN:

Step 1: If not already done, install cutadapt in a python virtual environment. After activating the virtual environment this can be done with the following command:

python3 -m pip install --user --upgrade cutadapt

Change the source ~/virtual_env/bin/activate in the 01-DADA2_job.sh script to the correct name of the virtual environment with cutadapt installed.


Step 2: Change the file path /home/kfgibson/scripts/DADA2_workflow_part1.R in the 01-DADA2_job.sh script to the correct file path on your machine. 


Step 3: Create the reads.fofn file. This should have the following format:

Sample_1	Sample_1_fwd_reads.fq	Sample_1_rev_reads.fq

Sample_2	Sample_2_fwd_reads.fq	Sample_2_rev_reads.fq


Step 4: Run 01-DADA2_job.sh:

sbatch ~/scripts/01-DADA2_job.sh <fwd primer> <rev primer> <min cutadapt read length> <max cutadapt read length>

For our data this can typically be run as: 

sbatch ~/scripts/01-DADA2_job.sh CCAGCASCCGCGGTAATWCC ACTTTCGTTCTTGATYRR 210 250


Step 5: Change the file paths /home/kfgibson/scripts/DADA2_workflow_part2.R and /home/kfgibson/scripts/DADA2_workflow_part3.R in the 02-DADA2_job.sh script to the correct file paths on your machine.


Step 6: Change the database locations in DADA2_workflow_part3.R to the correct location of these files on your machine.  


Step 7: Run 02-DADA2_job.sh:

sbatch ~/scripts/02-DADA2_job.sh <min DADA2 read length> <fed read truncation length> <rev read truncation length> <merge overlap length>

For our data this can typically be run as:   

sbatch ~/scripts/02-DADA2_job.sh 200 220 220 12 

However, you should check the read quality plots output from 01-DADA2_job.sh to determine the best parameters to run this with. 
