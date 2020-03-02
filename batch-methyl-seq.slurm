#!/bin/bash

# user email address
#SBATCH --mail-user=sxf279@case.edu

# mail is sent to you when the job starts and when it terminates or aborts
#SBATCH --mail-type=END,FAIL

# name of job
#SBATCH --job-name=methylseq

# standard output file
#SBATCH -o methylseq-pipe.%a.out
#SBATCH --time=100:30:00       # 20 minutes of wall time
#SBATCH --cpus-per-task=4      # 1 processor
#SBATCH --mem=50gb             # Assign 2gb memory; default in 1gb

# create array
#SBATCH --array=1-14

module load gcc/6.3.0
module load samtools/1.9
dirData="/scratch/users/sxf279/HIPCHBV/outs"
sampleID=$(find $dirData -name "*.bam")
sampleID=( $(echo $sampleID | \
		 tr ' ' '\n' | \
		 sort | \
		 uniq) )
bash script-methyl-seq.sh ${sampleID[$(($SLURM_ARRAY_TASK_ID - 1))]}
