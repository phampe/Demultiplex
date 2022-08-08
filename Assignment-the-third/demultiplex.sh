#!/bin/bash

#SBATCH --account=bgmp                 ### change this to your actual account for charging
#SBATCH --partition=bgmp               ### queue to submit to
#SBATCH --job-name=demultiplex         ### job name
#SBATCH --output=peter_output_%j.out   ### file in which to store job stdout. j holds the job id
#SBATCH --error=hostname_%j.err        ### file in which to store job stderr
#SBATCH --nodes=1                      ### number of nodes to use
#SBATCH --cpus-per-task=12             ### number of cores for each task



python demultiplexing.py -ic 35 -rc 0 -o "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" -t "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" -r "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" 
