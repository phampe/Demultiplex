#!/bin/bash

#SBATCH --account=bgmp                 ### change this to your actual account for charging
#SBATCH --partition=bgmp               ### queue to submit to
#SBATCH --job-name=Neverland           ### job name
#SBATCH --output=peter_output_%j.out   ### file in which to store job stdout. j holds the job id
#SBATCH --error=hostname_%j.err        ### file in which to store job stderr
#SBATCH --time=0-09:00:00              ### wall-clock time limit, in minutes
#SBATCH --nodes=1                      ### number of nodes to use
#SBATCH --cpus-per-task=8              ### number of cores for each task



python first_part1.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" -l 101 -o "read_1"
### read 1
#python first_part1.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" -l 8 -o "index_1"
### index 1
#python first_part1.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -l 8 -o "index_2"
### index 2
#python first_part1.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -l 101 -o "read_2"
### read 2
