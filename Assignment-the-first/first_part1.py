#!/usr/bin/env python

#SBATCH --account=bgmp                 ### change this to your actual account for charging
#SBATCH --partition=bgmp               ### queue to submit to
#SBATCH --job-name=Neverland           ### job name
#SBATCH --output=peter_output_%j.out   ### file in which to store job stdout. j holds the job id
#SBATCH --error=hostname_%j.err        ### file in which to store job stderr
#SBATCH --time=0-03:00:00              ### wall-clock time limit, in minutes
#SBATCH --nodes=1                      ### number of nodes to use
#SBATCH --cpus-per-task=8              ### number of cores for each task

### This is to make a histogram for our 4 files


import argparse
import gzip

def get_args():
    parser=argparse.ArgumentParser(description="Getting files and length")
    parser.add_argument("-l", help= "length of the read" )
    parser.add_argument("-f", help= "filename" )
    parser.add_argument("-o", help= "filename" )
    return parser.parse_args()

args=get_args() ### gets you your function 

f= args.f
l= int(args.l)
o = args.o

dir = "/projects/bgmp/shared/2017_sequencing/"
read1 = "$dir/1294_S1_L008_R1_001.fastq.gz"
read2 = "$dir/1294_S1_L008_R4_001.fastq.gz"
index1 = "$dir/1294_S1_L008_R2_001.fastq.gz"
index2 = "$dir/1294_S1_L008_R3_001.fastq.gz"

import bioinfo

def init_list(lst:list,value:float=0.0)->list:
    for i in range(l):
        lst.append(value)
    return lst

my_list: list=[]
my_list = init_list(my_list)
### sets up my list

#print(my_list)

with gzip.open(f,"rt") as fq: #### allows you to read the zip file
    line_count=0

    for file_q in fq:
        file_line=file_q.strip("\n")
        line_count+=1
        
        if line_count%4==0:
            #print(file_line)

            phred_score=file_line

            for count, pc in enumerate(phred_score):
                #print(f'debug d6 - {count=} {len(my_list)=} {my_list=}')
                my_list[count]+=bioinfo.convert_phred(pc)
                
                    #print(my_list)
#print(line_count/4)
#print(my_list)
### At this point we have the summation of each base read 

for position, count in enumerate(my_list):
    mean=count/(line_count/4)
    my_list[position]=mean

# print(my_list)
### we now have our means at the given positions

import matplotlib.pyplot as plt
from matplotlib import pyplot as plt

x= range(len(my_list))
y=(my_list)

plt.bar(x,y)
plt.xlabel("Position")
plt.ylabel("Mean Quality Score")
plt.title("Mean distribution")
plt.savefig(o)

