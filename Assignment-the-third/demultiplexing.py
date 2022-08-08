#!/usr/bin/env python

### Demultiplexing 8/3/2022

import bioinfo
import argparse
import itertools
import gzip
import numpy
import matplotlib

def get_args():
    parser=argparse.ArgumentParser(description="Getting files and length")
    # parser.add_argument("-l", help= "length of the read" )
    parser.add_argument("-ic", help= "quality score cutoff index" )
    parser.add_argument("-rc", help= "quality score cutoff reads" )
    parser.add_argument("-o", help= "read one" )
    parser.add_argument("-t", help= "read two" )
    parser.add_argument("-r", help= "read three" )
    parser.add_argument("-f", help= "read four" )
    return parser.parse_args()

args=get_args() ### gets you your function 

# l = int(args.l)
ic = int(args.ic)
rc = int(args.rc)

file1 = args.o
file2 = args.t
file3 = args.r
file4 = args.f

#print(l,ic,rc)

### here we are creating a set for our indexes for use

dir="/projects/bgmp/shared/2017_sequencing/"

index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

index_set = set()

with open (index_file,"r") as i:
    for lines in i:
        if lines[0] == "s":
            continue
        index_set.add(lines.strip("\n").split("\t")[4])
        
# print(index_set)

################################################################################################################################################
### creating our function to write out the reverse complementary

def reverse_complementary(dna):
    reverse = ""
    comp_bases_dict = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    #print(comp_bases_dict)

    for base in dna:
        reverse += (comp_bases_dict[base])
        reverse_comp = reverse[::-1]
    #print(reverse_comp)
    return(reverse_comp)

###############################################################################################################################################
### creating count dictionary for indexes that do work and unknown

match_counts = {}
for index in index_set:
    match_counts[index]=0

unknown_count_dict= {"unknown" : 0}
# print(unknown_count_dict)
# print(match_counts)

###############################################################################################################################################
###creating our other dictionary for our complementary and our set for our complementary

reverse_set=set()
reverse_comp_dict = {}
# print(index_set)
for index in index_set:
    # print(index)
    reverse=reverse_complementary(index)
    reverse_set.add(reverse)
    # print(reverse,index)
    reverse_comp_dict[index] = reverse

# print(reverse_comp_dict)
# print(reverse_set)

###############################################################################################################################################
### create possible permutation of the possible hopping indexes
hopping_count = {}
permutation_set = set()
a=list(itertools.permutations(index_set,2))
for pair in a:
    a=pair[0]+"+"+pair[1]
    permutation_set.add(a)
    hopping_count[a] = 0

#print(hopping_count)

### the format for the dictionary is index1+index2
### this is very important because this is what we are going to need to figure out 
###############################################################################################################################################
### creating function to open up all 52 files
### these are the first four files that we would be making for our unknown and index hopped list

unknown_r1 = open("unknown_r1.fq","w")
# unknown_r1 = gzip.open("unknown_r1.fq","wt")
unknown_r2 = open("unknown_r2.fq","w")
# unknown_r2 = gzip.open("unknown_r2.fq","wt")

hopped_r1 = open("hopped_r1.fq","w")
# hopped_r1 = gzip.open("hopped_r1.fq","wt")
hopped_r2 = open("hopped_r2.fq","w")
# hopped_r2 = gzip.open("hopped_r2.fq","wt")

test_set = ("one","two")
open_file = {}
for i in index_set:
# for i in test_set:
    open_file[i] = open(i + "_R1.fq","wt"), open(i + "_R2.fq","wt")

#print(open_file)

### make sure you aren't an idiot and the file exists
###############################################################################################################################################
### opening up our four files
### we are using a while loop that goes through each of the file
### we are setting each of the four lines in the file to their own tuple which we will be able to access

### don't forget to change this before you work on the gzip file
# with open (file1,"rt") as r1, open (file2,"rt") as i1, open (file3,"rt") as i2, open(file4,"rt") as r2:
with gzip.open(file1,"rt") as r1, gzip.open(file2,"rt") as i1, gzip.open(file3,"rt") as i2, gzip.open(file4,"rt") as r2:
    while True:
        header = r1.readline().strip("\n")
        if header == "":
            ### end of file
            break
        seq = r1.readline().strip("\n")
        plus = r1.readline().strip("\n")
        score = r1.readline().strip("\n")
        r1_record = (header,seq,plus,score)
        #print(r1_record)
        header = r2.readline().strip("\n")
        seq = r2.readline().strip("\n")
        plus = r2.readline().strip("\n")
        score = r2.readline().strip("\n")
        r2_record = (header,seq,plus,score)
        # print(r2_record)
        header = i1.readline().strip("\n")
        seq = i1.readline().strip("\n")
        plus = i1.readline().strip("\n")
        score = i1.readline().strip("\n")
        i1_record = (header,seq,plus,score)
        # print(i1_record)
        header = i2.readline().strip("\n")
        seq = i2.readline().strip("\n")
        plus = i2.readline().strip("\n")
        score = i2.readline().strip("\n")
        i2_record = (header,seq,plus,score)
        # print(i2_record)

        if "N" in i1_record[1] or "N" in i2_record[1]:
            # print(i1_record[1])
            # print(i2_record[1])
            unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
            unknown_r2.write(r2_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r2_record[1]+ "\n"  + r2_record[2]+ "\n"  +  r2_record[3]+"\n")
            unknown_count_dict["unknown"]+=1
            ### Checks to see if there is a N in the line, if there is then we put into unknown
            
        # elif "N" in i2_record[1]:
        #     print(i1_record[1])
        #     print(i2_record[1]) 
        #     unknown_count_dict["unknown"]+=1
        #     ### Checks to see if there is a N in the line, if there is then we put into unknown
            
        elif i1_record[1] not in index_set:
            # print(i1_record[1])
            # print(i2_record[1])
            unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
            unknown_r2.write(r2_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r2_record[1]+ "\n"  + r2_record[2]+ "\n"  +  r2_record[3]+"\n")
            unknown_count_dict["unknown"]+=1  
            ### If the index isn't one of our listen indexes then we put the who thing into unknown 

        elif bioinfo.qual_score(i1_record[3]) < ic:
            # print(i1_record[1])
            # print(i2_record[1])
            unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
            unknown_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r2_record[1]+ "\n"  + r2_record[2]+ "\n"  +  r2_record[3]+"\n")
            unknown_count_dict["unknown"] += 1
            ### If it doesn't meet our quality score, then we don't know what it is so we put into unknown

### at this point we have checked for all N's, not in index,  does not meet our quality score

        elif i1_record[1] in index_set:
            ### if there is a match to our set of indexes then we move ahead

            if reverse_comp_dict[i1_record[1]] == i2_record[1]:
                if bioinfo.qual_score(i2_record[3]) >= ic: ### index score greater than given threshold

                    if bioinfo.qual_score(r2_record[3]) >= rc: ### read score is greater than given threshold
                        ### If the file is good and meets all of our thresholds then we save it
                        ### these are the good files
                        # print(i1_record[1])
                        # print(i2_record[1])
                        open_file[i1_record[1]][0].write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                        open_file[i1_record[1]][1].write(r2_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r2_record[1]+ "\n"  + r2_record[2]+ "\n"  +  r2_record[3]+"\n")
                        match_counts[i1_record[1]] += 1

                    else:  
                        # print(i1_record[1])
                        # print(i2_record[1])
                        unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                        unknown_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                        unknown_count_dict["unknown"] += 1
                        ### these files didn't meet our threshold for our read cutoff and will be put into unknown
                        
                else: 
                    # print(i1_record[1])
                    # print(i2_record[1])
                    unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                    unknown_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                    unknown_count_dict["unknown"] += 1    
                    ### these file didn't meet the minimun requirement for index cutoff and are put into unknown             
            
            elif reverse_comp_dict[i1_record[1]] != i2_record[1]:
                ### these are files that don't have a matched index 

                if i2_record[1] in reverse_set:
                    ### if the complementary index is in another set implying index hopping

                    if bioinfo.qual_score(i2_record[3]) >= ic:      
                        ### index score threshold met

                        if bioinfo.qual_score(r2_record[3]) >= rc:   
                            # print(i1_record[1])
                            # print(i2_record[1])
                            hopped_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                            hopped_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                            hopping_count[i1_record[1]+"+"+reverse_complementary(i2_record[1])] += 1
                            ### the other index that is also present in our set, implying index hopping, met the critera
                            ### and is added to the hopped    
                        else:
                            # print(i1_record[1])
                            # print(i2_record[1])
                            unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                            unknown_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                            unknown_count_dict["unknown"] += 1
                            ### put into unknown as it didn't meet read cutoff threshold

                    else:
                        # print(i1_record[1])
                        # print(i2_record[1])
                        unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                        unknown_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                        unknown_count_dict["unknown"] += 1    
                        ### If the index threshold isn't met then we add it to the unknown                    

                else: 
                    # print(i1_record[1])
                    # print(i2_record[1])
                    unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                    unknown_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
                    unknown_count_dict["unknown"] += 1
                    ### if it not in the reverse set, then we can put it into the unknown because we don't 
                    ### know where it came from
                    
        else:
            # print(i1_record[1])
            # print(i2_record[1])
            unknown_r1.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
            unknown_r2.write(r1_record[0]+"_"+i1_record[1]+"+"+reverse_complementary(i2_record[1])+"\n"  + r1_record[1]+ "\n"  + r1_record[2]+ "\n"  +  r1_record[3]+"\n")
            print("WTF")
            unknown_count_dict["unknown"] += 1
            ### everything else we toss into the unknown
            
### close files

# print(unknown_count_dict)
# print(match_counts)
# print(hopping_count)
######################################################################################################################################
### we are going to zip up all our files

unknown_r1.close()
unknown_r2.close()

hopped_r1.close()
hopped_r2.close()

for i in index_set:
    open_file[i][0].close
    open_file[i][1].close

######################################################################################################################################
### in here we are going to output a file with all the stats of our fq reads
### we want percentage of read from each sample
### overall amount of index swapping
total_hopped = sum(hopping_count.values())
total_matched = sum(match_counts.values())
total_unknown = sum(unknown_count_dict.values())
total= total_hopped + total_matched + total_unknown

match_counts_sorted = sorted(match_counts, key = match_counts.get, reverse = True)
hopping_count_sorted = sorted(hopping_count, key = hopping_count.get, reverse = True)

with open("read_statistics.md","w") as stats:
    stats.write("Input file 1"+" = "+file1+"\n")
    stats.write("Input file 2"+" = "+file2+"\n")
    stats.write("Input file 3"+" = "+file3+"\n")
    stats.write("Input file 4"+" = "+file4+"\n")
    stats.write("Read cutoff"+" = "+str(rc)+"\n")
    stats.write("Index cutoff"+" = "+str(ic)+"\n")

    stats.write("Total Index Hopping =" + " " + str(total_hopped) + "\n")
    stats.write("Total Matched =" + " "+str(total_matched) + "\n")
    stats.write("Total Number Of Reads= " + str(total)+ "\n")
    stats.write("Index" + "\t" + "Count" + "\t" + "Percentage" + "\n")

    # for n in index_set:
    #         stats.write(n+"\t"+str(match_counts[n])+"\t"+str(match_counts[n]/total*100)+"%"+"\n")

    # stats.write("Unknown"+"\t"+str(total_unknown)+"\t"+str(total_unknown/total*100)+"%"+"\n")

    # for i in permutation_set:
    #     stats.write(i+"\t"+str(hopping_count[i])+"\t"+str(hopping_count[i]/total*100)+"%"+"\n")

    # stats.write("Total Index Hopping ="+" "+str(total_hopped)+"\n")
    # stats.write("Total Matched ="+" "+str(total_matched)+"\n")
############################################################# TEST
    for n in match_counts_sorted:
            stats.write(n+"\t"+str(match_counts[n])+"\t"+str(match_counts[n]/total*100)+"%"+"\n")

    stats.write("Unknown"+"\t"+str(total_unknown)+"\t"+str(total_unknown/total*100)+"%"+"\n")

    for i in hopping_count_sorted:
        stats.write(i+"\t"+str(hopping_count[i])+"\t"+str(hopping_count[i]/total*100)+"%"+"\n")

    stats.write("Total Index Hopping ="+" "+str(total_hopped)+"\n")
    stats.write("Total Matched ="+" "+str(total_matched)+"\n")


