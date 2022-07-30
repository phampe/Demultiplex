# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1  | 363246735 | -33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 363246735 | -33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 363246735 | -33 |
| 1294_S1_L008_R4_001.fastq.gz | read2  | 363246735 | -33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.


    2. We care about our indexes more than our biological data. There are less indexes which means less room for error. If there is more room
       for error when we are working with our read lines because we have more of it. As long as we can align the read to something then we are
       going to be fine in our downstream applications. In relationship to this data for our indexes I would set the lower limit to a qscore of 35. 
       With the best score being 40, we have some averages of 30 which based on the distrubition of are averages are bad. 
       With relationship to our read scores as there are a ton of scores that are nearing 40 I would say we could accept any of the data that is above 
       the averages at their respective positions. These are less of an issue to have a bad score on as we are just lining them up with our genomic template,
       but I do still want them to be above average scores.


    3. zcat 1294_S1_L008_R2_001.fastq.gz |grep -A1 "^@"| grep -v
       "^--"| grep -v "^@"|grep "N"|wc -l
       this gives use a count of 3976613 for index 1

       zcat 1294_S1_L008_R3_001.fastq.gz |grep -A1 "^@"| grep -v
       "^--"| grep -v "^@"|grep "N"|wc -l
       This will give us a count of 3328051 for index 2
    
## Part 2
1. Define the problem

   We have four files. Each file has a property of the read that we want to extrapulate and read them in parallel
   There are going to be issues with the reads and we need to separate them accordingly.
   There will be good reads that match with their proper indexes, those that have missmatches that are caused by index hopping
   and finally those that we don't know where they came from as unknown index.
   Finally we need to cut off those quality scores that fall below a certain threshold and then output this all into 4 different files

2. Describe output

   There needs to be four files that are output from our 4 files. 
   The first file is read_1 that have been accepted along with a file corresponding as read_2 that is the other end of the paired read
   two more files, one for read_1 that is unmatched and one for read_2 that is unmatched
   two more files on top of this so that we see those that are unknown or low quality

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

4. Pseudocode

### There is more written on the actual psuedo_code.txt file
###
Make an argparse so we can select our files and cutoff for our qscore

input files
   read 1
   read 2
   index 1 
   index 2

output files
   matched read_1
   matched read_2

   unmatched read_1
   unmatched read_2

   unknown read_1
   unknown read_2

list of sequences from file given by leslie
list of comp

set up a set for our complementary sequence

dictionary of sequences as keys and and reverse complements as values

Dictionary of counts for each index matched, unmatched and unknown

First we would need to open all the files so that they are running in parallel as well as open all the files we are going to output to
and set them to own variable

set up dictionaries of all teh read of individual indexes 48 indexes.
   read 1 gets one index while read 2 gets tehe second index.


### everything below here will be looped through

We now need to read index_1 and match it with our list of known indexes

   First we are going to go through the the sequence rather than the individual bases

   if there is an index in index_1 that does not exist we put it in the unknown
      grab the corresponding index_2 and concatenate 
      add the index to the header of our read_1 and read_2, then put them into unknown read_1 and read_2
         put a counter on the unknown

### for this we are going to go through each of the nucleotide for index_1

### leslie said to put this higher up check the n's first
   if there is a N for nucleotide we immeidately put the files into unknown
      grab the corresponding index_2 and concatenate
      add the index to the header of our read_1 and read_2, then put them into unknown read_1 and read_2
         put a counter on the unknown

   if the quality score is too low then we skip and put in the unknown
      grab the corresponding index_2 and concatenate
      add the index to the header of our read_1 and read_2, then put them into unknown read_1 and read_2
         put counter on the unknown

### for this we are looking through the sequence for index_2

   if it matches what we have in our list we go onto the index_2

      if the index in index_2 is the correlated with the value in the dictionary for index_1
         now we check to see if it meets the threshold for our quality scores
      
            if it meets our standards then we concatenate the files and add the index pair to the header
               add the file to our matched
               add counter for matched

            if it doesn't then we put it into our unknown
               add the file to our unknown
               add counter for unknown
      
      if it the two indexes don't match
         check to see if it is a possible complementary from another index
            it is is we check the quality score of it
               if it is over the threshold 
                  add files to unmatched
                  add counter to unmatched
                  
               if it doesn't reach our threshold
                  add files to unknown
                  add counter to unknown

      anything else
         we add it to files to unknown
         add counter to unknown
      
Close all our files at the end

iter-tools dictionary

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
