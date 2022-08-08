# Author: <Peter Pham> <peterpham.ghs@email.address>


'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
More can be added here to help with anything you need'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = "GCTAgcta"
RNA_bases = "GCAUgcau"

def convert_phred(letter:str)->int:
    '''docstring'''

    return ord(letter)-33

    # if __name__ == "__main__":
    # assert convert_phred("I") == 40, "wrong phred score for 'I'"
    # assert convert_phred("C") == 34, "wrong phred score for 'C'"
    # assert convert_phred("2") == 17, "wrong phred score for '2'"
    # assert convert_phred("@") == 31, "wrong phred score for '@'"
    # assert convert_phred("$") == 3, "wrong phred score for '$'"
    # print("Your convert_phred function is working! Nice job")
    # pass


### this then continues over the the next line of code.

def qual_score(phred_score):
    '''docstring'''
    total_length=len(phred_score) ##length of the phred score
    total_score=0
    
    for letter in phred_score:  ##individual characters in the phred score
        
        phred_score=convert_phred(letter)  ## gets their score
        
        total_score+=phred_score  ##sum of the whole line
        
    return total_score/total_length  ## gets you your average






def validate_base_seq(seq,RNAflag:bool=False):
    '''Checks to see if sequence is a sequence'''
    seq=seq.upper()
    total_sequence=seq.count("A")+seq.count("G")+seq.count("C")+seq.count("U" if RNAflag else "T")
    return seq==total_sequence

    # if __name__ == "__main__":
    # assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    # assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    # assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    # assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    # print("Passed DNA and RNA tests")
    # pass

def gc_content(dna):

    assert validate_base_seq(dna)

    dna = dna.upper()
    g = dna.count("G")
    c = dna.count("C")

    return ((g+c)/len(dna))


# def oneline_fasta(fasta_file,new_fasta):
#     '''docstring'''
#     with open(fasta_file,"r") as old_fa, open(new_fasta,"w") as new_fa:
#         ### old_fa has sequences with multiple lines and not one
#         ### our new file will have one sequence line to make it easier for downstream applications

#         for line in old_fa:
#             if line[0] != ">":
#                 seq = line.strip("\n")
#                 ### this grabs the sequence line and ignores the header 
#                 ### in this case the line with the >
            
#             else:
#                 if seq !="":
#                     new_fa.write(seq + "\n")
#                 seq = ""
#                 new_fa.write(line)
#                 ###writes out the first > to indicate the read
#         else:
#             new_fa.write(seq +"\n")
#             ### writes otu the new sequence with new line



# if __name__ == "__main__":
#     # write tests for functions above
#     pass