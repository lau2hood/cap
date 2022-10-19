#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from thefuzz import fuzz, process
import pandas as pd
import editdistance
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')
from matplotlib.animation import FuncAnimation
from IPython.display import display, clear_output
import time
import gzip


# In[ ]:


#This function opens the parameters file
def param():
    name = input("File name containing parameters: ")
    file = open(name, 'r')
    contents = file.readlines()
    return contents

#This function opens the fastq file that you want to read
def openfile(parameters):
    name = parameters[0].rstrip() #input("File name (including fastq.gz): ")
    file = open(name, 'r')
    file_contents = file.readlines()
    return file_contents

#This function opens the file containing the region of interest
def openregion(parameters):
    test = parameters[1].rstrip() #input("File name for region of interest: ")
    file2 = open(test,'r')
    needle = file2.readlines()
    return needle

#This function makes a list of the sequence names
def make_seq(old):
    i = 0
    seq = []

    while i < len(old):
        #seq.append(fastq[i])
        seq.append(old[i].rstrip())
        i += 4
    return seq

#This function makes a list of the DNA sequences
def make_read(old):
    i = 1
    read = []

    while i < len(old):
        #read.append(fastq[i+1])
        read.append(old[i].rstrip())
        i += 4
    return read

#This function makes a list of the pkus signs
def make_plus(old):
    i = 2
    plus = []

    while i < len(old):
        #read.append(fastq[i+1])
        plus.append(old[i].rstrip())
        i += 4
    return plus

#This function makes a list of the quality scores
def make_quality(old):
    i = 3
    quality = []

    while i < len(old):
        #read.append(fastq[i+1])
        quality.append(old[i].rstrip())
        i += 4
    return quality

#This function makes a list of the DNA sequences
def make_region(old):
    i = 1
    region = []

    while i < len(old):
        #read.append(fastq[i+1])
        region.append(old[i].rstrip())
        i += 2
    return region

#This function splits the sequence we are looking for into chunks of a specified number of bases
def extract(mylist, parameters):
    mystring = ''
    for x in mylist:
        mystring += ''+x
    
    length = int(parameters[2]) #int(input("How many bases should I look for? "))
    regions = []
    for i in range(len(mystring)):
        region = mystring[i: (i + length)]
        if (len(region) != length):
            region = mystring [(-length):]
            #print('last region: ', region)
            regions.append(region)
        else:
            regions.append(region)
    return regions

#This function takes the sequences and the chunks that we are searching for and returns the sequences that contain what we are looking for
def string(strings, sub, goal):
    strings_with_substring = []
    for i in range(len(sub)):
        #print("score iteration #", i)
        substring = sub[i]
        for i in range(len(strings)):
            #if substring in string:
            score = fuzz.partial_ratio(substring, strings[i])
            #print(substring)
            #print(score)
            if (score >= goal):
                strings_with_substring.append(strings[i])
    return strings_with_substring

#This function returns the indexes of all occurrences of given element in the list
def get_index_positions(list_of_reads, element):
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_reads.index(element, index_pos)
            # Add the index position in list
            index_pos_list.extend(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list

#This function returns a list of the indecies where the chunk we are looking for occurs, returns nested list
def index_list(filename, seqlist):
    indexlist = []
    for i in range(0, len(seqlist)):
        pos = get_index_positions(make_read(filename), seqlist[i])
        indexlist.append(pos)
    return indexlist

#This function un-nests the list output by index_list function
#def access(oldlist):
#    new = []
#    for i in range(len(oldlist)):
#        int1 = oldlist[i][0]
#        new.append(int1)
#    return new

#This function takes the indecies from the access function and finds the sequence IDs at those indecies 
import numpy as np
def findElements(lst1, lst2):
    return [lst1[i] for i in lst2]

#This function creates a new file with the sequence IDs (on new lines) found in the findElements function 
def newfasta(sequence_ids, idname):
    idname = idname + '.fasta'
    with open(idname, 'w') as f:
        for line in sequence_ids:
            f.write(f"{line}\n")
            
#This function creates a new file with the sequence IDs (on new lines) found in the findElements function 
def newfastq(sequence_ids, idname):
    idname = idname + '.fastq'
    with open(idname, 'w') as f:
        for line in sequence_ids:
            f.write(f"{line}\n")


# In[ ]:


def main():
    #Open parameters file
    parameters = param()

    #Open fastq file
    read_file = openfile(parameters)
    print("The file has been read. Currently creating list of IDs.")
    
    #make a list of read IDs
    identification = str(make_seq(read_file))
    print("The list of IDs has been created. Currently creating list of sequences.")
    
    #make a list of reads
    reads = str(make_read(read_file))
    print("The list of sequences has been created. Currently creating list of plus signs.")
    
    #make a list of the plus signs
    plus = str(make_plus(read_file))
    print("The list of plus signs has been created. Currently creating list of quality scores.")
    
    #make a list of the quality scores
    qual = str(make_quality(read_file))
    print("The list of quality scores has been created. Currently reading file.")
    
    #Open file containing region of interest
    read_region = openregion(parameters)
    print("The file has been read. Currently creating region of interest.")
    
    region = str(make_region(read_region))
    print("The region has been created.")
    
    #create a list of the region of interest divided into specified sections
    region_list = extract(region, parameters)
    print('Region list has been created. Moving into loop.')
    
    #create loop to run through the list of reads to search for any that contain the various sections of the region of interest
    found_fasta = []
    found_fastq = []
    print("Empty lists created.")

    #allowed = int(input("How many errors to allow? "))
    percent_match = int(parameters[3]) #int(input("% match goal: "))
    print("Percent found. Moving into loop 2 .")
    for i in range(len(region_list)):
        #find the reads that have the region of interest in the i-th position of the list of sectioned region of interest
        #region_list[i] = [region_list[i]]
        
        read_list = string(reads, region_list[i], percent_match)
        print("Read list made.")
        
        

        #create a nested list of the indicies of the reads containing the region of interest
        matches = index_list(read_file, read_list) #results_list)
        print("Matches list created.")
        
        #matches = access(ilist)
        
        #find the sequence IDs at the indicies 
        findID = findElements(identification, matches)
        findseq = findElements(reads, matches)
        findplus = findElements(plus, matches)
        findquality = findElements(qual, matches)
        #print('find', find, '\n\n')
        print("List of sequences created.")
        
        for j in range(len(findID)):
            if findID[j] not in found:
                found_fasta.append(findID[j])
                found_fastq.append(findID[j])
                
                found_fasta.append(findseq[j])
                found_fastq.append(findseq[j])
                
                found_fastq.append(findplus[j])
                
                found_fastq.append(findquality[j])
 
 
    #create a new file with the list of sequences
    new_name = parameters[4] #input('New file name: ')
    sequence_list = newfasta(found_fasta, new_name)
    fastq_list = newfastq(found_fastq, new_name)

The required inputs for the code are:
1. The name of the file (including file type - .txt, .fastq, etc.)

2. The name of a file containing the seqeunce/region of interest (including file type - .txt, .fastq, etc.)
    a. You may need to create this file. It can be a txt file created within the Anaconda application

3. The number of bases you want to search for at a time (e.g. if you select 20 bp for a 200 bp region of interest, your region would be split into ten 20 bp long "chunks")

4. The goal for how much of the sequence must match the read. This is a little fickle, so we may need to play around with the optimal number. It may also be dependent on the number of bases in the read - if you have 10bp, you may want a higher value vs 100bp)

5. New file name. The program will create a fasta-like file containing the Read IDs and reads that match the region of interest. You will give it a unique name and it will be saved in the same folder as this Notebook file. You can then download it locally. 
# In[ ]:


main()

The commented out section of the main() function creates a graph showing the number of matches foud at each iteration. It makes the program take a very long time to run, but it can be a cool visual if you don't mind the extra-loong run time.

The normal program (without graphics) takes approximately 5-10 minutes to complete for a 24 kB file, though that time varies drastically depending on the file size.  