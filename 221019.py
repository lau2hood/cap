#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from thefuzz import fuzz, process
import pandas as pd
import editdistance


# In[ ]:


#This function opens the fastq file that you want to read
def openfile():
    name = input("File name: ")
    file = open(name,'r')
    fastq = file.readlines()
    return fastq

#This function opens the file containing the region of interest
def openregion():
    test = input("File name for region of interest: ")
    file2 = open(test,'r')
    needle = file2.readlines()
    return needle

#This function makes a list of the sequence names
def make_seq(old):
    i = 0
    seq = []

    while i < len(old):
        #seq.append(fastq[i])
        seq.append(old[i].rstrip('\n'))
        i += 4
    return seq

#This function makes a list of the DNA sequences
def make_read(old):
    i = 1
    read = []

    while i < len(old):
        #read.append(fastq[i+1])
        read.append(old[i].rstrip('\n'))
        i += 4
    return read

#This function makes a list of the sequence we are searching for
def make_needle(old):
    i = 0
    sub = []

    while i < len(old):
        #read.append(fastq[i+1])
        sub.append(old[i].rstrip('\n'))
        i += 1
    return sub

#This function splits the sequence we are looking for into chunks of a specified number of bases
def extract(mylist):
    mystring = ''
    for x in mylist:
        mystring += ''+x
    
    length = int(input("How many bases should I look for? "))
    
    i = 0 
    j = i + length
    regions = []
    while (i < len(mystring)):
        region = mystring[i:j]
        regions.append(region)
        i=j
        j = i + length
    return regions

#This function takes the sequences and the chunks that we are searching for and returns the sequences that contain what we are looking for
def string(strings, sub):
    i = 0
    while (i < len(sub)):
        substring = sub[i]
        strings_with_substring = [string for string in strings if substring in string]
        i += 1
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
            index_pos_list.append(index_pos)
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
def access(oldlist):
    i = 0
    new = []
    while i < len(oldlist):
        int1 = oldlist[i][0]
        new.append(int1)
        i += 1
    return new

#This function takes the indecies from the access function and finds the sequence IDs at those indecies 
import numpy as np
def findElements(lst1, lst2):
    return [lst1[i] for i in lst2]

#This function creates a new file with the sequence IDs (on new lines) found in the findElements function 
def newfile(sequence_ids):
    idname = input("New file name: ")
    with open(idname, 'w') as f:
        for line in sequence_ids:
            f.write(f"{line}\n")


# In[ ]:


def distance(strings, sub, allow):
    j = 0
    lists = []
    while (j < len(strings)):
        difference = len(strings[j]) - len(sub)
        distance = editdistance.eval(strings[j], sub)
        errors = distance - difference
        if (errors <= allow):
            lists.append(strings[j])
    j+=1
    return lists
    
def removeDups(sequences):
    res = [*set(sequences)]


# In[ ]:


def main():
    #Open fastq file
    read_file = openfile()
    
    #Open file containing region of interest
    read_region = openregion()
    
    #make a list of read IDs
    identification = make_seq(read_file)
    
    #make a list of reads
    reads = make_read(read_file)
    
    #make a list of the region of interest
    region = make_needle(read_region)
    
    #create a list of the region of interest divided into specified sections
    region_list = extract(region)
    
    #create loop to run through the list of reads to search for any that contain the various sections of the region of interest
    i = 0
    allowed = int(input("How many errors to allow? "))
    while (i < len(region_list)):
        #find the reads that have the region of interest in the i-th position of the list of sectioned region of interest
        #read_list = string(reads, region_list[i])
        distance_list = distance(reads, region_list[i], allowed)
        removed = removeDups(distance_list)
        
        
        #create a nested list of the indicies of the reads containing the region of interest
        ilist = index_list(read_file, removed) #read_list
        
        #un-nest the list of indicies
        matches = access(ilist)
        
        #find the sequence IDs at the indicies 
        find = findElements(identification, matches)
        
        #move forward one position
        i += 1
    #create a new file with the list of sequences
    sequence_list = newfile(find)
    print(read_list, "\n\n")
    print(removed)


# In[ ]:


main()


# In[ ]:





# In[ ]:




