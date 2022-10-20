#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from thefuzz import fuzz, process
import pandas as pd
import editdistance
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')
from matplotlib.animation import FuncAnimation
from IPython.display import display, clear_output
import time
import gzip


def param():
    f = input('File containing parameters: ')
    file = open(f, 'r')
    return file.read()


def readReads(parameters):
    reads = parameters[0]
    readFile = open(reads,'r')
    return readFile.read(), 


def readInterest(parameters):
    interest = parameters[1]
    readInterest = open(interest,'r')
    return readInterest.read()


def makeSeq(readFile):
    IDs = []
    for i in range(0, len(reads), 4):
        IDs.append(readFile[i].rstrip())
    return IDs


def makeReads(readFile):
    readList = []
    for i in range(1, len(reads), 4):
        readList.append(readFile[i].rstrip())
    return readList


def makePlus(readFile):
    plus = []
    for i in range(2, len(reads), 4):
        plus.append(readFile[i].rstrip())
    return plus


def makeQuality(readFile):
    quality = []
    for i in range(3, len(reads), 4):
        quality.append(readFile[i].rstrip())
    return quality


def makeInterest(interest):
    region = []
    for i in range(1, len(interest)):
        region.append(interest[i].rstrip())
    return region


def extractParts(interest, parameters):
    length = int(parameters[2])
    parts = [] 
    i = 0
    j = length
    while i < len(interest):
        section = interest[i:j]
        if (len(section) != length):
            section = interest[(-length):]
            parts.append(section)
        else:
            parts.append(section)
        i += 1
        j += 1
    return parts


def stringMatching(readList, parts):
    goal = int(parameters[3])
    stringMatches = [] 
    for i in range(0, len(parts)):
        substring = parts[i]
        for i in range(0, len(readList)):
            score = fuzz.partial_ratio(substring, reads[i])
            if (score >= goal):
                stringMatches.append(readList[i])
    return stringMatches


def indexPositions(readList, stringMatches):
    positions = []
    i = 0
    while True:
        try:
            indexPosition = readList.index(stringMatches, i)
            positions.extend(i)
            i += 1
        except ValueError as e:
            break
    return positions


def findElements(list1, positions):
    matches = []
    for i in positions:
        matches.append(list1[i])
    return matches

def newFasta(parameters, inputs):
    fileName = parameters[4] + '.fasta'
    with open(fileName, 'w') as file:
        for line in inputs:
            file.write(f"{line}\n")
    return file

            
def newFastq(parameters, inputs):
    fileName = parameters[4] + '.fastq'
    with open(fileName, 'w') as file:
        for line in inputs:
            file.write(f"{line}\n")
    return file


def main():
    parameters = param()
    
    readFile = readReads(parameters)
    ids = makeSeq(readFile)
    readList = makeReads(readFile)
    plus = makePlus(readFile)
    quality = makeQuality(readFile)
    
    readRegion = readInterest(parameters)
    interest = makeInterest(readRegion)
    
    parts = extractParts(interest, parameters)
    
    foundFasta = []
    foundFastq = []
    for i in range(0, len(parts)):
    
        matching = stringMatching(readList, parts)
    
        positions = indexPositions(readList, stringMatching)
    
        sequenceIDs = findElements(ids, positions)
        sequenceReads = findElements(readList, positions)
        sequencePlus = findElements(plus, positions)
        sequenceQuality = findElements(quality, positions)
    
        for i in range(0, len(sequenceIDs)):
            if sequenceIDs[i] not in foundFasta:
                foundFasta.append(sequenceIDs[i])
                foundFasta.append(sequenceReads[i])
                
                foundFastq.append(sequenceIDs[i])
                foundFastq.append(sequenceReads[i])
                foundFastq.append(sequencePlus[i])
                foundFastq.append(sequenceQuality[i])

    fasta = newFasta(parameters, foundFasta)
    fastq = newFastq(parameters, foundFastq)


# In[ ]:


main()

