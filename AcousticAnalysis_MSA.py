import os
import shutil
import csv
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.interpolate import make_interp_spline, BSpline

print('Input path name. Empty response uses current working directory.\n')
workingDirectory = input('Folder/path of file to analyze: ') or os.getcwd()

file = input('Name of .csv file to analyze: ')

csvfile = os.path.join(workingDirectory, file)

f = open(csvfile) #can use input later to ask for file name
csv_f = csv.reader(f, delimiter= ",")
csv_f = list(csv_f)

dataOutFolder = file[:-4] + 'MSAUSVAnalysis'
newDataFolder = os.path.join(workingDirectory, 'MSAUSVAnalysis')
if not os.path.exists(newDataFolder):
    os.makedirs(newDataFolder)
newnewDataFolder = os.path.join(newDataFolder, file[:-4])
if not os.path.exists(newnewDataFolder):
    os.makedirs(newnewDataFolder)

batchName = []
syllableType = []
harmonic = []
duration = []
isi = []
fqVariance = []
spectralPurity = []
amplitude = []
fqMin = []
fqMean = []
fqMax = []
fqStart = []
fqEnd = []
totalBandwith = []

##############################################################################################
isiCutoff = .25 #can use input later to ask for specific isicutoff
recordLen = 5 #5 minutes per recording, can use input later
##############################################################################################
for row in range(1,len(csv_f)):
    batchName.append(csv_f[row][0])
    syllableType.append(csv_f[row][1])
    harmonic.append(float(csv_f[row][2]))
    duration.append(float(csv_f[row][3]))
    isi.append(float(csv_f[row][4]))
    fqVariance.append(float(csv_f[row][5]))
    spectralPurity.append(float(csv_f[row][6]))
    amplitude.append(float(csv_f[row][7]))
    fqMin.append(float(csv_f[row][8]))
    fqMean.append(float(csv_f[row][9]))
    fqMax.append(float(csv_f[row][10]))
    fqStart.append(float(csv_f[row][11]))
    fqEnd.append(float(csv_f[row][12]))
    totalBandwith.append(float(csv_f[row][13]))

print('There are ' + str(len(csv_f)) +' USVs in this csv file')

usvVars = {}
usvVars['Batch Name'] = batchName
usvVars['Syllable Type'] = syllableType
usvVars['Harmonic'] = harmonic
usvVars['Duration'] = duration
usvVars['ISI'] = isi
usvVars['fq Variance'] = fqVariance
usvVars['Spectral Purity'] = spectralPurity
usvVars['Amplitude'] = amplitude
usvVars['fq Min'] = fqMin
usvVars['fq Mean'] = fqMean
usvVars['fq Max'] = fqMax
usvVars['fq Start'] = fqStart
usvVars['fq End'] = fqEnd
usvVars['Total Bandwith'] = totalBandwith


listTot = [batchName, syllableType, harmonic, duration, isi, fqVariance, spectralPurity, amplitude, fqMin, fqMean, fqMax, fqStart, fqEnd, totalBandwith]

#Removing the Unclassified
for x in range(len(syllableType)):
    if syllableType[x] == 'notIDd':
        isi[x-1] = isi[x-1] + isi[x]

indices = [i for i, x in enumerate(syllableType) if x == 'notIDd']

for smallList in listTot:
    for index in range(len(indices)):
        smallList.pop(indices[index]-index)

#####################################################################################
''' DATA'''

#syl
#Based on the number of pitch jumps in a single syllable
#Categorization of the 4 syllable types
syl = []
for x in syllableType: #goes through each of the syllable types in the orginal list
    if len(x) >= 2: #want to know how many pitch jumps there were in a syllable
        syl.append('M') #if more than one - labelled as M, add to new list we made
    else:
        syl.append(x) #if one pitch/pitch jump - keep it labelled the same 

'''NEED TO FIX CODE! WHAT IF THE FIRST SYLLABLE = STOP?'''    

#sequence1
#Based on the ISI cutoff value - should be 250 ms 
#anything equal to or greater than ISI cutoff = STOP 
#will start a new sequence for the next one
seq1 = []
for x in isi: #goes through the isi values of each syllable
    if x <= isiCutoff: #syllable is part of the same sequence if isi value less than the cutoff
        seq1.append('1') #labelled as 1 in the new list that we made
    else: #means that isi is higher than the cutoff - we stop the sequence here
        seq1.append('STOP') #labelled as STOP since we stop the sequence here
        
        
#sequence2
#just keeps count of syllables 
#based on where the syllable is within the sequence
seq2 = []
string = ''
for stopOrNot in range(len(seq1)): #going through the 1's and stops in seq ; cant also just use usvData[4] and cutoff if the functions need to run together
    if stopOrNot == 0 and seq1[stopOrNot] == 'STOP':
        seq2.append('STOP')
    elif seq1[stopOrNot] =='1': #if x equals 1 from seq1
        string = string + 'a' #add an a to the string
        seq2.append(string) #add that string to the list
    else: #if not 1, must be stop
        string = string + 'a' #add final a to the string, end of seq
        seq2.append(string) #add that string to the list
        string = '' #reset the string to an empty string for the new sequence starting with the next value


#length
#telling you the length of the sequences (how many syllables are in them)
length = []
for x in range(len(seq1)): #going through the values in seq 1 to see where the stops are
    if x == 0 and seq1[0] == 'STOP': length.append(1)
    elif seq1[x] == "STOP": #if that val is a stop
        length.append(len(str(seq2[x]))) #want to get the length of that sequence, add that number to the list
    else: #if not a stop
        length.append(' ') #then just add a blank space to the list to keep the order ok


#1st call
#determines the first syllable type is for the start of a sequence 
firstCall = []
for x in range(len(seq1)): #goes through seq1 to see where the stops are since new seq starts right after
    if x == 0 and seq2[0] == 'STOP': #first value is always a new sequence 
        firstCall.append('') #add sylalble type to the list
    else: #for the rest of the values, we can look back to see if there was a stop right before
        if seq1[x-1] == 'STOP': #if the last value was stop
            firstCall.append(syl[x]) #add the type of syllable to the list
        else: #if the last one was not a stop aka this syllable is not the start of a sequence
            firstCall.append(' ') #just adds a blank space to the list to keep things in order


#trans/sequence
#specifies current and next syllable type
#if last, use x for the next syllable type (silence)
transSeq = []
for x in range(len(seq1)): #goes through seq1 to see where sequences stop
    if x == (len(seq1)-1): #first checks to see if at the last syllable 
        transSeq.append(syl[x] + 'X') #just add its own syllable type to the list
    else: #if not the last syl
        if seq1[x] == 'STOP': #if the current syl is a stop aka there is silence after
            transSeq.append(syl[x] + 'X') #add a string of the current syl type and X (silence) to the list
        else: #if not a stop syl
            transSeq.append(syl[x] + syl[x+1]) #add a string of current and next syl type to the list

#trans from X
#telling us at the last syllable of a sequence and what the first syllable of the next sequence is
#except at end of file 
transFromX = []
for x in range(len(seq1)): #goes through the seq1 vals to see where stops are
    if seq1[x] == 'STOP': #if seq stops at that syl
    	if x != (len(seq1)-1):
    		transFromX.append('X' + syl[x+1]) #add a string of X (silence) and next syl type to the list
    	else:
    		transFromX.append('X') #add a string of X (silence) and next syl type to the list
    else: #if not at a stop syl
        transFromX.append(' ') #add empty space to list to keep things in order 


#Number of C
#at a specific instance within the sequence it is telling you how many syllables have been categorized as M so far
#The last one in a sequence should say "SEQ OK" if there's more than 1 (>=2) type M syllables
numC = []
num = 0
for x in range(len(seq1)): #going through seq1 for stops
    if seq1[x] == 'STOP': #checks for stop syl first
        if num >= 2: #AND there's more than one M in the whole sequence
            numC.append('SEQ OK') #label as SEQ OK in list
        else: #if syl is stop and there's <2 M's in the whole sequence
            numC.append(num) #label as current number of M's (can only be 0 or 1)
        num = 0  #reset num value to 0 for new sequence
    else: #if not stop syl
        if syl[x] == 'M': #if that syl type is M
            num = num + 1 #num goes up by one 
        numC.append(num) #if none of the above, then just append current num; want to know how many M's at a specific point in a seq

leifEriksonDay = []
for x in range(len(syl)):
	newl = [syl[x], seq1[x], seq2[x], length[x], firstCall[x], transSeq[x], transFromX[x], numC[x]]
	leifEriksonDay.append(newl)



colNames = ['syl', 'seq1', 'seq2', 'length', 'firstCall', 'trans/seq', 'transFromX', 'numC']
cols = [syl, seq1, seq2, length, firstCall, transSeq, transFromX, numC]


with open(newnewDataFolder + '/BasicTransofmrations.csv', 'w', newline = '') as results:
	writer = csv.writer(results)
	writer.writerow(colNames)
	for hoopla in leifEriksonDay:
		writer.writerow(hoopla)

#####################################################################################
'''FEATURES'''

#unique function for lists
batchNew = []
batchNew.append(batchName[0]) #add the first name to the list
for x in range(1,len(batchName)): #go through each of the names
    if batchName[x] != batchName[x-1]: #if the current name is not the same as the last, add to list
        batchNew.append(batchName[x])
        

#Nb USV
#total number of USVs for each individual mouse recording
#num syllables = num USVs
nbUSV = []
for x in batchNew: #for each id in the list of different ids
    usvTotal = 0 #total starting USV is 0
    for y in batchName: #for each id in the full list of ids
        if x == y: #if the id from the list of diff ids = id from full list
            usvTotal = usvTotal + 1 #add one to the number of USVs
    nbUSV.append(usvTotal) #at the end, append the total number of usvs for each id


#call rate
#nb USV/ length of recording (in minutes)
#catch - reference call rate to recording time (change if needed)
callRate = []
for x in nbUSV: #for each set of total usvs in a recording
    callRate.append(x/recordLen) #call rate = that total/recording length of each recording


#% starting syllable
#% of starting syllables in total sequences of each recording
#number of starting syllables = proxy for num of sequences
listS = [] #list for % of starting s in each recording
listD = [] #list for % of starting d in each recording
listU = [] #list for % of starting u in each recording
listM = [] #list for % of starting m in each recording
for ident in range(len(nbUSV)): #for each of the separate recordings
    s = 0 #number of starting s in a recording
    d = 0 #number of starting d in a recording
    u = 0 #number of starting u in a recording
    m = 0 #number of starting m in a recording
    t = 0 
    if ident == 0: #for the first recordings
        #bla = 0
        for firstSyl in range(0, nbUSV[ident]): #first --> x char in firstcall
            #bla = bla + 1
            if isinstance(length[firstSyl], str):
                if firstCall[firstSyl] == 's': #if char = s
                    s = s + 1 #add 1 to num of starting s in that recording
                    t = t + 1
                elif firstCall[firstSyl] == 'd': #if char = d
                    d = d + 1 #add 1 to num of starting d in that recording
                    t = t + 1
                elif firstCall[firstSyl] == 'u': #if char = u
                    u = u + 1 #add 1 to num of starting u in that recording
                    t = t + 1
                elif firstCall[firstSyl] == 'M': #if char = M
                    m = m + 1 #add 1 to num of starting m in that recording
                    t = t + 1
    else: #for all other recordings after the first
        for firstSyl in range((np.sum(nbUSV[0:ident])), (np.sum(nbUSV[0:ident+1]))): #for each char in firstcall from the end of the last recording to (last + current total usvs)
            if isinstance(length[firstSyl], str):
                if firstCall[firstSyl] == 's': #if char = s
                    s = s + 1 #add 1 to num of starting s in that recording
                    t = t + 1
                elif firstCall[firstSyl] == 'd': #if char = d
                    d = d + 1 #add 1 to num of starting d in that recording
                    t = t + 1
                elif firstCall[firstSyl] == 'u': #if char = u
                    u = u + 1 #add 1 to num of starting u in that recording
                    t = t + 1
                elif firstCall[firstSyl] == 'M': #if char = M
                    m = m + 1 #add 1 to num of starting M in that recording
                    t = t + 1
    if s > 0:
        listS.append(s/t*100) #add percentage of starting s for that recording to the list of s%
    elif s == 0:
       listS.append(0)
    if d > 0:
       listD.append(d/t*100) #add percentage of starting d for that recording to the list of d%
    elif d == 0:
       listD.append(0)
    if u > 0:
       listU.append(u/t*100) #add percentage of starting u for that recording to the list of u%
    elif u == 0:
       listU.append(0)
    if m > 0:   
        listM.append(m/t*100) #add percentage of starting M for that recording to the list of M%
    elif m == 0:
       listM.append(0)

#repertoire composition (%)
#all types/categories of USVs for a mouse
#% of all 4 types of syllables for each recording
listSt = [] #list for % of starting s in each recording
listDt = [] #list for % of starting d in each recording
listUt = [] #list for % of starting u in each recording
listMt = [] #list for % of starting m in each recording
for ident in range(len(nbUSV)): #for each of the separate recordings
    st = 0 #number of total s in a recording
    dt = 0 #number of total d in a recording
    ut = 0 #number of total u in a recording
    mt = 0 #number of total m in a recording
    tt = 0 #total number of syllables in a recording
    if ident == 0: #for the first recordings
        #bla = 0
        for currentSyl in syl[0:nbUSV[ident]]: #first --> x char in syl
            #bla = bla+1
            if currentSyl == 's': #if char = s
                st = st + 1 #add 1 to num of s in that recording
                tt = tt + 1 #add 1 to total num of syls
            elif currentSyl == 'd': #if char = d
                dt = dt + 1 #add 1 to num of d in that recording
                tt = tt + 1 #add 1 to total num of syls
            elif currentSyl == 'u': #if char = u
                ut = ut + 1 #add 1 to num of u in that recording
                tt = tt + 1 #add 1 to total num of syls
            elif currentSyl == 'M': #if char = M
                mt = mt + 1 #add 1 to num of m in that recording
                tt = tt + 1 #add 1 to total num of syls
    else: #for all other recordings after the first
        for currentSyl in syl[(np.sum(nbUSV[0:ident])):(np.sum(nbUSV[0:ident+1]))]: #for each char in syl from the end of the last recording to (last + current total usvs)
            if currentSyl == 's': #if char = s
                st = st + 1 #add 1 to num of s in that recording
                tt = tt + 1 #add 1 to total num of syls
            elif currentSyl == 'd': #if char = d
                dt = dt + 1 #add 1 to num of d in that recording
                tt = tt + 1 #add 1 to total num of syls
            elif currentSyl == 'u': #if char = u
                ut = ut + 1 #add 1 to num of u in that recording
                tt = tt + 1 #add 1 to total num of syls
            elif currentSyl == 'M': #if char = M
                mt = mt + 1 #add 1 to num of m in that recording
                tt = tt + 1 #add 1 to total num of syls
    listSt.append(st/tt*100) #add percentage of total s for that recording to the list of s%
    listDt.append(dt/tt*100) #add percentage of total d for that recording to the list of d%
    listUt.append(ut/tt*100) #add percentage of total u for that recording to the list of u%
    listMt.append(mt/tt*100) #add percentage of total M for that recording to the list of M%

#seq length
#1 value- average of all sequence lengths for each recording
seqLen = []
for x in range(len(nbUSV)): #for each of the separate recordings
    numSeq = 0 #number of sequences in that recording is 0 at start
    lenT = 0 #total length to average the sequences is 0 at start
    if x == 0: #for the first recording
        for y in length[0:nbUSV[x]]: #for first --> x char in length (first recording sequences end there)
            if y != ' ' and y > 1: #if y is not ' ' aka it is a number
                lenT = lenT + y #add that length to the total of the length
                numSeq = numSeq + 1 #add 1 to the number of sequences in that recording
    else: #for the recordings after the first one
        start = 0
        sCount = x-1
        while sCount >= 0:
            start = start + nbUSV[sCount]
            sCount = sCount - 1
        end = 0
        eCount = x
        while eCount >= 0:
            end = end + nbUSV[eCount]
            eCount = eCount - 1
        for y in length[start:end]: #for each char from last end point --> last + current num to get to next end index
            if y != ' ' and y > 1: #if y is not ' ' aka it is a number
                lenT = lenT + y #add that length to the total of the length
                numSeq = numSeq + 1 #add 1 to the number of sequences in that recording
    if numSeq>0:
        seqLen.append(lenT/numSeq) #divide the total length by the number of sequences in the recording, add that to the list
    elif numSeq==0:
         seqLen.append(0)

#fq Mean
#mean of fq means per recording
fqM = []
for x in range(len(nbUSV)): #for each of the separate recordings
    numFq = 0 #number of fq means in that recording is 0 at start
    meanT = 0 #total means to average the sequences is 0 at start
    if x == 0: #for the first recording
        for y in fqMean[0:nbUSV[x]]: #for first --> x char in length (first recording sequences end there)
            if y != ' ': #if y is not ' ' aka it is a number
                meanT = meanT + y #add that length to the total of the length
                numFq = numFq + 1 #add 1 to the number of sequences in that recording
    else: #for the recordings after the first one
        start = 0
        sCount = x-1
        while sCount >= 0:
            start = start + nbUSV[sCount]
            sCount = sCount - 1
        end = 0
        eCount = x
        while eCount >= 0:
            end = end + nbUSV[eCount]
            eCount = eCount - 1
        for y in fqMean[start:end]: #for each char from last end point --> last + current num to get to next end index
            if y != ' ': #if y is not ' ' aka it is a number
                meanT = meanT + y #add that length to the total of the length
                numFq = numFq + 1 #add 1 to the number of sequences in that recording
    fqM.append(meanT/numFq) #divide the total length by the number of sequences in the recording, add that to the list


#bandwidth
#average of total bandwith for each recording
bwM = []
for x in range(len(nbUSV)): #for each of the separate recordings
    numBw = 0 #number of total bandwidth in that recording is 0 at start
    bwT = 0 #total bandwidth to average the sequences is 0 at start
    if x == 0: #for the first recording
        for y in totalBandwith[0:nbUSV[x]]: #for first --> x char in length (first recording sequences end there)
            if y != ' ': #if y is not ' ' aka it is a number
                bwT = bwT + y #add that total bandwidth to the total of the length
                numBw = numBw + 1 #add 1 to the number of bandwidths in that recording
    else: #for the recordings after the first one
        start = 0
        sCount = x-1
        while sCount >= 0:
            start = start + nbUSV[sCount]
            sCount = sCount - 1
        end = 0
        eCount = x
        while eCount >= 0:
            end = end + nbUSV[eCount]
            eCount = eCount - 1
        for y in totalBandwith[start:end]: #for each char from last end point --> last + current num to get to next end index
            if y != ' ': #if y is not ' ' aka it is a number
                bwT = bwT + y #add that total bandwidth to the total of the length
                numBw = numBw + 1 #add 1 to the number of bandwidths in that recording
    bwM.append(bwT/numBw) #divide the total bandwidth value by the number of bandwidths in the recording, add that to the list



#amplitude
#average of amplitude for each recording
ampM = []
for x in range(len(nbUSV)): #for each of the separate recordings
    numAmp = 0 #number of amp values in that recording is 0 at start
    ampT = 0 #total amp to average the sequences is 0 at start
    if x == 0: #for the first recording
        for y in amplitude[0:nbUSV[x]]: #for first --> x char in length (first recording sequences end there)
            if y != ' ': #if y is not ' ' aka it is a number
                ampT = ampT + y #add that amp val to the total of the length
                numAmp = numAmp + 1 #add 1 to the number of amps in that recording
    else: #for the recordings after the first one
        start = 0
        sCount = x-1
        while sCount >= 0:
            start = start + nbUSV[sCount]
            sCount = sCount - 1
        end = 0
        eCount = x
        while eCount >= 0:
            end = end + nbUSV[eCount]
            eCount = eCount - 1
        for y in amplitude[start:end]: #for each char from last end point --> last + current num to get to next end index
            if y != ' ': #if y is not ' ' aka it is a number
                ampT = ampT + y #add that amp val to the total of the length
                numAmp = numAmp + 1 #add 1 to the number of amps in that recording
    ampM.append(ampT/numAmp) #divide the total amp value by the number of amps in the recording, add that to the list

#duration
#average duration for each recording
durM = []
for x in range(len(nbUSV)): #for each of the separate recordings
    numDur = 0 #number of duration values in that recording is 0 at start
    durT = 0 #total duration to average the sequences is 0 at start
    if x == 0: #for the first recording
        for y in duration[0:nbUSV[x]]: #for first --> x char in length (first recording sequences end there)
            if y != ' ': #if y is not ' ' aka it is a number
                durT = durT + y #add that duration val to the total of the length
                numDur = numDur + 1 #add 1 to the number of duration values in that recording
    else: #for the recordings after the first one
        start = 0
        sCount = x-1
        while sCount >= 0:
            start = start + nbUSV[sCount]
            sCount = sCount - 1
        end = 0
        eCount = x
        while eCount >= 0:
            end = end + nbUSV[eCount]
            eCount = eCount - 1
        for y in duration[start:end]: #for each char from last end point --> last + current num to get to next end index
            if y != ' ': #if y is not ' ' aka it is a number
                durT = durT + y #add that duration val to the total of the length
                numDur = numDur + 1 #add 1 to the number of duration values in that recording
    durM.append(durT/numDur) #divide the total duration value by the number of duration values in the recording, add that to the list

#fq variance
#avg fq variance for each recording
fqvM = []
for x in range(len(nbUSV)): #for each of the separate recordings
    numFqv = 0 #number of fqv values in that recording is 0 at start
    fqvT = 0 #total fqv to average the sequences is 0 at start
    if x == 0: #for the first recording
        for y in fqVariance[0:nbUSV[x]]: #for first --> x char in length (first recording sequences end there)
            if y != ' ': #if y is not ' ' aka it is a number
                fqvT = fqvT + y #add that fqv val to the total of the length
                numFqv = numFqv + 1 #add 1 to the number of fqv values in that recording
    else: #for the recordings after the first one
        start = 0
        sCount = x-1
        while sCount >= 0:
            start = start + nbUSV[sCount]
            sCount = sCount - 1
        end = 0
        eCount = x
        while eCount >= 0:
            end = end + nbUSV[eCount]
            eCount = eCount - 1
        for y in fqVariance[start:end]: #for each char from last end point --> last + current num to get to next end index
            if y != ' ': #if y is not ' ' aka it is a number
                fqvT = fqvT + y #add that fqv val to the total of the length
                numFqv = numFqv + 1 #add 1 to the number of fqv values in that recording
    fqvM.append(fqvT/numFqv) #divide the total fqv value by the number of fqv values in the recording, add that to the list


#spectral purity
#average spectral pruity for each recording
spM = []
for x in range(len(nbUSV)): #for each of the separate recordings
    numSp = 0 #number of spectral purity values in that recording is 0 at start
    spT = 0 #total spectral purity to average the sequences is 0 at start
    if x == 0: #for the first recording
        for y in spectralPurity[0:nbUSV[x]]: #for first --> x char in length (first recording sequences end there)
            if y != ' ': #if y is not ' ' aka it is a number
                spT = spT + y #add that spectral purity val to the total of the length
                numSp = numSp + 1 #add 1 to the number of spectral purity values in that recording
    else: #for the recordings after the first one
        start = 0
        sCount = x-1
        while sCount >= 0:
            start = start + nbUSV[sCount]
            sCount = sCount - 1
        end = 0
        eCount = x
        while eCount >= 0:
            end = end + nbUSV[eCount]
            eCount = eCount - 1
        for y in spectralPurity[start:end]: #for each char from last end point --> last + current num to get to next end index
            if y != ' ': #if y is not ' ' aka it is a number
                spT = spT + y #add that spectral purity val to the total of the length
                numSp = numSp + 1 #add 1 to the number of spectral purity values in that recording
    spM.append(spT/numSp) #divide the total spectral purity value by the number of fqv values in the recording, add that to the list


#ratio compldata/simple
#num seq oks/(total num sequences - seq ok)
csRatio = []
for ident in range(len(nbUSV)):
    if ident == 0: #for the first recording
        numOK = 0 #number of seq oks is 0
        numS = 0 #number of simple endings is 0
        for y in range(nbUSV[x]): #for first --> x char in nbUSV of the 1st recording
            if isinstance(length[y], int) and length[y] > 1: #if that char of seq 1 is STOP
                if numC[y] == 'SEQ OK': #if that char of numC is SEQ OK
                    numOK = numOK + 1 #add 1 to total num of seq ok
                else:
                    numS = numS + 1 #else add 1 to non seq ok numbers because it is the end of a seq but not seq ok
    else: #for the recordings after the first one
        numOK = 0 #number of seq oks is 0
        numS = 0 #number of simple endings is 0
        for y in range((np.sum(nbUSV[0:ident])), (np.sum(nbUSV[0:ident+1]))): #for each char from last end point --> last + current num to get to next end index
            if isinstance(length[y], int) and length[y] > 1: #if that char of seq1 is STOP
                if numC[y] == 'SEQ OK': #if that char of numC is SEQ OK
                    numOK = numOK + 1 #add 1 to total num of seq ok
                else:
                    numS = numS + 1 #else add 1 to non seq ok numbers because it is the end of a seq but not seq ok
    if numS > 0:
        csRatio.append(numOK/numS) #ratio = total seq ok/ total non seq ok in a recording?
    elif numS == 0:
        csRatio.append(0)

#s,d,u,m data below

sfqMin = []
sfqMean = []
sfqMax = []
sfqStart = []
sfqEnd = []
sBw = []
sAmp  =[]
sfqVar = []
sPurity = []
sDur = []


dfqMin = []
dfqMean = []
dfqMax = []
dfqStart = []
dfqEnd = []
dBw = []
dAmp  =[]
dfqVar = []
dPurity = []
dDur = []


ufqMin = []
ufqMean = []
ufqMax = []
ufqStart = []
ufqEnd = []
uBw = []
uAmp  =[]
ufqVar = []
uPurity = []
uDur = []


mfqMin = []
mfqMean = []
mfqMax = []
mfqStart = []
mfqEnd = []
mBw = []
mAmp  =[]
mfqVar = []
mPurity = []
mDur = []


for ident in range(len(nbUSV)):
    sfqMinT = 0
    sfqMeanT = 0
    sfqMaxT = 0
    sfqStartT = 0
    sfqStartT = 0
    sfqEndT = 0
    sBwT = 0
    sAmpT  =0
    sfqVarT = 0
    sPurityT = 0
    sDurT = 0
    dfqMinT = 0
    dfqMeanT = 0
    dfqMaxT = 0
    dfqStartT = 0
    dfqStartT = 0
    dfqEndT = 0
    dBwT = 0
    dAmpT  =0
    dfqVarT = 0
    dPurityT = 0
    dDurT = 0
    ufqMinT = 0
    ufqMeanT = 0
    ufqMaxT = 0
    ufqStartT = 0
    ufqStartT = 0
    ufqEndT = 0
    uBwT = 0
    uAmpT  =0
    ufqVarT = 0
    uPurityT = 0
    uDurT = 0
    mfqMinT = 0
    mfqMeanT = 0
    mfqMaxT = 0
    mfqStartT = 0
    mfqStartT = 0
    mfqEndT = 0
    mBwT = 0
    mAmpT  =0
    mfqVarT = 0
    mPurityT = 0
    mDurT = 0
    sCounter = 0
    dCounter = 0
    uCounter = 0 
    mCounter = 0
    if ident == 0: #for the first recording
        for y in range(nbUSV[x]): #for first --> x char in nbUSV of the 1st recording
            if syl[y] == 's': #if that char of seq 1 is STOP
                sfqMinT = sfqMinT + fqMin[y]
                sfqMeanT = sfqMeanT + fqMean[y]
                sfqMaxT = sfqMaxT + fqMax[y]
                sfqStartT = sfqStartT + fqStart[y]
                sfqEndT = sfqEndT + fqEnd[y]
                sBwT = sBwT + totalBandwith[y]
                sAmpT = sAmpT + amplitude[y]
                sfqVarT = sfqVarT + fqVariance[y]
                sPurityT = sPurityT + spectralPurity[y]
                sDurT = sDurT + duration[y]
                sCounter = sCounter + 1
            elif syl[y] == 'd':
                dfqMinT = dfqMinT + fqMin[y]
                dfqMeanT = dfqMeanT + fqMean[y]
                dfqMaxT = dfqMaxT + fqMax[y]
                dfqStartT = dfqStartT + fqStart[y]
                dfqEndT = dfqEndT + fqEnd[y]
                dBwT = dBwT + totalBandwith[y]
                dAmpT = dAmpT + amplitude[y]
                dfqVarT = dfqVarT + fqVariance[y]
                dPurityT = dPurityT + spectralPurity[y]
                dDurT = dDurT + duration[y]
                dCounter = dCounter + 1
            elif syl[y] == 'u':
                ufqMinT = ufqMinT + fqMin[y]
                ufqMeanT = ufqMeanT + fqMean[y]
                ufqMaxT = ufqMaxT + fqMax[y]
                ufqStartT = ufqStartT + fqStart[y]
                ufqEndT = ufqEndT + fqEnd[y]
                uBwT = uBwT + totalBandwith[y]
                uAmpT = uAmpT + amplitude[y]
                ufqVarT = ufqVarT + fqVariance[y]
                uPurityT = uPurityT + spectralPurity[y]
                uDurT = uDurT + duration[y]
                uCounter = uCounter + 1
            elif syl[y] == 'M':
                mfqMinT = mfqMinT + fqMin[y]
                mfqMeanT = mfqMeanT + fqMean[y]
                mfqMaxT = mfqMaxT + fqMax[y]
                mfqStartT = mfqStartT + fqStart[y]
                mfqEndT = mfqEndT + fqEnd[y]
                mBwT = mBwT + totalBandwith[y]
                mAmpT = mAmpT + amplitude[y]
                mfqVarT = mfqVarT + fqVariance[y]
                mPurityT = mPurityT + spectralPurity[y]
                mDurT = mDurT + duration[y]
                mCounter = mCounter + 1
    else:
        for y in range((np.sum(nbUSV[0:ident])), (np.sum(nbUSV[0:ident+1]))):
            if syl[y] == 's': #if that char of seq 1 is STOP
                sfqMinT = sfqMinT + fqMin[y]
                sfqMeanT = sfqMeanT + fqMean[y]
                sfqMaxT = sfqMaxT + fqMax[y]
                sfqStartT = sfqStartT + fqStart[y]
                sfqEndT = sfqEndT + fqEnd[y]
                sBwT = sBwT + totalBandwith[y]
                sAmpT = sAmpT + amplitude[y]
                sfqVarT = sfqVarT + fqVariance[y]
                sPurityT = sPurityT + spectralPurity[y]
                sDurT = sDurT + duration[y]
                sCounter = sCounter + 1
            elif syl[y] == 'd':
                dfqMinT = dfqMinT + fqMin[y]
                dfqMeanT = dfqMeanT + fqMean[y]
                dfqMaxT = dfqMaxT + fqMax[y]
                dfqStartT = dfqStartT + fqStart[y]
                dfqEndT = dfqEndT + fqEnd[y]
                dBwT = dBwT + totalBandwith[y]
                dAmpT = dAmpT + amplitude[y]
                dfqVarT = dfqVarT + fqVariance[y]
                dPurityT = dPurityT + spectralPurity[y]
                dDurT = dDurT + duration[y]
                dCounter = dCounter + 1
            elif syl[y] == 'u':
                ufqMinT = ufqMinT + fqMin[y]
                ufqMeanT = ufqMeanT + fqMean[y]
                ufqMaxT = ufqMaxT + fqMax[y]
                ufqStartT = ufqStartT + fqStart[y]
                ufqEndT = ufqEndT + fqEnd[y]
                uBwT = uBwT + totalBandwith[y]
                uAmpT = uAmpT + amplitude[y]
                ufqVarT = ufqVarT + fqVariance[y]
                uPurityT = uPurityT + spectralPurity[y]
                uDurT = uDurT + duration[y]
                uCounter = uCounter + 1
            elif syl[y] == 'M':
                mfqMinT = mfqMinT + fqMin[y]
                mfqMeanT = mfqMeanT + fqMean[y]
                mfqMaxT = mfqMaxT + fqMax[y]
                mfqStartT = mfqStartT + fqStart[y]
                mfqEndT = mfqEndT + fqEnd[y]
                mBwT = mBwT + totalBandwith[y]
                mAmpT = mAmpT + amplitude[y]
                mfqVarT = mfqVarT + fqVariance[y]
                mPurityT = mPurityT + spectralPurity[y]
                mDurT = mDurT + duration[y]
                mCounter = mCounter + 1
    if sCounter > 0:
        sfqMin.append(sfqMinT/sCounter)
        sfqMean.append(sfqMeanT/sCounter)
        sfqMax.append(sfqMaxT/sCounter)
        sfqStart.append(sfqStartT/sCounter) 
        sfqEnd.append(sfqEndT/sCounter) 
        sBw.append(sBwT/sCounter)
        sAmp.append(sAmpT/sCounter)
        sfqVar.append(sfqVarT/sCounter)
        sPurity.append(sPurityT/sCounter)
        sDur.append(sDurT/sCounter)
    elif sCounter == 0:
        sfqMin.append(np.nan)
        sfqMean.append(np.nan)
        sfqMax.append(np.nan)
        sfqStart.append(np.nan) 
        sfqEnd.append(np.nan) 
        sBw.append(np.nan)
        sAmp.append(np.nan)
        sfqVar.append(np.nan)
        sPurity.append(np.nan)
        sDur.append(np.nan)
    if dCounter > 0:
        dfqMin.append(dfqMinT/dCounter)
        dfqMean.append(dfqMeanT/dCounter)
        dfqMax.append(dfqMaxT/dCounter)
        dfqStart.append(dfqStartT/dCounter) 
        dfqEnd.append(dfqEndT/dCounter) 
        dBw.append(dBwT/dCounter)
        dAmp.append(dAmpT/dCounter)
        dfqVar.append(dfqVarT/dCounter)
        dPurity.append(dPurityT/dCounter)
        dDur.append(dDurT/dCounter)
    elif dCounter == 0:
        dfqMin.append(np.nan)
        dfqMean.append(np.nan)
        dfqMax.append(np.nan)
        dfqStart.append(np.nan) 
        dfqEnd.append(np.nan) 
        dBw.append(np.nan)
        dAmp.append(np.nan)
        dfqVar.append(np.nan)
        dPurity.append(np.nan)
        dDur.append(np.nan)
    
    if uCounter > 0:
        ufqMin.append(ufqMinT/uCounter)
        ufqMean.append(ufqMeanT/uCounter)
        ufqMax.append(ufqMaxT/uCounter)
        ufqStart.append(ufqStartT/uCounter) 
        ufqEnd.append(ufqEndT/uCounter) 
        uBw.append(uBwT/uCounter)
        uAmp.append(uAmpT/uCounter)
        ufqVar.append(ufqVarT/uCounter)
        uPurity.append(uPurityT/uCounter)
        uDur.append(uDurT/uCounter)
    elif uCounter == 0:
        ufqMin.append(np.nan)
        ufqMean.append(np.nan)
        ufqMax.append(np.nan)
        ufqStart.append(np.nan) 
        ufqEnd.append(np.nan) 
        uBw.append(np.nan)
        uAmp.append(np.nan)
        ufqVar.append(np.nan)
        uPurity.append(np.nan)
        uDur.append(np.nan)
    
    if mCounter > 0:
        mfqMin.append(mfqMinT/mCounter)
        mfqMean.append(mfqMeanT/mCounter)
        mfqMax.append(mfqMaxT/mCounter)
        mfqStart.append(mfqStartT/mCounter) 
        mfqEnd.append(mfqEndT/mCounter) 
        mBw.append(mBwT/mCounter)
        mAmp.append(mAmpT/mCounter)
        mfqVar.append(mfqVarT/mCounter)
        mPurity.append(mPurityT/mCounter)
        mDur.append(mDurT/mCounter)
    elif mCounter == 0:
        mfqMin.append(np.nan)
        mfqMean.append(np.nan)
        mfqMax.append(np.nan)
        mfqStart.append(np.nan) 
        mfqEnd.append(np.nan) 
        mBw.append(np.nan)
        mAmp.append(np.nan)
        mfqVar.append(np.nan)
        mPurity.append(np.nan)
        mDur.append(np.nan)

sAvgs = [sfqMin, sfqMean, sfqMax, sfqStart, sfqEnd, sBw, sAmp, sfqVar, sPurity, sDur]
dAvgs = [dfqMin, dfqMean, dfqMax, dfqStart, dfqEnd, dBw, dAmp, dfqVar, dPurity, dDur]
uAvgs = [ufqMin, ufqMean, ufqMax, ufqStart, ufqEnd, uBw, uAmp, ufqVar, uPurity, uDur]
mAvgs = [mfqMin, mfqMean, mfqMax, mfqStart, mfqEnd, mBw, mAmp, mfqVar, mPurity, mDur]


leedleLEE = []
for x in range(len(batchNew)):
	newl = [batchNew[x], nbUSV[x], callRate[x], \
         listS[x], listD[x], listU[x], listM[x], \
         listSt[x], listDt[x], listUt[x], listMt[x], \
         seqLen[x], fqM[x], bwM[x], ampM[x], durM[x], fqvM[x], spM[x], csRatio[x], \
         sfqMin[x], sfqMean[x], sfqMax[x], sfqStart[x], sfqEnd[x], sBw[x], sAmp[x], sfqVar[x], sPurity[x], sDur[x], \
         dfqMin[x], dfqMean[x], dfqMax[x], dfqStart[x], dfqEnd[x], dBw[x], dAmp[x], dfqVar[x], dPurity[x], dDur[x], \
         ufqMin[x], ufqMean[x], ufqMax[x], ufqStart[x], ufqEnd[x], uBw[x], uAmp[x], ufqVar[x], uPurity[x], uDur[x], \
         mfqMin[x], mfqMean[x], mfqMax[x], mfqStart[x], mfqEnd[x], mBw[x], mAmp[x], mfqVar[x], mPurity[x], mDur[x]]
	leedleLEE.append(newl)


featureNames = ['ID', 'Number of USVs', 'Call Rate', \
                'percent starting s', 'percent starting d', 'percent starting u', 'percent starting m', \
                'percent composition s', 'percent composition u', 'percent composition d', 'percent composition m', \
                'sequence Length', 'Fq Mean', 'Bandwidth Mean', 'Amplitude Mean', 'Duraiton Mean',\
                    'Fq Variance Mean', 'Spectral Purity Mean', 'Complex to Simple Ratio', \
                'sfqMin', 'sfqMean', 'sfqMax', 'sfqStart', 'sfqEnd', 'sBw', 'sAmp', 'sfqVar', 'sPurity', 'sDur', \
                'dfqMin', 'dfqMean', 'dfqMax', 'dfqStart', 'dfqEnd', 'dBw', 'dAmp', 'dfqVar', 'dPurity', 'dDur', \
                'ufqMin', 'ufqMean', 'ufqMax', 'ufqStart', 'ufqEnd', 'uBw', 'uAmp', 'ufqVar', 'uPurity', 'uDur', \
                'mfqMin', 'mfqMean', 'mfqMax', 'mfqStart', 'mfqEnd', 'mBw', 'mAmp', 'mfqVar', 'mPurity', 'mDur']

with open(newnewDataFolder + '/Features.csv', 'w', newline = '') as results:
	writer = csv.writer(results)
	writer.writerow(featureNames)
	for hoopla in leedleLEE:
		writer.writerow(hoopla)
        
########################################################################################
'''OVERALL PORBABILITIES'''

overallProbSs = []
overallProbSd = []
overallProbSu = []
overallProbSm = []
overallProbSx = []
overallProbDs = []
overallProbDd = []
overallProbDu = []
overallProbDm = []
overallProbDx = []
overallProbUs = []
overallProbUd = []
overallProbUu = []
overallProbUm = []
overallProbUx = []
overallProbMs = []
overallProbMd = []
overallProbMu = []
overallProbMm = []
overallProbMx = []
overallProbXs = []
overallProbXd = []
overallProbXu = []
overallProbXm = []


for ident in range(len(nbUSV)):
    ss = 0
    sd = 0
    su = 0
    sm = 0
    sx = 0
    ds = 0
    dd = 0
    du = 0
    dm = 0
    dx = 0
    us = 0
    ud = 0
    uu = 0
    um = 0
    ux = 0 
    ms = 0
    md = 0
    mu = 0
    mm = 0
    mx = 0
    xs = 0
    xd = 0
    xu = 0
    xm = 0
    tot = 0
    if ident == 0:
        #takes every element from index 0 to the last index of the first recording
        for y in transSeq[0:nbUSV[0]]: 
            tot = tot + 1
            if y == 'ss':
                ss = ss + 1
            elif y == 'sd':
                sd = sd + 1
            elif y == 'su':
                su = su + 1
            elif y == 'sM':
                sm = sm + 1
            elif y == 'sX':
                sx = sx + 1
            elif y == 'ds':
                ds = ds + 1
            elif y == 'dd':
                dd = dd + 1
            elif y == 'du':
                du = du + 1
            elif y == 'dM':
                dm = dm + 1
            elif y == 'dX':
                dx = dx + 1
            elif y == 'us':
                us = us + 1
            elif y == 'ud':
                ud = ud + 1
            elif y == 'uu':
                uu = uu + 1
            elif y == 'uM':
                um = um + 1
            elif y == 'uX':
                ux = ux + 1
            elif y == 'Ms':
                ms = ms + 1
            elif y == 'Md':
                md = md + 1
            elif y == 'Mu':
                mu = mu + 1
            elif y == 'MM':
                mm = mm + 1
            elif y == 'MX':
                mx = mx + 1
        for trans in transFromX[0:nbUSV[0]]:
            if trans == 'Xs':
                xs = xs + 1
                tot = tot + 1
            elif trans == 'Xd':
                xd = xd + 1
                tot = tot + 1
            elif trans == 'Xu':
                xu = xu + 1
                tot = tot + 1
            elif trans == 'XM':
                xm = xm + 1
                tot = tot + 1
        print('%s: %s USVs'%(batchNew[ident],nbUSV[ident]))
    else:
        #-1 here to make sure for loop does not include last transition/end of recording
        for y in range((np.sum(nbUSV[0:ident])), ((np.sum(nbUSV[0:ident+1])))): 
            if transSeq[y] == 'ss':
                ss = ss + 1
            elif transSeq[y] == 'sd':
                sd = sd + 1
            elif transSeq[y] == 'su':
                su = su + 1
            elif transSeq[y] == 'sM':
                sm = sm + 1
            elif transSeq[y] == 'sX':
                sx = sx + 1
            elif transSeq[y] == 'ds':
                ds = ds + 1
            elif transSeq[y] == 'dd':
                dd = dd + 1
            elif transSeq[y] == 'du':
                du = du + 1
            elif transSeq[y] == 'dM':
                dm = dm + 1
            elif transSeq[y] == 'dX':
                dx = dx + 1
            elif transSeq[y] == 'us':
                us = us + 1
            elif transSeq[y] == 'ud':
                ud = ud + 1
            elif transSeq[y] == 'uu':
                uu = uu + 1
            elif transSeq[y] == 'uM':
                um = um + 1
            elif transSeq[y] == 'uX':
                ux = ux + 1
            elif transSeq[y] == 'Ms':
                ms = ms + 1
            elif transSeq[y] == 'Md':
                md = md + 1
            elif transSeq[y] == 'Mu':
                mu = mu + 1
            elif transSeq[y] == 'MM':
                mm = mm + 1
            elif transSeq[y] == 'MX':
                mx = mx + 1
            elif transSeq[y] == 'Xs':
                xs = xs + 1
            elif transSeq[y] == 'Xd':
                xd = xd + 1
            elif transSeq[y] == 'Xu':
                xu = xu + 1
            elif transSeq[y] == 'XM':
                xm = xm + 1
            tot = tot + 1

        for trans in range((np.sum(nbUSV[0:ident])), ((np.sum(nbUSV[0:ident+1])))):
            if transFromX[trans] == 'Xs':
                xs = xs + 1
                tot = tot + 1
            elif transFromX[trans] == 'Xd':
                xd = xd + 1
                tot = tot + 1
            elif transFromX[trans] == 'Xu':
                xu = xu + 1
                tot = tot + 1
            elif transFromX[trans] == 'XM':
                xm = xm + 1
                tot = tot + 1
        print('%s: %s USVs'%(batchNew[ident],nbUSV[ident]))
        
    if ss > 0:
        overallProbSs.append(ss/tot)
    elif ss == 0:
        overallProbSs.append(0)
    if sd > 0:
        overallProbSd.append(sd/tot)
    elif sd == 0:
        overallProbSd.append(0)
    if su > 0:
        overallProbSu.append(su/tot)
    elif su == 0:
        overallProbSu.append(0)
    if sm > 0:
        overallProbSm.append(sm/tot)
    elif sm == 0:
        overallProbSm.append(0)
    if sx > 0:
        overallProbSx.append(sx/tot)
    elif sx == 0:
        overallProbSx.append(0)
    if ds > 0:
        overallProbDs.append(ds/tot)
    elif ds == 0:
        overallProbDs.append(0)
    if dd > 0:
        overallProbDd.append(dd/tot)
    elif dd == 0:
        overallProbDd.append(0)
    if du > 0:
        overallProbDu.append(du/tot)
    elif du == 0:
        overallProbDu.append(0)
    if dm > 0:
        overallProbDm.append(dm/tot)
    elif dm == 0:
        overallProbDm.append(0)
    if dx > 0:
        overallProbDx.append(dx/tot)
    elif dx == 0:
        overallProbDx.append(0)
    if us > 0:
        overallProbUs.append(us/tot)
    elif us == 0:
        overallProbUs.append(0)
    if ud > 0:
        overallProbUd.append(ud/tot)
    elif ud == 0:
        overallProbUd.append(0)
    if uu > 0:
        overallProbUu.append(uu/tot)
    elif uu == 0:
        overallProbUu.append(0)
    if um > 0:
        overallProbUm.append(um/tot)
    elif um == 0:
        overallProbUm.append(0)
    if ux > 0:
        overallProbUx.append(ux/tot)
    elif ux == 0:
        overallProbUx.append(0)
    if ms > 0:
        overallProbMs.append(ms/tot)
    elif ms == 0:
        overallProbMs.append(0)
    if md > 0:
        overallProbMd.append(md/tot)
    elif md == 0:
        overallProbMd.append(0)
    if mu > 0:
        overallProbMu.append(mu/tot)
    elif mu == 0:
        overallProbMu.append(0)
    if mm > 0:
        overallProbMm.append(mm/tot)
    elif mm == 0:
        overallProbMm.append(0)
    if mx > 0:
        overallProbMx.append(mx/tot)
    elif mx == 0:
        overallProbMx.append(0)
    if xs > 0:
        overallProbXs.append(xs/tot)
    elif xs == 0:
        overallProbXs.append(0)
    if xd > 0:
        overallProbXd.append(xd/tot)
    elif xd == 0:
        overallProbXd.append(0)
    if xu > 0:
        overallProbXu.append(xu/tot)
    elif xu == 0:
        overallProbXu.append(0)
    if xm > 0:
        overallProbXm.append(xm/tot)
    elif xm == 0:
        overallProbXm.append(0)


isMayonnaiseAnInstrument = []
for x in range(len(batchNew)):
	newl = [batchNew[x],
         overallProbSs[x], overallProbSd[x], overallProbSu[x], overallProbSm[x], overallProbSx[x], \
         overallProbDs[x], overallProbDd[x], overallProbDu[x], overallProbDm[x], overallProbDx[x], \
         overallProbUs[x], overallProbUd[x], overallProbUu[x], overallProbUm[x], overallProbUx[x], \
         overallProbMs[x], overallProbMd[x], overallProbMu[x], overallProbMm[x], overallProbMx[x], \
         overallProbXs[x], overallProbXd[x], overallProbXu[x], overallProbXm[x]]
	isMayonnaiseAnInstrument.append(newl)


overallNames = ['ID', \
         'overallProb Ss', 'overallProb Sd', 'overallProb Su', 'overallProb Sm', 'overallProb Sx', \
         'overallProb Ds', 'overallProb Dd', 'overallProb Du', 'overallProb Dm', 'overallProb Dx', \
         'overallProb Us', 'overallProb Ud', 'overallProb Uu', 'overallProb Um', 'overallProb Ux', \
         'overallProb Ms', 'overallProb Md', 'overallProb Mu', 'overallProb Mm', 'overallProb Mx', \
         'overallProb Xs', 'overallProb Xd', 'overallProb Xu', 'overallProb Xm']

with open(newnewDataFolder + '/TransitionProbs_overall.csv', 'w', newline = '') as results:
	writer = csv.writer(results)
	writer.writerow(overallNames)
	for hoopla in isMayonnaiseAnInstrument:
		writer.writerow(hoopla)


########################################################################################
'''CONDITIONAL PROBABILITIES'''

condProbSs = []
condProbSd = []
condProbSu = []
condProbSm = []
condProbSx = []
condProbDs = []
condProbDd = []
condProbDu = []
condProbDm = []
condProbDx = []
condProbUs = []
condProbUd = []
condProbUu = []
condProbUm = []
condProbUx = []
condProbMs = []
condProbMd = []
condProbMu = []
condProbMm = []
condProbMx = []
condProbXs = []
condProbXd = []
condProbXu = []
condProbXm = []


for ident in range(len(nbUSV)):
    ssc = 0
    sdc=0
    suc = 0
    smc = 0
    sxc = 0
    dsc = 0
    ddc = 0
    duc = 0
    dmc = 0
    dxc = 0
    usc = 0
    udc = 0
    uuc = 0
    umc = 0
    uxc = 0 
    msc = 0
    mdc = 0
    muc = 0
    mmc = 0
    mxc = 0
    xsc = 0
    xdc = 0
    xuc = 0
    xmc = 0
    sTot = 0
    dTot = 0
    uTot = 0
    mTot = 0
    xTot = 0
    if ident == 0:
        for y in transSeq[0:nbUSV[0]]:
            if y == 'ss':
                ssc = ssc + 1
                sTot =sTot + 1
            elif y == 'sd':
                sdc = sdc + 1
                sTot = sTot + 1
            elif y == 'su':
                suc = suc + 1
                sTot = sTot + 1
            elif y == 'sM':
                smc = smc + 1
                sTot = sTot + 1
            elif y == 'sX':
                sxc = sxc + 1
                sTot = sTot + 1
            elif y == 'ds':
                dsc = dsc + 1
                dTot = dTot + 1
            elif y == 'dd':
                ddc = ddc + 1
                dTot = dTot + 1
            elif y == 'du':
                duc = duc + 1
                dTot = dTot + 1
            elif y == 'dM':
                dmc = dmc + 1
                dTot = dTot + 1
            elif y == 'dX':
                dxc = dxc + 1
                dTot = dTot + 1
            elif y == 'us':
                usc = usc + 1
                uTot = uTot + 1
            elif y == 'ud':
                udc = udc + 1
                uTot = uTot + 1
            elif y == 'uu':
                uuc = uuc + 1
                uTot = uTot + 1
            elif y == 'uM':
                umc = umc + 1
                uTot = uTot + 1
            elif y == 'uX':
                uxc = uxc + 1
                uTot = uTot + 1
            elif y == 'Ms':
                msc = msc + 1
                mTot = mTot + 1
            elif y == 'Md':
                mdc = mdc + 1
                mTot = mTot + 1
            elif y == 'Mu':
                muc = muc + 1
                mTot = mTot + 1
            elif y == 'MM':
                mmc = mmc + 1
                mTot = mTot + 1
            elif y == 'MX':
                mxc = mxc + 1
                mTot = mTot + 1
        for y in transFromX[0:nbUSV[0]]:
            if y == 'Xs':
                xsc = xsc + 1
                xTot = xTot + 1
            elif y == 'Xd':
                xdc = xdc + 1
                xTot = xTot + 1
            elif y == 'Xu':
                xuc = xuc + 1
                xTot = xTot + 1
            elif y == 'XM':
                xmc = xmc + 1
                xTot = xTot + 1
    else:
        for y in range((np.sum(nbUSV[0:ident])), ((np.sum(nbUSV[0:ident+1])))):
            if transSeq[y] == 'ss':
                ssc = ssc + 1
                sTot =sTot + 1
            elif transSeq[y] == 'sd':
                sdc = sdc + 1
                sTot = sTot + 1
            elif transSeq[y] == 'su':
                suc = suc + 1
                sTot = sTot + 1
            elif transSeq[y] == 'sM':
                smc = smc + 1
                sTot = sTot + 1
            elif transSeq[y] == 'sX':
                sxc = sxc + 1
                sTot = sTot + 1
            elif transSeq[y] == 'ds':
                dsc = dsc + 1
                dTot = dTot + 1
            elif transSeq[y] == 'dd':
                ddc = ddc + 1
                dTot = dTot + 1
            elif transSeq[y] == 'du':
                duc = duc + 1
                dTot = dTot + 1
            elif transSeq[y] == 'dM':
                dmc = dmc + 1
                dTot = dTot + 1
            elif transSeq[y] == 'dX':
                dxc = dxc + 1
                dTot = dTot + 1
            elif transSeq[y] == 'us':
                usc = usc + 1
                uTot = uTot + 1
            elif transSeq[y] == 'ud':
                udc = udc + 1
                uTot = uTot + 1
            elif transSeq[y] == 'uu':
                uuc = uuc + 1
                uTot = uTot + 1
            elif transSeq[y] == 'uM':
                umc = umc + 1
                uTot = uTot + 1
            elif transSeq[y] == 'uX':
                uxc = uxc + 1
                uTot = uTot + 1
            elif transSeq[y] == 'Ms':
                msc = msc + 1
                mTot = mTot + 1
            elif transSeq[y] == 'Md':
                mdc = mdc + 1
                mTot = mTot + 1
            elif transSeq[y] == 'Mu':
                muc = muc + 1
                mTot = mTot + 1
            elif transSeq[y] == 'MM':
                mmc = mmc + 1
                mTot = mTot + 1
            elif transSeq[y] == 'MX':
                mxc = mxc + 1
                mTot = mTot + 1
        for trans in range((np.sum(nbUSV[0:ident])), ((np.sum(nbUSV[0:ident+1])))):
            if transFromX[trans] == 'Xs':
                xsc = xsc + 1
                xTot = xTot + 1
            elif transFromX[trans] == 'Xd':
                xdc = xdc + 1
                xTot = xTot + 1
            elif transFromX[trans] == 'Xu':
                xuc = xuc + 1
                xTot = xTot + 1
            elif transFromX[trans] == 'XM':
                xmc = xmc + 1
                xTot = xTot + 1
    if sTot > 0:
        condProbSs.append(ssc/sTot)
        condProbSd.append(sdc/sTot)
        condProbSu.append(suc/sTot)
        condProbSm.append(smc/sTot)
        condProbSx.append(sxc/sTot)
    elif sTot == 0:
        condProbSs.append(0)
        condProbSd.append(0)
        condProbSu.append(0)
        condProbSm.append(0)
        condProbSx.append(0)
    if dTot > 0:
        condProbDs.append(dsc/dTot)
        condProbDd.append(ddc/dTot)
        condProbDu.append(duc/dTot)
        condProbDm.append(dmc/dTot)
        condProbDx.append(dxc/dTot)
    elif dTot == 0:
        condProbDs.append(0)
        condProbDd.append(0)
        condProbDu.append(0)
        condProbDm.append(0)
        condProbDx.append(0)
    if uTot > 0:
        condProbUs.append(usc/uTot)
        condProbUd.append(udc/uTot)
        condProbUu.append(uuc/uTot)
        condProbUm.append(umc/uTot)
        condProbUx.append(uxc/uTot)
    elif uTot == 0:
        condProbUs.append(0)
        condProbUd.append(0)
        condProbUu.append(0)
        condProbUm.append(0)
        condProbUx.append(0)
    if mTot > 0:
        condProbMs.append(msc/mTot)
        condProbMd.append(mdc/mTot)
        condProbMu.append(muc/mTot)
        condProbMm.append(mmc/mTot)
        condProbMx.append(mxc/mTot)
    elif mTot == 0:
        condProbMs.append(0)
        condProbMd.append(0)
        condProbMu.append(0)
        condProbMm.append(0)
        condProbMx.append(0)
    if xTot > 0:
        condProbXs.append(xsc/xTot)
        condProbXd.append(xdc/xTot)
        condProbXu.append(xuc/xTot)
        condProbXm.append(xmc/xTot)
    elif xTot == 0: 
        condProbXs.append(0)
        condProbXd.append(0)
        condProbXu.append(0)
        condProbXm.append(0)


def interval_round(value):
    return round(value*40)/40

intervals = np.arange(0.0, 1.025, 0.025)
isi_filt = [x for x in isi if x == x] #removes NaN
isi_c = isi_filt

for v in range(0, len(isi_filt)):
    if isi_filt[v] > 1:
        isi_c[v] = 1
isi_r = [interval_round(val) for val in isi_c]

histISI = []

for inter in intervals:
    indexer = 0
    countz = 0
    for y in isi_r:
        if y == inter:
            countz = countz + 1
    histISI.append(countz/len(isi_r))

fig, ax = plt.subplots()
xnew = np.linspace(intervals.min(), intervals.max(), 300) 
spl = make_interp_spline(intervals, histISI, k=3)  # type: BSpline
ISIsmooth = spl(xnew)
ax.plot(xnew, ISIsmooth, color = 'k', linewidth = 0.5)
ax.set_xlabel('ISI Duration (s)')
ax.set_ylabel('Relative\nFrequency')
ax.set_title('ISI Density')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(newnewDataFolder + '/ISIDensity.png', dpi = 300)

smittywerbenjagermanjensen = []
for x in range(len(batchNew)):
    newl = [batchNew[x],
         condProbSs[x], condProbSd[x], condProbSu[x], condProbSm[x], condProbSx[x], \
         condProbDs[x], condProbDd[x], condProbDu[x], condProbDm[x], condProbDx[x], \
         condProbUs[x], condProbUd[x], condProbUu[x], condProbUm[x], condProbUx[x], \
         condProbMs[x], condProbMd[x], condProbMu[x], condProbMm[x], condProbMx[x], \
         condProbXs[x], condProbXd[x], condProbXu[x], condProbXm[x]]
    smittywerbenjagermanjensen.append(newl)


condNames = ['ID', \
         'condProb Ss', 'condProb Sd', 'condProb Su', 'condProb Sm', 'condProb Sx', \
         'condProb Ds', 'condProb Dd', 'condProb Du', 'condProb Dm', 'condProb Dx', \
         'condProb Us', 'condProb Ud', 'condProb Uu', 'condProb Um', 'condProb Ux', \
         'condProb Ms', 'condProb Md', 'condProb Mu', 'condProb Mm', 'condProb Mx', \
         'condProb Xs', 'condProb Xd', 'condProb Xu', 'condProb Xm']

with open(newnewDataFolder + '/TransitionProbs_conditional.csv', 'w', newline = '') as results:
    writer = csv.writer(results)
    writer.writerow(condNames)
    for hoopla in smittywerbenjagermanjensen:
        writer.writerow(hoopla)
print('Done!')

           