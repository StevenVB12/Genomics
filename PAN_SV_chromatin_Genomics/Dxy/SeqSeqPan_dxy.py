#!/usr/bin/env python

from __future__ import division
#from Bio.Seq import Seq
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--infile", help="extract regions from seq-seq-pan", required = True)
parser.add_argument("-g", "--genomeList", help="comma separated list of genomes in order of seq-seq-pan", required = True)
parser.add_argument("-w", "--windowSize", help="number", required = True)
parser.add_argument("-o", "--outfile", help="txt", required = True)
args = parser.parse_args()

infile = args.infile
genomes = args.genomeList.split(",")
window = args.windowSize
outfile = args.outfile

count = 0

genome1 = genomes[0]
genome2 = genomes[1]

identified = 0
sequence1 = str()

count +=1
fileIn = open(infile, "r")
fileLine = fileIn.readline()
fileLine = fileLine.strip()

print('searching genome ' + str(count) + ' ' + str(genome1))

while fileLine:

    if '> ' + str(genome1) + ':' in fileLine:
           
        identified += 1
        
        fileLine = fileIn.readline()
        fileLine = fileLine.strip()
        seq = fileLine
            
        sequence1 = sequence1 + seq

    if '=' not in fileLine and '>' not in fileLine:
        lenSeq = len(fileLine)

    if '=' in fileLine and identified == 0:
        sequence1 = sequence1 + '-' * lenSeq

    if '=' in fileLine and identified == 1:
        identified = 0
    
    fileLine = fileIn.readline()
    fileLine = fileLine.strip()

print('total length: ' + str(len(sequence1)))        

identified = 0
sequence2 = str()

count +=1
fileIn = open(infile, "r")
fileLine = fileIn.readline()
fileLine = fileLine.strip()

print('searching genome ' + str(count) + ' ' + str(genome2))

###
while fileLine:

    if '> ' + str(genome2) + ':' in fileLine:
           
        identified += 1
        
        fileLine = fileIn.readline()
        fileLine = fileLine.strip()
        seq = fileLine
            
        sequence2 = sequence2 + seq

    if '=' not in fileLine and '>' not in fileLine:
        lenSeq = len(fileLine)

    if '=' in fileLine and identified == 0:
        sequence2 = sequence2 + '-' * lenSeq

    if '=' in fileLine and identified == 1:
        identified = 0
    
    fileLine = fileIn.readline()
    fileLine = fileLine.strip()

print('total length: ' + str(len(sequence2)))        

f = open(outfile + '.txt',"a+")
f.write("\t".join(("start", "end", "good_sites", "nonSNP", "SNP", "\n")))

start = 0
end = int(window)
windcount = 0

      
while int(end) <= len(sequence1):
    
    windcount += 1

    if windcount % 100 == 0:
         print(str(windcount) + ' windows done')
    
    seq1 = sequence1[start:end]
    seq2 = sequence2[start:end]
    
    good = 0
    countSame = 0
    countDiff = 0
    for i in range(0,len(seq1)):

        base1 = seq1[i]
        base2 = seq2[i]

        if base1 != '-' and base1 != 'N' and base2 != '-' and base2 != 'N':
            good += 1
            
            if base1 == base2:
                countSame +=1
            if base1 != base2:
                countDiff +=1
                
    f.write("\t".join((str(start),
                       str(end),
                       str(good),
                       str(countSame),
                       str(countDiff), "\n")))   
            
    start += int(window)
    end += int(window)

print('all done')    
f.close()

