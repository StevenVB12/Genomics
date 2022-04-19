#!/usr/bin/env python

from __future__ import division
from Bio.Seq import Seq
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--infile", help="extract regions from seq-seq-pan", required = True)
parser.add_argument("-g", "--genomeList", help="comma seperated list of genomes in order of seq-seq-pan", required = True)
args = parser.parse_args()

infile = args.infile
genomes = args.genomeList.split(",")


count = 0
for genome in genomes:

    identified = 0
    sequence = str()

    
    f= open(genome + '_missing_intervals.txt',"a+")
    f.write("\t".join(("start","end","\n")))
    count +=1
    fileIn = open(infile, "r")
    fileLine = fileIn.readline()
    fileLine = fileLine.strip()
    
    print('searching genome ' + str(genome) + ' ' + str(genome))
    
    while fileLine:
    
        if '> ' + str(genome) + ':' in fileLine:
            
            if ' + ' in fileLine:
                inverse = False
            if ' - ' in fileLine:
                inverse = True
            
            identified += 1
            
            fileLine = fileIn.readline()
            fileLine = fileLine.strip()
            
            seq = Seq(fileLine)
            
            if inverse == True:
                seq = seq.reverse_complement()
                
            sequence = sequence + seq

        if '=' not in fileLine and '>' not in fileLine:
            lenSeq = len(fileLine)

        if '=' in fileLine and identified == 0:
            sequence = sequence + '-' * lenSeq

        if '=' in fileLine and identified == 1:
            identified = 0
        
        fileLine = fileIn.readline()
        fileLine = fileLine.strip()
        #print('total length: ' + str(len(sequence)))    
    
    print('total length: ' + str(len(sequence)))   
    start = 0 
    for i in range(len(sequence)):
        
        if i > 0 and i < len(sequence)-1:
            if sequence[i-1] != '-':
                if sequence[i] == '-':
                    start = i
            if sequence[i+1] != '-':
                if sequence[i] == '-':
                    end = i
                    f.write("\t".join((str(start),str(end),"\n")))
    
    sequence = str()
                    

