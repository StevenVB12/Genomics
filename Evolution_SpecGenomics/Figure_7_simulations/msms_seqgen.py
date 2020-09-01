#!/usr/bin/env python

### September 2017 ### Borrowed from Simon Martin ###
### This script evolved from one that ran simulations and analysed the data.
### now it is a generic script for running msms and seq-gen and outputing the seqs
### However, it does include functions to converte these into an alignment using my genomics library
from __future__ import division
import subprocess, tempfile, sys, os
import argparse
import egglib
import genomics


def run_msms(Nsam, Nreps=1, pops=None, theta = None, joins=[], splits=[], migs=[], migChanges=[],
             trees = True, threads = 1, 
             thetaIsProp = False, length = 100, recRate = 0, Ne=None,
             selection = False, s=0, h=0, pos=0.5, selStartTime = None, selStartFreqs=[0], Smu = None,
             selChanges=[], alpha=0):
    
    msms_par = make_msms_params(Nsam, Nreps, pops, theta, joins, splits, migs, migChanges,
                     trees, threads,  
                     thetaIsProp, length, recRate, Ne,
                     selection, s, h, pos, selStartTime, selStartFreqs, Smu, selChanges, alpha)
    
    sys.stderr.write(msms_par + "\n")
    
    result = subprocess.check_output(["/work/rpapa/sbelleghem/Programs/msms/bin/msms", msms_par]).split("\n")
    
    return result


def run_seqGen(trees,length,seqGenModel="HKY",seqGenScaler=0.01):
        
    with tempfile.NamedTemporaryFile(mode="w",delete=False) as tempTrees:tempTrees.write("\n".join(trees))
    
    #run seqgen, suppressing stderr
    with open(os.devnull, 'w') as devnull:
        seqStr = subprocess.check_output(["/work/rpapa/sbelleghem/Programs/Seq-Gen.v1.3.3/source/seq-gen", "-m" + seqGenModel, "-l", str(length), "-s", str(seqGenScaler), "-p" + str(len(trees)), tempTrees.name], stderr=devnull)
    
    os.remove(tempTrees.name)
    
    return seqStr




def make_msms_params(Nsam, Nreps=1, pops=None, theta = None, joins=[], splits=[], migs=[], migChanges=[],
                     trees = True, threads = 1, 
                     thetaIsProp = False, length = 100, recRate = 0, Ne=None,
                     selection = False, s=0, h=0, pos=0.5, selStartTime = None, selStartFreqs=[0], Smu = None,
                     selChanges=[], alpha=0):

    
    joinsText = " ".join(["-ej " + " ".join([str(i) for i in j]) for j in joins])
    splitsText = " ".join(["-es " + " ".join([str(i) for i in j]) for j in splits])
    migsText =  " ".join(["-m " + " ".join([str(i) for i in j]) for j in migs])
    migChangesText = " ".join(["-em " + " ".join([str(i) for i in j]) for j in migChanges])
    
    recText = " ".join(["-r",str(1.*length*recRate),str(length)])
    
    popsText = "-I " + str(len(pops)) + " " + " ".join([str(x) for x in pops]) if pops else ""
    
    if thetaIsProp: theta = theta*length
    thetaText = "-t " + str(theta) if theta else ""
    
    treesText = "-T" if trees else ""
    
    threadsText = "-threads " + str(threads)
    
    alphaText = "-G " + str(alpha)
    
    parText = " ".join([str(Nsam), str(Nreps), thetaText, treesText, popsText, joinsText, splitsText, migsText, migChangesText, recText, alphaText])
    
    if selection:
        assert Ne is not None
        #sAA_text = "-SAA " + str(s*2.*Ne)
        #sAa_text = "-SAa " + str(s*h*2.*Ne)
        posText = "-Sp " + str(pos)
        SmuText = "-Smu " + str(Smu) if Smu else ""
        selStartText = " ".join(["-SI ",str(selStartTime), str(len(selStartFreqs)), " ".join([str(x) for x in selStartFreqs])]) if selStartTime else ""
        NeText = "-N " + str(Ne)
        selChangesText = " ".join(["-Sc " + " ".join([str(i) for i in j]) for j in selChanges])
        
        #parText += " " + " ".join([sAA_text, sAa_text, posText, SmuText, selStartText, NeText, selChangesText, alphaText])
        parText += " " + " ".join([posText, SmuText, selStartText, NeText, selChangesText, alphaText])
    
    return parText


def simSeq(paramDict, returnGroupDict=False):
    
    msmsOutput = run_msms(Nsam=paramDict["Nsam"], Nreps=paramDict["Nreps"], theta=paramDict["theta"], pops=paramDict["pops"], Ne=paramDict["Ne"],
                          joins = paramDict["joins"],splits = paramDict["splits"], migs=paramDict["migs"], migChanges=paramDict["migChanges"],
                          length = paramDict["length"], recRate=paramDict["recRate"], pos=paramDict["pos"], thetaIsProp=paramDict["thetaIsProp"],
                          selection=paramDict["selection"], s = paramDict["s"], h = paramDict["h"], 
                          selStartTime=paramDict["selStartTime"], selStartFreqs=paramDict["selStartFreqs"], Smu = paramDict["Smu"], 
                          selChanges=paramDict["selChanges"], alpha=paramDict["alpha"])
    
    trees = [i for i in msmsOutput if i.endswith(";")]
    
    seqGenStr = run_seqGen(trees,length=paramDict["length"],seqGenScaler=paramDict["seqGenScaler"])
    
    if returnGroupDict:
        groupDict = dict(zip([str(x) for x in range(1,paramDict["Nsam"]+1)],
                        [x for pop in [[paramDict["popNames"][i]]*paramDict["pops"][i]
                                    for i in range(len(paramDict["pops"]))] for x in pop]))
        
        return (seqGenStr, groupDict,)

    return seqGenStr


#a function to add the necessary parser arguments. This is so that you can import this function in other scripts and it'll automatically add the required arguments
def addSimArgsToParser(parser):
    parser.add_argument("--Nsam", help="Number of samples", type=int, action = "store", required=True)
    parser.add_argument("--theta", help="theta per site", type=float, action = "store", required=True)
    parser.add_argument("--thetaIsProp", help="theta per site", action = "store_true")
    parser.add_argument("--pops", help="Samples per population", type=int, nargs="+", action = "store")
    parser.add_argument("--join", help="Join events", nargs=3, action = "append")
    parser.add_argument("--split", help="Split events", nargs=3, action = "append")
    parser.add_argument("--mig", help="Migratin rates", nargs=3, action = "append")
    parser.add_argument("--migChange", help="Migratin rate changes", nargs=4, action = "append")
    parser.add_argument("--timesInGenerations", help="Times are specified in generations (not 4Ne*Gen)", action='store_true')
    parser.add_argument("--migsInProportions", help="Migrations are specified in proportions (not 4Ne*prop)", action='store_true')
    parser.add_argument("--length", help="Sequence length", type=int, action = "store", required=True)
    parser.add_argument("--recRate", help="Recombination rate, rho = 4Ne*r", type=float, action = "store", default=0.0)
    parser.add_argument("--selection", help="Do selection?", action = "store_true")
    parser.add_argument("--Ne", help="Population size", type=float, action = "store")
    parser.add_argument("--mu", help="Mutation rate for seq-gen scaling (if not providing seqGenScaler)", type=float, action = "store")
    parser.add_argument("--seqGenScaler", help="Scaler for seq-gen", type=float, action = "store")
    parser.add_argument("--selCoefficient" , help="Selection coefficient", type=float, action = "store", default = 0.)
    parser.add_argument("--hetEffect" , help="Selection heterozygote effect", type=float, action = "store", default = 0.)
    parser.add_argument("--selPos" , help="Selected site position", type=float, action = "store", default = 0.5)
    parser.add_argument("--selStartTime" , help="Selection start time", type=float, action = "store")
    parser.add_argument("--selStartFreqs" , help="Selection start frequencies", type=float, nargs="+", action = "store", default = [0.])
    parser.add_argument("--selMu" , help="Selected allele mutation rate, 4Ne*mu", type=float, action = "store", default = [0.])
    parser.add_argument("--selChange", help="Selection changes", nargs=5, action = "append")
    parser.add_argument("--msmsThreads", help="Threads for msms", type=int, action = "store", default = 1)
    parser.add_argument("--popNames", help="Pop names (only applies if converting to Aln object)", type=str, nargs="+",action = "store")
    parser.add_argument("--alpha", help="Growth rate", type=float, action = "store", default=0.0)

def getParamDict(args):
    paramDict={'Nsam': args.Nsam,
            'Nreps': 1,
            'pops': args.pops,
            'theta': args.theta,
            'thetaIsProp': args.thetaIsProp,
            'joins': args.join if args.join else [],
            'splits': args.split if args.split else [],
            'migs': args.mig if args.mig else [],
            'migChanges': args.migChange if args.migChange else [],
            'trees': True,
            'threads': args.msmsThreads,
            'timesInGenerations': args.timesInGenerations,
            'migsInProportions': args.migsInProportions,
            'thetaIsProp': True,
            'length': args.length,
            'recRate': args.recRate,
            'Ne': args.Ne,
            'mu': args.mu,
            'seqGenScaler': args.seqGenScaler,
            'selection': args.selection,
            's': args.selCoefficient,
            'h': args.hetEffect,
            'pos': args.selPos,
            'selStartTime' : args.selStartTime,
            'selStartFreqs': args.selStartFreqs,
            'Smu': args.selMu,
            'selChanges': args.selChange if args.selChange else [],
            'alpha': args.alpha,
            }
        #convert times
    if args.timesInGenerations:
        assert args.Ne is not None
        paramDict["selStartTime"] = args.selStartTime/(4.*args.Ne)
        for j in paramDict["joins"]: j[0] = float(j[0])/(4.*args.Ne)
        for s in paramDict["splits"]: s[0] = float(s[0])/(4.*args.Ne)
        for m in paramDict["migChanges"]: m[0] = float(m[0])/(4.*args.Ne)
    if args.migsInProportions:
        assert args.Ne is not None
        for m in paramDict["migs"]: m[2] = float(m[2])*(4.*args.Ne)
        #paramDict["migs"][0][2] = paramDict["migs"][0][2]*(paramDict["joins"][1][0]/paramDict["migChanges"][0][0])
        #paramDict["migs"][1][2] = paramDict["migs"][1][2]*(paramDict["joins"][1][0]/paramDict["migChanges"][0][0])
        for m in paramDict["migChanges"]: m[3] = float(m[3])*(4.*args.Ne)
        
    if args.pops is not None:
        paramDict["popNames"] = args.popNames if args.popNames is not None else [int(i+1) for i in range(len(args.pops))]
    
    if paramDict["seqGenScaler"] is None:
        assert paramDict["Ne"] is not None and paramDict["mu"] is not None
        paramDict["seqGenScaler"] = 4*paramDict["Ne"]*paramDict["mu"]
    
    return paramDict


def dxy(align): # "align" if the egglib alignment object, this consists of sequences, sequence names and "groups". If the object contains two groups, the function will consider only the first two.
    
    # retrieve all the positions of sequences in group 1
    P1 = [i for i in range(len(align)) if int(align[i][2])==2]
    # retrieve all the positions of sequences in group 2
    P2 = [i for i in range(len(align)) if int(align[i][2])==3]
    
    pairwiseSum = 0     #total of pairwise Pis
    totalPairs = 0      #haplotype pairs considered
    
    totalCalculated = []
    for i in P1:        #for each sequence in pop1...
      for j in P2:      #for sequence in pop2...
        seqA = align[i][1]
        seqB = align[j][1]
        zippedSeqs = zip(seqA,seqB)
        diffs = sum(sA != sB for sA, sB in zippedSeqs if sA != "N" and sB != "N")
        #sites = sum(sA != "N" and sB != "N" for sA, sB in zippedSeqs)
        sites = len([site for site in zippedSeqs if site[0] != "N" and site[1] != "N"])
        parDiff = 1.0 * diffs/sites
        totalCalculated.append(parDiff)
        
    #after considering all positions for each pair of haplotypes, return the average pairwise pi
    return sum(totalCalculated) / float(len(totalCalculated))# should be set to sites when monomorphic sites are included!!!
    #return 1.0 * diffs/windSize

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outFile", help="Output csv file", action = "store")
    parser.add_argument("--Nsim", help="Number of simulations", type=int, action = "store", default=1)
    
    addSimArgsToParser(parser)
    
    args = parser.parse_args()
    
    #args = parser.parse_args("--Nsam 40 --pops 10 10 10 10 --join 0.5 2 1 --join 1. 3 1 --join 1.5 4 1 --length 1000".split())
    
    if args.outFile:
        outFile = gzip.open(args.outFile,"r") if args.outFile.endswith(".gz") else open(args.outFile,"w")
    else: outFile=sys.stdout
    
    paramDict = getParamDict(args)
    

    #paramDict={'Nsam': 16,
            #'Nreps': 1,
            #'pops': [4,4,4,4],
            #'theta': None,
            #'joins': [(0.5,2,1),(1,3,1),(1.5,4,1)],
            #'splits': [],
            #'migs': [(2,3,4.)],
            #'migChanges': [(0.1,2,3,0)],
            #'trees': True,
            #'threads': 1,
            #'timesInGenerations': False,
            #'migsInProportions': False,
            #'thetaIsProp': False,
            #'length': 1000,
            #'recRate': 0.02,
            #'Ne': 100000,
            #'selection': True,
            #'s': 0.01,
            #'h': 0.5,
            #'pos': 0.5,
            #'selStartTime' : 0.1,
            #'selStartFreqs': [0,0,0,0],
            #'Smu': 0.001,
            #'selChanges': []
            #}

    s= 0
    while s < args.Nsim:

        sys.stderr.write("Doing simulation " + str(s) + "\n")
        
        try:
            seqs = simSeq(paramDict, returnGroupDict = True)
            
            popS = args.pops
    
            lines = seqs[0].split('\n')
            lines = lines[1:popS[0]+popS[1]+popS[2]+popS[3]+1]
    
            sequencesALL = []
            
            for i in lines:
                line = i.split()
                sequencesALL.append((line[0], line[1], seqs[1][line[0]]))
            #print sequences
                
            # seqPop1 = []
            # seqPop2 = []
            # sampleNrPop1 = 0
            # sampleNrPop2 = 0
            # for x in sequences:
            #     if x[2] == paramDict["popNames"][0]:
            #         sampleNrPop1 += 1
            #         seqPop1.append(x)
            #     if x[2] == paramDict["popNames"][1]:
            #         sampleNrPop2 += 1
            #         seqPop2.append(x)
            
            sequencesS = []
            sampleNrPop = 0
            for p in sequencesALL:
                if p[2] == 3 or p[2] == 2:
                    sampleNrPop += 1
                    sequencesS.append(p)
    
            
            alignSeq = egglib.Align.create(sequencesS)
            # alignSeqPop1 = egglib.Align.create(seqPop1)
            # alignSeqPop2 = egglib.Align.create(seqPop2)  
                 
            eggStats = egglib.stats.ComputeStats()
            # eggStats1 = egglib.stats.ComputeStats()
            # eggStats2 = egglib.stats.ComputeStats()
            
            eggStats.add_stats('lseff', 'Fst', 'Pi','S')
            # eggStats1.add_stats('lseff', 'Fst', 'Pi', 'D','S')
            # eggStats2.add_stats('lseff', 'Fst', 'Pi', 'D','S')
            
            structure = egglib.stats.get_structure(alignSeq, lvl_pop=0)
            # print structure.as_dict()
            
            stats = eggStats.process_align(alignSeq, struct=structure)
            # stats1 = eggStats.process_align(alignSeqPop1)
            # stats2 = eggStats.process_align(alignSeqPop2)
            # print stats
            
            if s == 0:
                outHeaderList = ['rec','pos','m23','mchange','split','selStart','Fst', 'Pi', 'S', 'Dxy', 'D', 'fd', 'ABBA', 'BABA']  
                outFile.write('\t'.join([str(i) for i in outHeaderList]) + '\n')
                
                # outFile = gzip.open(args.outFile,"r") if args.outFile.endswith(".gz") else open(args.outFile,"a")
            
            Fst = stats['Fst']
            Pi = round(int(stats['Pi'])/int(stats['lseff']),6)
            # Pi1 = round(int(stats1['Pi'])/int(stats1['lseff']),6)
            # Pi2 = round(int(stats2['Pi'])/int(stats2['lseff']),6)
            # D = stats['D']
            # D1 = stats1['D']
            # D2 = stats2['D']
            S = stats['S']
            # S1 = stats1['S']
            # S2 = stats2['S']
            Dxy = round(dxy(sequencesS),6)
            
            
            sequencesRAW = []
            seqNames = []
            groupNames = []
            for i in sequencesALL:
                seqNames.append(i[0])
                sequencesRAW.append(i[1])
                groupNames.append(i[2])
                
            Aln = genomics.Alignment(sequences = sequencesRAW, names = seqNames, groups = groupNames, positions = None, sampleNames = groupNames)
            
            P1=4
            P2=3
            P3=2
            O=1
            
            statsDict = genomics.fourPop(Aln, P1, P2, P3, O, minData = 0, polarize = False, fixed = False)
            
            fd = statsDict['fd']
            D = statsDict['D']
            ABBA= statsDict['ABBA']
            BABA= statsDict['BABA']
            outList = [args.recRate, args.selPos, args.mig[0][2], paramDict["migChanges"][0][0], paramDict["joins"][0][0], paramDict["selStartTime"], Fst, Pi, S, Dxy, D, fd, ABBA, BABA]
            
    
            outFile.write('\t'.join([str(i) for i in outList]) + '\n')
            
            s+=1
            
        except:
            print "Lets try again"
            continue
            
            # outFile.write(seqs[0])

    outFile.close()

