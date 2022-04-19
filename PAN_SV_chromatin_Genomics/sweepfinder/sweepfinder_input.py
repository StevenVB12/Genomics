###***********************************************************************************************
###
###modified from Simon Martin's createSweepFinderInput.py
###
### VERSION 0.1
### 15.05.2018
###***********************************************************************************************
###
### This script generates the input for SweepFinder.
### An outgroup sequence/population can be defined as a reference for the ancestral state,
###
### Only biallelic sites will be considered.
### If multiple out-groups exist and they are not fixed for a single allele,
### a consensus can be taken. 
### For each population specified, and each scaffold indicated, a unique input file will be created
### with the name POP_SCAF.sweepfinder.input
###
###***********************************************************************************************
###
### Input "calls file" format:
### scaffold  position  ind1  ind2  ind3  etc...
### scf1  1 T T C Y
###
### NOTE a header row with unique names for all individuals is essential.
###
###
###***********************************************************************************************
###
### USAGE
###
### python createSweepFinderInput.py -i <input file> -p <population_string> -O <out-group population name> -c <chromosome/scaffold name> [or --file <file containing chrom or scaf names]
##
### the populations are specified in a single string as follows:
### -p "pop1Name[ind1,ind2,ind3,ind4];pop2Name[ind1,ind2,ind3,ind4]" (quotation marks must be present to avoid conflict with unix)
### the pout-group must match one of the population names. e.g. -O pop2Name
###
### scaffolds must be spaced by commas only (or one per line if a file is provided) and when given on command line add quotes
###
### USE ONLY ONE FOCAL POP with or without OUTGROUP
###
### --addMono adds sites monomorphic in the ingroup but substitutions compared to outgroup (=substitutions), e.g. x=10, n=10, folded=0
### --addMonoall adds sites that are fixed for the same allele  in both in- and outgroup and polarized, e.g. x=0, n=10, folded=0, and if there are not enough outgroup samples for polarization but all alleles present are identical, it adds unpolarized monomorphic sites, i.e. x=0, n=10 folded=1
### --polarized outputs only sites that can be polarized
###
### in order to output all sites monomorphic for the ingroup (ancestral monomorphic and substitutions and unpolarized monomorphic sites), both -addMono and --addMonoall need to be set
###
###
### for polarized polymorphic sites, substitutions and ancestral monomorphic: --addMono, --addMonoall, --polarized
### for polarized polymorphic sites and substitutions: --addMono, --polarized
### for polarized polymorphic sites and substitutions where the ancestral state is randomly selected from all OG alleles considering : --addMono, --polarized, --randomOGbase
###
###***********************************************************************************************
###***********************************************************************************************


import sys, random

### Functions
## for parsing the pops and individuals
def get_intv(string,borders = "()",inc = False):
  if len(borders) != 2:
    print "WARNING: borders must contain two characters"
  starts = []
  ends = []
  output = []
  for x in range(len(string)):
    if string[x] == borders[0]:
      starts.append(x)
    if string[x] == borders[1]:
      ends.append(x+1)
  if len(starts) <= len(ends):
    for n in range(len(starts)):
      if inc:
        output.append(string[starts[n]:ends[n]])
      else:
        output.append(string[starts[n]+1:ends[n]-1])
  else:
    for n in range(len(ends)):
      if inc:
        output.append(string[starts[n]:ends[n]])
      else:
        output.append(string[starts[n]+1:ends[n]-1])
  return output

##translate diplo calls
def haplo(calls):
  output = []
  for call in calls:
    if call in "ACGTN":
      output.append(call)
      output.append(call)
    elif call == "K":
      output.append("G")
      output.append("T")
    elif call == "M":
      output.append("A")
      output.append("C")
    elif call == "R":
      output.append("A")
      output.append("G")
    elif call == "S":
      output.append("C")
      output.append("G")
    elif call == "W":
      output.append("A")
      output.append("T")
    elif call == "Y":
      output.append("C")
      output.append("T")
    else:
      print "WARNING", call, "is not recognised as a valid base or ambiguous base"
      output.append("N")
      output.append("N")
  return output

##get option value from input (check what this does)
def getOptionValue(option):
  optionPos = [i for i,j in enumerate(sys.argv) if j == option][0]
  optionValue = sys.argv[optionPos + 1]
  return optionValue

##return unique elements
def unique(things):
  output = list(set(things)) # set returns a unordered collection of unique elements which are here coverted into a list
  output.sort() # orders the input
  return output

## exclude elements that equal x
def exclude(things, x):
  output = [i for i in things if i != x]
  return(output)

## returns unique alleles from the haplo function output
def uniqueAlleles(bases):
  haploBases = haplo(bases)
  output = unique([i for i in haploBases if i in "ACGT"])
  return output

## calculate base frequency
def baseFreq(bases,base):
  haploBases = [i for i in haplo(bases) if i in "ACGT"]
  freq = (float(haploBases.count(base))) / len(haploBases)
  return freq

## output most common element
def mostCommon(things):
  output = []
  counts = []
  uniqueThings = unique(things)
  for thing in uniqueThings:
    counts.append(things.count(thing))
  maxCount = max(counts)
  for n in range(len(counts)):
    if counts[n] == maxCount:
      output.append(uniqueThings[n])
  return output

## return consensus base
def consensus(bases):
  haploBases = [i for i in haplo(bases) if i in "ACGT"]
  major = mostCommon(haploBases)
  return major[0] ## here I'll probably need a random number generator somewhere (if major has more than two elements randomly select first or second)

## draw from all alleles (for drawing ancestral state from outgroup)
def drawFromOG(bases):
  haploBases = [i for i in haplo(bases) if i in "ACGT"]
  randBase = random.choice(haploBases)
  #print "all alleles (except for N's): ",haploBases
  #print "randomly selected base: ",randBase
  return randBase

### get files

if "-i" in sys.argv:
  fileName = getOptionValue("-i")
  file = gzip.open(fileName , "r") if fileName .endswith(".gz") else open(fileName , "r")
else:
  file = sys.stdin

#file = open(fileName, "rU")
line = file.readline()
names = line.split()
line= file.readline()


if "-o" in sys.argv:
  outPref = getOptionValue("-o")
else:
  outPref = "out"

# read the population string - for dealing with multiple populations, the script needs to be re-worked, works fine for one pop and 1 OG
if "-p" in sys.argv:
  popString = getOptionValue("-p")
else:
  print "\nplease specify populations using -p\n"
  sys.exit()

if "-O" in sys.argv:
  outGroup = getOptionValue("-O")
  OG = True
else:
  OG = False



if "-M" in sys.argv:
  minimum = int(getOptionValue("-M"))
else:
  minimum = 5

if "-OM" in sys.argv:
  OGminimum = int(getOptionValue("-OM"))
else:
  OGminimum = 1



pops = []
#for each population, store the name and individual names
for popData in popString.strip("\"").split(";"): # strip the quotes (necessary for linux command line) and split at ";"
  currentPop = popData.split("[")[0]
  pops.append(currentPop)
  vars()[currentPop + "Inds"] = get_intv(popData,"[]")[0].split(",") #vars() creates a name for an object consisting of a variable and a string here (useful to create a variablename)
  for ind in vars()[currentPop + "Inds"]:
    if ind not in names:
      print ind, "not found in header line."
      sys.exit()

if OG:
  if outGroup not in pops:
    print "the specified outgroup, ", outGroup, ", was not a specified population."
    sys.exit()
  pops = [pop for pop in pops if pop!=outGroup]



if "-c" in sys.argv:
  scafs = getOptionValue("-c").split(",")
  allAsOne = False
else:
  print "\nScaffold or chromosome  not specified - outputting all sites to one file\n"
  allAsOne = True



if "--OGconsensus" in sys.argv:
  outgroupConsensus = True
else:
  outgroupConsensus = False

if "--addMono" in sys.argv:
  addMono = True
else:
  addMono = False

if "--addMonoall" in sys.argv:
  addMonoall = True
else:
  addMonoall = False

if "--polarized" in sys.argv:
  polarized = True
else:
  polarized = False

if "--randomOGbase" in sys.argv:
  randomOGbase = True
else:
  randomOGbase = False

if "--indStart" in sys.argv:
  indStart = int(getOptionValue("--indStart"))
else:
  indStart = 2

if "--chromCol" in sys.argv:
  chromCol = int(getOptionValue("--chromCol"))
else:
  chromCol = 0


###************************************************************************************************

### output files. For each population, for each scaffold, we create an output file.

for pop in pops:
  if allAsOne: # when no scaffolds are indicated
    vars()[pop] = open(outPref + "." + pop + ".sweepfinder.input", "w") #creates file
    vars()[pop].write("position\tx\tn\tfolded\n") # writes header to file
  else:
    for scaf in scafs: # when scafolds are indicated
      vars()[pop + scaf] = open(outPref + "." + pop + "_" + scaf + ".sweepfinder.input", "w")
      vars()[pop + scaf].write("position\tx\tn\tfolded\n")

p = 1

### for each line, check if its a scaf we want
while len(line) > 1:
  objects = line.split()
  if allAsOne or objects[chromCol] in scafs:
    #first check that the allele is biallelic over all samples (we don't stuff around with tri-allelic crap)
    allCalls = objects[indStart:]
    allAlleles = uniqueAlleles(allCalls) #check how many unique alleles are in a line (=outgroup + ingroups) - check this when correcting multiple pops error
    #print allCalls
    #print allAlleles
    if len(allAlleles) == 2:
      #then start by getting the ancestral state, if there's an outgroup
      if OG:  # when an outgroup is given
        ogCalls = [] # create empty list for Outgroup calls
        for ind in vars()[outGroup + "Inds"]:  # loop over all outgroup individuals
          ogCalls.append(objects[names.index(ind)]) # append diplo call of each individual in outgroup in this line
          #print "outgroup calls1:", ogCalls
        ogCalls = [i for i in ogCalls if i != "N"] #remove "N"s
        #print "outgroup calls:", ogCalls
        if len(ogCalls) >= OGminimum:           # check for sufficient data in OG saamples
          ogAlleles = uniqueAlleles(ogCalls)
          if len(ogAlleles) == 1:               # if OG is monomorphic (fixed), site is polarized, ancState is the fixed OG allele
            ancState = ogAlleles[0]
            folded = 0
          elif outgroupConsensus:               # if outgroupConsensus is specified and sufficient samples have data, site is polarized, ancState is outgroup consensus
            ancState = consensus(ogCalls)
            print(ancState)
            folded = 0
          elif randomOGbase == True:		# if randomOGbase is specified and sufficient samples have data, the site is polarized and ancStat is a randomly chosen allele from the outgroup with chanced of sampling the allele being proportional to its frequency             
            ancState = drawFromOG(ogCalls)
            folded = 0	
          else:                                 # if outgroup is polymorphic and outgroupConsensus not specified, sites are not polarized, ancState is consensus of allCalls, this may include sites that are monomorphic in the ingroup!
            ancState = consensus(allCalls)
            folded = 1
        else:                                   # if number of outgroup samples is below OGminimim, sites are not polarized, ancState is consensus of allCalls, this may also include sites that are monomorphic in the ingroup
          ancState = consensus(allCalls)
          folded = 1
      else:                                     # if no outgroup is given, sites are not polarized, ancState is consensus of all Calls
        ancState = consensus(allCalls)
        folded = 1
      derState = exclude(allAlleles, ancState)[0] # exclude ancestral alleles from allAlleles and pick the left derived allele
      #print "ancestral state:", ancState
      #print "derived state:", derState
      #for each pop, check sufficient sample size, then output
      for pop in pops:
        calls = []
        for ind in vars()[pop + "Inds"]:
          calls.append(objects[names.index(ind)])
        calls = [i for i in calls if i != "N"]
        if len(calls) >= minimum:
          bases = haplo(calls)
          derCount = len([i for i in bases if i == derState])
          if addMono == True or (derCount > 0 and derCount < len(bases)):
            if polarized == True and folded == 0 and derCount > 0:
                if allAsOne:
                  vars()[pop].write(str(p) + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")
                else:
                  vars()[pop + objects[chromCol]].write(objects[1] + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")
            elif polarized == False:
                if allAsOne:
                  vars()[pop].write(str(p) + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")
                else:
                  vars()[pop + objects[chromCol]].write(objects[1] + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")

    elif addMonoall == True and len(allAlleles) == 1: #returns all fixed sites (no matter how many OGminimum samples were assessed)
      if OG:  # when an outgroup is given
        ogCalls = [] # create empty list for Outgroup calls
        for ind in vars()[outGroup + "Inds"]:  # loop over all outgroup individuals
          ogCalls.append(objects[names.index(ind)]) # append diplo call of each individual in outgroup in this line
          #print "outgroup calls1:", ogCalls
        ogCalls = [i for i in ogCalls if i != "N"] #remove "N"s
        #print "outgroup calls:", ogCalls
        if len(ogCalls) >= OGminimum:           # check for sufficient data in OG samples and if that's ok this fixed and polarized site is added
            folded = 0
        else:
            folded = 1                          # if not sufficient OG samples, the fixes site is added as not polarized
        for pop in pops:
            calls = []
            for ind in vars()[pop + "Inds"]:
                calls.append(objects[names.index(ind)])
            calls = [i for i in calls if i != "N"]
            if len(calls) >= minimum:
                bases = calls*2
                derCount = 0
                if polarized == True and folded == 0:
                    if allAsOne:
                        vars()[pop].write(str(p) + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")
                    else:
                        vars()[pop + objects[chromCol]].write(objects[1] + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")
                elif polarized == False:
                    if allAsOne:
                        vars()[pop].write(str(p) + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")
                    else:
                        vars()[pop + objects[chromCol]].write(objects[1] + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")

    # elif addMono and len(allAlleles) == 1: #returns all fixed sites (no matter how many OGminimum samples were assessed) and returns them as unpolarized
    #       folded = 1
    #       for pop in pops:
    #         calls = []
    #         for ind in vars()[pop + "Inds"]:
    #           calls.append(objects[names.index(ind)])
    #         calls = [i for i in calls if i != "N"]
    #         if len(calls) >= minimum:
    #           bases = calls*2
    #           derCount = 0
    #           if allAsOne:
    #             vars()[pop].write(str(p) + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")
    #           else:
    #             vars()[pop + objects[chromCol]].write(objects[1] + "\t" + str(derCount) + "\t" + str(len(bases)) + "\t" + str(folded) + "\n")

    p += 1
  line = file.readline()

file.close