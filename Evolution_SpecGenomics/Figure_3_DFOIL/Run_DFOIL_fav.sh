#!/bin/bash
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=6:00:00
#SBATCH --job-name=DFOILfav
#SBATCH --error=DFOILfav
#SBATCH --output=DFOILfav
#SBATCH --partition=mpi
#SBATCH --ntasks=1

module load python2

ID=$((SLURM_ARRAY_TASK_ID -1))

#for i in {1..21}; 
#i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
#python ~/work/scripts/genomics_general-master/filterGenotypes.py -i /work/rpapa/share/calls_erato/chromosomes/Herato_040518.Herato$(echo "${i[ID]}").calls.gz -of diplo -s himera_001,himera_002,himera_003,himera_006,himera_030,BC2565,BC2566,BC2567,BC2570,BC2563,BC2577,BC2578,BC2579,GS020redo,GS021redo,NCS_1671,NCS_1672,NCS_1673,NCS_1674,NCS_1675,cyrCAM040545,cyrCAM040578,cyrCAM040584,cyrCAM040673,cyrCAM040692,cyrCAM040695,cyrCAM040842,cyrCAM040853,cyrCAM040860,cyrCAM040861,cyrbia_004,cyrbia_005,cyrbia_023,cyrbia_024,BC2635,BC2637,BC2638,BC2639,GS012redo,NCS_0471,NCS_0473,NCS_0476,NCS_0478,NCS_0479,NCS_2554,NCS_2555,hermathena_13,hermathena_14,hermathena_15 --minCalls 5 --minAlleles 2 | python ~/work/scripts/genomics_general-master/genoToSeq.py --format fasta --windType coordinate --windSize 50000 --stepSize 50000 --minSites 1 --separateFiles --seqFile ~/work/DFOIL_out/all --mode windows;


###########
#him_emm 560
#P1=(himera_001 himera_002 himera_003 himera_006 himera_030)
#P2=(BC2565 BC2566 BC2567 BC2570)
#P3=(BC2563 BC2577 BC2578 BC2579)
#P4=(GS020redo GS021redo NCS_1671 NCS_1672 NCS_1673 NCS_1674 NCS_1675)

#him_fav 640
P1=(himera_001 himera_002 himera_003 himera_006 himera_030)
P2=(BC2565 BC2566 BC2567 BC2570)
P3=(BC2635 BC2637 BC2638 BC2639)
P4=(GS012redo NCS_0471 NCS_0473 NCS_0476 NCS_0478 NCS_0479 NCS_2554 NCS_2555)

#him_cyr
#P1=(himera_001 himera_002 himera_003 himera_006 himera_030)
#P1=(himera_003 himera_006 himera_030)
#P2=(BC2565 BC2566 BC2567 BC2570)
#P3=(cyrbia_004 cyrbia_005 cyrbia_023 cyrbia_024)
#P4=(cyrCAM040545 cyrCAM040578 cyrCAM040584 cyrCAM040673 cyrCAM040692 cyrCAM040695 cyrCAM040842 cyrCAM040853 cyrCAM040860 cyrCAM040861)


# create list with combinations
LIST=()
for a in ${P1[@]}; do for b in ${P2[@]}; do for c in ${P3[@]}; do for d in ${P4[@]}; do 
LIST=("${LIST[@]}" "$a,$b,$c,$d,hermathena_13"); 
LISTa=("${LISTa[@]}" "$a"); 
LISTb=("${LISTb[@]}" "$b"); 
LISTc=("${LISTc[@]}" "$c"); 
LISTd=("${LISTd[@]}" "$d"); 
done; done; done; done

#if [ ! -f ~/work/DFOIL_out/him_fav_DFOIL/$(echo "${LISTa[ID]}")\_$(echo "${LISTb[ID]}")\_$(echo "${LISTc[ID]}")\_$(echo "${LISTd[ID]}")\_hermathena_13_DFOIL ]; 
#then

for f in ~/work/DFOIL_out/all*; 

do 
fbname=$(basename "$f" .fa | cut -f2 -d'.')
echo grepping file $fbname $(echo "${LISTa[ID]}") $(echo "${LISTb[ID]}") $(echo "${LISTc[ID]}") $(echo "${LISTd[ID]}")

grep --no-group-separator -A 1 -e $(echo "${LISTa[ID]}") -e $(echo "${LISTb[ID]}") -e $(echo "${LISTc[ID]}") -e $(echo "${LISTd[ID]}") -e 'hermathena_13' $f > ~/work/DFOIL_out/him_fav/$fbname\_$(echo "${LISTa[ID]}")\_$(echo "${LISTb[ID]}")\_$(echo "${LISTc[ID]}")\_$(echo "${LISTd[ID]}")\_hermathena_13.fa

done;

echo creating conts file

FILES=""
for f in ~/work/DFOIL_out/him_fav/*$(echo "${LISTa[ID]}")\_$(echo "${LISTb[ID]}")\_$(echo "${LISTc[ID]}")\_$(echo "${LISTd[ID]}")\_hermathena_13*; do FILES="$FILES $f"; done

python ~/work/Programs/dfoil/fasta2dfoil.py $FILES -o ~/work/DFOIL_out/him_fav/$(echo "${LISTa[ID]}")\_$(echo "${LISTb[ID]}")\_$(echo "${LISTc[ID]}")\_$(echo "${LISTd[ID]}")\_hermathena_13_count -n $(echo "${LIST[ID]}")

echo calculating DFOIL
python ~/work/Programs/dfoil/dfoil.py --infile ~/work/DFOIL_out/him_fav/$(echo "${LISTa[ID]}")\_$(echo "${LISTb[ID]}")\_$(echo "${LISTc[ID]}")\_$(echo "${LISTd[ID]}")\_hermathena_13_count --out ~/work/DFOIL_out/him_fav/$(echo "${LISTa[ID]}")\_$(echo "${LISTb[ID]}")\_$(echo "${LISTc[ID]}")\_$(echo "${LISTd[ID]}")\_hermathena_13_DFOIL

echo removing fasta files
for f in ~/work/DFOIL_out/him_fav/*$(echo "${LISTa[ID]}")\_$(echo "${LISTb[ID]}")\_$(echo "${LISTc[ID]}")\_$(echo "${LISTd[ID]}")\_hermathena_13.fa; do rm $f; done

#fi


############




# create list with combinations
#LIST=()
#for a in ${P1[@]}; do for b in ${P2[@]}; do for c in ${P3[@]}; do for d in ${P4[@]}; do LIST=("${LIST[@]}" "$a,$b,$c,$d,hermathena_13"); done; done; done; done


#for i in {1..21}; 
#do python ~/work/scripts/genomics_general-master/filterGenotypes.py -i /work/rpapa/share/calls_erato/chromosomes/Herato_040518.Herato$i.calls.gz -of diplo -s $(echo "${LIST[ID]}") --minCalls 5 --minAlleles 2 | python ~/work/scripts/genomics_general-master/genoToSeq.py --format fasta --windType coordinate --windSize 50000 --stepSize 50000 --minSites 1 --separateFiles --seqFile ~/work/DFOIL_out/him_fav/$(echo "${LIST[ID]}") --mode windows;
#done

#FILES=""
#for f in ~/work/DFOIL_out/him_fav/$(echo "${LIST[ID]}")*; do FILES="$FILES $f"; done

#python Programs/dfoil/fasta2dfoil.py $FILES -o ~/work/DFOIL_out/him_fav/$(echo "${LIST[ID]}")_count -n $(echo "${LIST[ID]}")

#python Programs/dfoil/dfoil.py --infile ~/work/DFOIL_out/him_fav/$(echo "${LIST[ID]}")_count --out ~/work/DFOIL_out/him_fav/$(echo "${LIST[ID]}")_DFOIL

#for x in $FILES; do rm $x; done


###
#for a in ${P1[@]}; 
#do for b in ${P2[@]}; 
#do for c in ${P3[@]}; 
#do for d in ${P4[@]};
#do for f in ~/work/DFOIL_out/all*; 

#do 
#fbname=$(basename "$f" .fa | cut -f2 -d'.')
#echo grepping file $fbname $a $b $c $d

#grep --no-group-separator -A 1 -e $a -e $b -e $c -e $d -e 'hermathena_13' $f > ~/work/DFOIL_out/him_fav/$fbname\_$a\_$b\_$c\_$d\_hermathena_13.fa

#done;

#echo creating conts file

#FILES=""
#for f in ~/work/DFOIL_out/him_fav/*$a\_$b\_$c\_$d\_hermathena_13*; do FILES="$FILES $f"; done
#python ~/work/Programs/dfoil/fasta2dfoil.py $FILES -o ~/work/DFOIL_out/him_fav/$a\_$b\_$c\_$d\_hermathena_13_count -n $a,$b,$c,$d,hermathena_13

#echo calculating DFOIL
#python ~/work/Programs/dfoil/dfoil.py --infile ~/work/DFOIL_out/him_fav/$a\_$b\_$c\_$d\_hermathena_13_count --out ~/work/DFOIL_out/him_fav/$a\_$b\_$c\_$d\_hermathena_13_DFOIL

#echo removing fasta files
#for f in ~/work/DFOIL_out/him_fav/*$a\_$b\_$c\_$d\_hermathena_13.fa; do rm $f; done


#done; done; done; done

###


