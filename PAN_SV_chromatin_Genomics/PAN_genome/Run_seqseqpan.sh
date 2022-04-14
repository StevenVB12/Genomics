#!/bin/bash
#SBATCH --mem-per-cpu=120gb
#SBATCH --time=144:00:00
#SBATCH --job-name=seqpan
#SBATCH --error=seqpan
#SBATCH --output=seqpan
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=8


module load seq-seq-pan/current

cd /work/rpapa/sbelleghem/SeqSeqPan_erato_melp_char/

seq-seq-pan-wga --config genomefile=/work/rpapa/sbelleghem/SeqSeqPan_erato_melp_char/genome_list.txt outfilename=/work/rpapa/sbelleghem/SeqSeqPan_erato_melp_char/SeqSeqPan_erato_melp_char