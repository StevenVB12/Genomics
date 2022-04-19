#!/bin/bash
#SBATCH --mem-per-cpu=24gb
#SBATCH --time=48:00:00
#SBATCH --job-name=fasta
#SBATCH --error=fasta
#SBATCH --output=fasta
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load python2

python /work/rpapa/share/SVB_Maylin/SeqSeqPan_erato_melp_char/seq-seq-pan_toFasta.py -I /work/rpapa/share/SVB_Maylin/SeqSeqPan_erato_melp_char/SeqSeqPan_erato_melp_char_noNewline.xmfa -g 1,2,3,4
