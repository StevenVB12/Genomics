#!/bin/bash
#SBATCH --mem-per-cpu=24gb
#SBATCH --time=8:00:00
#SBATCH --job-name=trimm
#SBATCH --error=trimm
#SBATCH --output=trimm
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=8

module load bedtools

ID=$((SLURM_ARRAY_TASK_ID -1))

#5
#samples=(BR10_Demophoon_Brain BR12_Demophoon_Brain BR15_Demophoon_Brain E4-Head BR14_Demophoon_Brain)
#3
#samples=(BR11_Rosina_Brain BR5_Rosina_Brain M4-Head)
#2
samples=(CH_head2 CH_head6)

#REF=/work/rpapa/share/REF/H_erato_dem/Herato_final_bowtie
#REF=/work/rpapa/share/REF/Hmel2/Hmel2_bowtie
#REF=/work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1_bowtie
REF=/work/rpapa/share/REF/H_charithonia_PR/H_charithonia_10X.2_bowtie

#REFNAME=Herato
#REFNAME=Hmel2
#REFNAME=Hcha1
REFNAME=Hcha2

#SIZES=/work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes
#SIZES=/work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes
#SIZES=/work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1.fasta.sizes
SIZES=/work/rpapa/share/REF/H_charithonia_PR/H_charithonia_10X.2.fasta.sizes



#### get coverage data

/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov \
-ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.bam -bg \
> /work/rpapa/share/SVB_Maylin/ATAC/bw/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.bdg 

LC_COLLATE=C sort -k1,1 -k2,2n /work/rpapa/share/SVB_Maylin/ATAC/bw/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.bdg \
> /work/rpapa/share/SVB_Maylin/ATAC/bw/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.collate.bdg 

rm /work/rpapa/share/SVB_Maylin/ATAC/bw/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.bdg 

/work/rpapa/sbelleghem/scripts/bedGraphToBigWig \
/work/rpapa/share/SVB_Maylin/ATAC/bw/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.collate.bdg \
$SIZES \
/work/rpapa/share/SVB_Maylin/ATAC/bw/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.collate.bw


