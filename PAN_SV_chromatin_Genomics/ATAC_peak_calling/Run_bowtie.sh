#!/bin/bash
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=8:00:00
#SBATCH --job-name=trimm
#SBATCH --error=trimm
#SBATCH --output=trimm
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=8

module load bowtie2
module load samtools
module load jdk

ID=$((SLURM_ARRAY_TASK_ID -1))

#5
#samples=(BR10_Demophoon_Brain BR12_Demophoon_Brain BR15_Demophoon_Brain E4-Head BR14_Demophoon_Brain)
#3
samples=(BR11_Rosina_Brain BR5_Rosina_Brain M4-Head)
#2
#samples=(CH_head2 CH_head6)
#6
samples=(CH_FW2 CH_FW3 CH_FW7 CH_HW2 CH_HW3 CH_HW7)

#REF=/work/rpapa/share/REF/H_erato_dem/Herato_final_bowtie
#REF=/work/rpapa/share/REF/Hmel2/Hmel2_bowtie
REF=/work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1_bowtie
#REF=/work/rpapa/share/REF/H_charithonia_PR/H_charithonia_10X.2_bowtie

#REFNAME=Herato
#REFNAME=Hmel2
REFNAME=Hcha1
#REFNAME=Hcha2

#samtools faidx /work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1.fasta; cut -f1,2 /work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1.fasta.fai > /work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1.fasta.sizes
#samtools faidx /work/rpapa/share/REF/H_charithonia_PR/H_charithonia_10X.2.fasta; cut -f1,2 /work/rpapa/share/REF/H_charithonia_PR/H_charithonia_10X.2.fasta.fai > /work/rpapa/share/REF/H_charithonia_PR/H_charithonia_10X.2.fasta.sizes


#SIZES=/work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes
#SIZES=/work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes
SIZES=/work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1.fasta.sizes
#SIZES=/work/rpapa/share/REF/H_charithonia_PR/H_charithonia_10X.2.fasta.sizes

#### Run bowtie

gunzip /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${samples[ID]}")_trim_R1.fastq.gz
gunzip /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${samples[ID]}")_trim_R2.fastq.gz

bowtie2 -t -k 2 -p 8 --local -x $REF \
-1 /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${samples[ID]}")_trim_R1.fastq \
-2 /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${samples[ID]}")_trim_R2.fastq |\
samtools view -bS - > /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.bam

gzip /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${samples[ID]}")_trim_R1.fastq
gzip /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${samples[ID]}")_trim_R2.fastq


#### filter bam files

samtools view -f 0x02 -q 20 -b /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.bam > /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.bam

samtools sort /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.bam /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted

java -jar /work/rpapa/sbelleghem/Programs/picard-tools-2.5.0/picard.jar MarkDuplicates \
I=/work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.bam \
O=/work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.bam \
Remove_Duplicates=true  M=/work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_dup_metrics.txt ASSUME_SORTED=true

rm /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.bam
rm /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.bam

