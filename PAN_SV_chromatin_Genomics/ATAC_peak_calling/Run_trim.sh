#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=8:00:00
#SBATCH --job-name=trimm
#SBATCH --error=trimm
#SBATCH --output=trimm
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load trimmomatic
module load jdk

indID=$((SLURM_ARRAY_TASK_ID -1))

#10
#sampleList1=(BR10_Demophoon_Brain BR12_Demophoon_Brain BR15_Demophoon_Brain CH_head2 E4-Head BR11_Rosina_Brain BR14_Demophoon_Brain BR5_Rosina_Brain CH_head6 M4-Head)
#6
sampleList1=(CH_FW2 CH_FW3 CH_FW7 CH_HW2 CH_HW3 CH_HW7)

#java -jar /work/rpapa/share/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 /work/rpapa/share/SVB_Maylin/ATAC/raw/$(echo "${sampleList1[indID]}")_R1.fastq.gz /work/rpapa/share/SVB_Maylin/ATAC/raw/$(echo "${sampleList1[indID]}")_R2.fastq.gz \
#/work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_R1.fastq.gz /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_unpair_R1.fastq.gz \
#/work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_R2.fastq.gz /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_unpair_R2.fastq.gz \
#ILLUMINACLIP:/work/rpapa/share/scripts/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:32

java -jar /work/rpapa/share/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 /work/rpapa/share/reads_ATACseq/$(echo "${sampleList1[indID]}")_R1.fastq.gz /work/rpapa/share/reads_ATACseq/$(echo "${sampleList1[indID]}")_R2.fastq.gz \
/work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_R1.fastq.gz /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_unpair_R1.fastq.gz \
/work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_R2.fastq.gz /work/rpapa/share/SVB_Maylin/ATAC/trimmed/$(echo "${sampleList1[indID]}")_trim_unpair_R2.fastq.gz \
ILLUMINACLIP:/work/rpapa/share/scripts/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:32

