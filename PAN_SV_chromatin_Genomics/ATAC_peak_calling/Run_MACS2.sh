#!/bin/bash
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=36:00:00
#SBATCH --job-name=macs2
#SBATCH --error=macs2
#SBATCH --output=macs2
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load python3/3.6.10

ID=$((SLURM_ARRAY_TASK_ID -1))

#5
#samples=(BR14_Demophoon_Brain BR12_Demophoon_Brain) # BR10_Demophoon_Brain E4-Head BR15_Demophoon_Brain)
#3
samples=(BR11_Rosina_Brain M4-Head) # BR5_Rosina_Brain)
#2
#samples=(CH_head2 CH_head6)
#6
samples=(CH_FW2 CH_FW3 CH_FW7 CH_HW2 CH_HW3 CH_HW7)

#macs2 callpeak \
#-t /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_Herato.trim.filtered.sorted.nd.bam \
#-n $(echo "${samples[ID]}") --outdir /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/ \
#-f BAMPE -g 382844248 --nomodel --shift -100 --extsize 200

#macs2 callpeak \
#-t /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_Hmel2.trim.filtered.sorted.nd.bam \
#-n $(echo "${samples[ID]}") --outdir /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/ \
#-f BAMPE -g 275198613 --nomodel --shift -100 --extsize 200

macs2 callpeak \
-t /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_Hcha1.trim.filtered.sorted.nd.bam \
-n $(echo "${samples[ID]}")_Hcha1 --outdir /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/ \
-f BAMPE -g 355179777 --nomodel --shift -100 --extsize 200

#macs2 callpeak \
#-t /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_Hcha2.trim.filtered.sorted.nd.bam \
#-n $(echo "${samples[ID]}")_Hcha2 --outdir /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/ \
#-f BAMPE -g 357058320 --nomodel --shift -100 --extsize 200
