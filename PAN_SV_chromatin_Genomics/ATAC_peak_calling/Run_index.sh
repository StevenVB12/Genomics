#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=8:00:00
#SBATCH --job-name=trimm
#SBATCH --error=trimm
#SBATCH --output=trimm
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load bowtie2
module load samtools
module load jdk

ID=$((SLURM_ARRAY_TASK_ID -1))

#samples=(BR10_Demophoon_Brain BR12_Demophoon_Brain BR15_Demophoon_Brain E4-Head BR14_Demophoon_Brain)
#3
#samples=(BR11_Rosina_Brain BR5_Rosina_Brain M4-Head)
#2
#samples=(CH_head2 CH_head6)
#6
samples=(CH_FW2 CH_FW3 CH_FW7 CH_HW2 CH_HW3 CH_HW7)

#REFNAME=Herato
#REFNAME=Hmel2
REFNAME=Hcha1
#REFNAME=Hcha2


samtools index /work/rpapa/share/SVB_Maylin/ATAC/BAM/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.bam