#!/bin/bash
#SBATCH --mem-per-cpu=24gb
#SBATCH --time=72:00:00
#SBATCH --job-name=bed
#SBATCH --error=bed
#SBATCH --output=bed
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1


module load bedtools
module load samtools

#83
sampleList1=(BR14_Demophoon_Brain BR12_Demophoon_Brain) # BR10_Demophoon_Brain E4-Head BR15_Demophoon_Brain)
#3
sampleList2=(BR11_Rosina_Brain M4-Head) # BR5_Rosina_Brain)
#8
sampleList3=(CH_head2 CH_head6 CH_FW2 CH_FW3 CH_FW7 CH_HW2 CH_HW3 CH_HW7)

## make sample list command
ALL_LIST1=""
for FILE in ${sampleList1[*]}
do
ALL_LIST1="$ALL_LIST1 /work/rpapa/share/SVB_Maylin/ATAC/BAM/${FILE}_Herato.trim.filtered.sorted.nd.bam"
done
echo $ALL_LIST1

## make sample list command
ALL_LIST2=""
for FILE in ${sampleList2[*]}
do
ALL_LIST2="$ALL_LIST2 /work/rpapa/share/SVB_Maylin/ATAC/BAM/${FILE}_Hmel2.trim.filtered.sorted.nd.bam"
done
echo $ALL_LIST2

## make sample list command
ALL_LIST3=""
for FILE in ${sampleList3[*]}
do
ALL_LIST3="$ALL_LIST3 /work/rpapa/share/SVB_Maylin/ATAC/BAM/${FILE}_Hcha1.trim.filtered.sorted.nd.bam"
done
echo $ALL_LIST3

## make sample list command
ALL_LIST4=""
for FILE in ${sampleList3[*]}
do
ALL_LIST4="$ALL_LIST4 /work/rpapa/share/SVB_Maylin/ATAC/BAM/${FILE}_Hcha2.trim.filtered.sorted.nd.bam"
done
echo $ALL_LIST4


eval command1=\$$(echo ALL_LIST1)
eval command2=\$$(echo ALL_LIST2)
eval command3=\$$(echo ALL_LIST3)
eval command4=\$$(echo ALL_LIST4)

# count reads in peaks
#bedtools multicov -bams $(echo $command1) -bed /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/H_e_dem_peaks_intersect_25.merged.sort.bed > /work/rpapa/share/SVB_Maylin/ATAC/counts_sub/H_e_dem_peaks_intersect_25.counts
#bedtools multicov -bams $(echo $command2) -bed /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/H_m_ros_peaks_intersect_25.merged.sort.bed > /work/rpapa/share/SVB_Maylin/ATAC/counts_sub/H_m_ros_peaks_intersect_25.counts
bedtools multicov -bams $(echo $command3) -bed /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/H_cha1_peaks_ALL_m1.merged.sort.bed > /work/rpapa/share/SVB_Maylin/ATAC/counts_sub/H_cha1_peaks_ALL.counts
#bedtools multicov -bams $(echo $command4) -bed /work/rpapa/share/SVB_Maylin/ATAC/MACS2_sub/H_cha2_peaks_intersect_25.merged.sort.bed > /work/rpapa/share/SVB_Maylin/ATAC/counts_sub/H_cha2_peaks_intersect_25.counts

