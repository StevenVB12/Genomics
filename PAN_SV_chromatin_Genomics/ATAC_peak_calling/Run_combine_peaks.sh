#!/bin/bash
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=36:00:00
#SBATCH --job-name=macs2
#SBATCH --error=macs2
#SBATCH --output=macs2
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load bedtools

## make sample list command
#echo combining all bed files melpomene
#5
sampleList1=(BR14_Demophoon_Brain BR12_Demophoon_Brain) # BR10_Demophoon_Brain E4-Head BR15_Demophoon_Brain)
#3
sampleList2=(BR11_Rosina_Brain M4-Head) # BR5_Rosina_Brain)
#2
sampleList3=(CH_head2 CH_head6)

ALL_LIST=""
for FILE in ${sampleList1[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done

eval command=\$$(echo ALL_LIST)

## cat bed files
cat $(echo $command) > H_e_dem_peaks.bed

ALL_LIST=""
for FILE in ${sampleList2[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done

eval command=\$$(echo ALL_LIST)

## cat bed files
cat $(echo $command) > H_m_ros_peaks.bed

ALL_LIST=""
for FILE in ${sampleList3[*]}
do
ALL_LIST="$ALL_LIST $FILE"_Hcha1_peaks.narrowPeak""
done

eval command=\$$(echo ALL_LIST)

## cat bed files
cat $(echo $command) > H_cha1_peaks.bed

ALL_LIST=""
for FILE in ${sampleList3[*]}
do
ALL_LIST="$ALL_LIST $FILE"_Hcha2_peaks.narrowPeak""
done

eval command=\$$(echo ALL_LIST)

## cat bed files
cat $(echo $command) > H_cha2_peaks.bed


## sort bed files
#echo sorting bed files erato
sort -k1,1 -k2,2n H_e_dem_peaks.bed > H_e_dem_peaks.sort.bed
#echo sorting bed files melpomene
sort -k1,1 -k2,2n  H_m_ros_peaks.bed >  H_m_ros_peaks.sort.bed

sort -k1,1 -k2,2n  H_cha1_peaks.bed >  H_cha1_peaks.sort.bed
sort -k1,1 -k2,2n  H_cha2_peaks.bed >  H_cha2_peaks.sort.bed


## merge overlapping peaks
#echo merging overlapping peaks erato
bedtools merge -i H_e_dem_peaks.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > H_e_dem_peaks_m1.merged.sort.bed
#echo merging overlapping peaks melpomene
bedtools merge -i H_m_ros_peaks.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > H_m_ros_peaks_m1.merged.sort.bed

bedtools merge -i H_cha1_peaks.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > H_cha1_peaks_m1.merged.sort.bed
bedtools merge -i H_cha2_peaks.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > H_cha2_peaks_m1.merged.sort.bed


bedtools merge -i H_e_dem_peaks.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > H_e_dem_peaks_m0.merged.sort.bed
#echo merging overlapping peaks melpomene
bedtools merge -i H_m_ros_peaks.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > H_m_ros_peaks_m0.merged.sort.bed

bedtools merge -i H_cha1_peaks.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > H_cha1_peaks_m0.merged.sort.bed
bedtools merge -i H_cha2_peaks.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > H_cha2_peaks_m0.merged.sort.bed


# intersect 25%

bedtools intersect -wa -f 0.25 -a BR14_Demophoon_Brain_peaks.narrowPeak -b BR12_Demophoon_Brain_peaks.narrowPeak > H_e_dem_peaks_intersect_25_A.bed
bedtools intersect -wa -f 0.25 -a BR12_Demophoon_Brain_peaks.narrowPeak -b BR14_Demophoon_Brain_peaks.narrowPeak > H_e_dem_peaks_intersect_25_B.bed

bedtools intersect -wa -f 0.25 -a BR11_Rosina_Brain_peaks.narrowPeak -b M4-Head_peaks.narrowPeak > H_m_ros_peaks_intersect_25_A.bed
bedtools intersect -wa -f 0.25 -a M4-Head_peaks.narrowPeak -b BR11_Rosina_Brain_peaks.narrowPeak > H_m_ros_peaks_intersect_25_B.bed

bedtools intersect -wa -f 0.25 -a CH_head2_Hcha1_peaks.narrowPeak -b CH_head6_Hcha1_peaks.narrowPeak > H_cha1_peaks_intersect_25_A.bed
bedtools intersect -wa -f 0.25 -a CH_head6_Hcha1_peaks.narrowPeak -b CH_head2_Hcha1_peaks.narrowPeak > H_cha1_peaks_intersect_25_B.bed

bedtools intersect -wa -f 0.25 -a CH_head2_Hcha2_peaks.narrowPeak -b CH_head6_Hcha2_peaks.narrowPeak > H_cha2_peaks_intersect_25_A.bed
bedtools intersect -wa -f 0.25 -a CH_head6_Hcha2_peaks.narrowPeak -b CH_head2_Hcha2_peaks.narrowPeak > H_cha2_peaks_intersect_25_B.bed

cat H_e_dem_peaks_intersect_25_A.bed H_e_dem_peaks_intersect_25_B.bed | sort -k1,1 -k2,2n | bedtools merge -i - > H_e_dem_peaks_intersect_25.merged.sort.bed
cat H_m_ros_peaks_intersect_25_A.bed H_m_ros_peaks_intersect_25_B.bed | sort -k1,1 -k2,2n | bedtools merge -i - > H_m_ros_peaks_intersect_25.merged.sort.bed
cat H_cha1_peaks_intersect_25_A.bed H_cha1_peaks_intersect_25_B.bed | sort -k1,1 -k2,2n | bedtools merge -i - > H_cha1_peaks_intersect_25.merged.sort.bed
cat H_cha2_peaks_intersect_25_A.bed H_cha2_peaks_intersect_25_B.bed | sort -k1,1 -k2,2n | bedtools merge -i - > H_cha2_peaks_intersect_25.merged.sort.bed
