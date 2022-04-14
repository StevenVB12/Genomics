#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --job-name=trimm
#SBATCH --error=trimm
#SBATCH --output=trimm
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/BR10_Demophoon_Brain_Herato.trim.filtered.sorted.nd.bam -bg -scale 0.4983282 -g /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/BR10_Demophoon_Brain_Herato_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/BR12_Demophoon_Brain_Herato.trim.filtered.sorted.nd.bam -bg -scale 0.6296673 -g /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/BR12_Demophoon_Brain_Herato_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/BR15_Demophoon_Brain_Herato.trim.filtered.sorted.nd.bam -bg -scale 0.8011250 -g /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/BR15_Demophoon_Brain_Herato_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/BR14_Demophoon_Brain_Herato.trim.filtered.sorted.nd.bam -bg -scale 2.7586922 -g /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/BR14_Demophoon_Brain_Herato_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/E4-Head_Herato.trim.filtered.sorted.nd.bam -bg -scale 1.4405454 -g /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/E4-Head_Herato_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/BR11_Rosina_Brain_Hmel2.trim.filtered.sorted.nd.bam -bg -scale 0.7861225 -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/BR11_Rosina_Brain_Hmel2_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/BR5_Rosina_Brain_Hmel2.trim.filtered.sorted.nd.bam -bg -scale 0.5333467 -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/BR5_Rosina_Brain_Hmel2_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/M4-Head_Hmel2.trim.filtered.sorted.nd.bam -bg -scale 2.6462821 -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes > /work/rpapa/share/SVB_Maylin/ATAC/bg/M4-Head_Hmel2_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/CH_head2_Hcha1.trim.filtered.sorted.nd.bam -bg -scale 0.6003002 -g /work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1_bowtie > /work/rpapa/share/SVB_Maylin/ATAC/bg/CH_head2_Hcha1_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/CH_head6_Hcha1.trim.filtered.sorted.nd.bam -bg -scale 1.6658331 -g /work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1_bowtie > /work/rpapa/share/SVB_Maylin/ATAC/bg/CH_head6_Hcha1_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/CH_head2_Hcha2.trim.filtered.sorted.nd.bam -bg -scale 0.6011422 -g /work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1_bowtie > /work/rpapa/share/SVB_Maylin/ATAC/bg/CH_head2_Hcha2_normalised.bg
#/work/rpapa/sbelleghem/Programs/bedtools2/bin/bedtools genomecov -ibam /work/rpapa/share/SVB_Maylin/ATAC/BAM/CH_head6_Hcha2.trim.filtered.sorted.nd.bam -bg -scale 1.6634998 -g /work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1_bowtie > /work/rpapa/share/SVB_Maylin/ATAC/bg/CH_head6_Hcha2_normalised.bg

cd /work/rpapa/share/SVB_Maylin/ATAC/bg

module load wiggletools
wiggletools write_bg brain_H_e_dem_normalized_mean.bg mean BR10_Demophoon_Brain_Herato_normalised.bg BR12_Demophoon_Brain_Herato_normalised.bg BR14_Demophoon_Brain_Herato_normalised.bg BR15_Demophoon_Brain_Herato_normalised.bg E4-Head_Herato_normalised.bg
wiggletools write_bg brain_H_e_ros_normalized_mean.bg mean BR11_Rosina_Brain_Hmel2_normalised.bg BR5_Rosina_Brain_Hmel2_normalised.bg M4-Head_Hmel2_normalised.bg
wiggletools write_bg brain_H_cha1_normalized_mean.bg mean CH_head2_Hcha1_normalised.bg CH_head6_Hcha1_normalised.bg
wiggletools write_bg brain_H_cha2_normalized_mean.bg mean CH_head2_Hcha2_normalised.bg CH_head6_Hcha2_normalised.bg