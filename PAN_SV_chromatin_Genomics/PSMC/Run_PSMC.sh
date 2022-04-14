#!/bin/bash
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=48:00:00
#SBATCH --job-name=psmc
#SBATCH --error=psmc
#SBATCH --output=psmc
#SBATCH --ntasks=16
#SBATCH --partition=rpapa


module load samtools
module load python3
module load bcftools
module load bwa
module load jdk

indID=$((SLURM_ARRAY_TASK_ID -1))

REF=/work/rpapa/share/REF/H_charithonia_PR/Hcha.assembly.v1.1.fasta

#zcat /work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3851_HTTKNDSXX_L1_1.fq.gz \
#/work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3852_HTTKNDSXX_L1_1.fq.gz \
#/work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3853_HTTKNDSXX_L1_1.fq.gz \
#/work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3854_HTTKNDSXX_L1_1.fq.gz | gzip > /work/rpapa/share/SVB_Maylin/PSMC/SVB003.R1.fastq.gz

#zcat /work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3851_HTTKNDSXX_L1_2.fq.gz \
#/work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3852_HTTKNDSXX_L1_2.fq.gz \
#/work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3853_HTTKNDSXX_L1_2.fq.gz \
#/work/rpapa/share/10X_SVB/H_charithonia_10X_1/H202SC19100641/Rawdata/SVB_003/SVB_003_FDHX192247700-1a-AK3854_HTTKNDSXX_L1_2.fq.gz | gzip > /work/rpapa/share/SVB_Maylin/PSMC/SVB003.R2.fastq.gz

#bwa index $REF

#bwa mem -t 16 -M $REF /work/rpapa/share/SVB_Maylin/PSMC/SVB003.R1.fastq.gz /work/rpapa/share/SVB_Maylin/PSMC/SVB003.R2.fastq.gz | samtools view -bS - > /work/rpapa/share/SVB_Maylin/PSMC/SVB003.bam
#bwa mem -t 16 -M $REF /work/rpapa/share/reads_erato/charithonia/char_001.read1.gz /work/rpapa/share/reads_erato/charithonia/char_001.read2.gz | samtools view -bS - > /work/rpapa/share/SVB_Maylin/PSMC/char_001.bam
bwa mem -t 16 -M $REF /work/rpapa/share/reads_erato/charithonia/char_002.read1.gz /work/rpapa/share/reads_erato/charithonia/char_002.read2.gz | samtools view -bS - > /work/rpapa/share/SVB_Maylin/PSMC/char_002.bam

#samtools view -f 0x02 -q 20 -b /work/rpapa/share/SVB_Maylin/PSMC/SVB003.bam > /work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.bam
#samtools view -f 0x02 -q 20 -b /work/rpapa/share/SVB_Maylin/PSMC/char_001.bam > /work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.bam
samtools view -f 0x02 -q 20 -b /work/rpapa/share/SVB_Maylin/PSMC/char_002.bam > /work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.bam

#samtools sort /work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.bam /work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.sorted
#samtools sort /work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.bam /work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.sorted
samtools sort /work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.bam /work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.sorted

#java -jar /work/rpapa/sbelleghem/Programs/picard-tools-2.5.0/picard.jar MarkDuplicates \
#I=/work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.sorted.bam \
#O=/work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.sorted.nd.bam \
#Remove_Duplicates=true  M=/work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.sorted_dup_metrics.txt ASSUME_SORTED=true

#java -jar /work/rpapa/sbelleghem/Programs/picard-tools-2.5.0/picard.jar MarkDuplicates \
#I=/work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.sorted.bam \
#O=/work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.sorted.nd.bam \
#Remove_Duplicates=true  M=/work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.sorted_dup_metrics.txt ASSUME_SORTED=true

java -jar /work/rpapa/sbelleghem/Programs/picard-tools-2.5.0/picard.jar MarkDuplicates \
I=/work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.sorted.bam \
O=/work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.sorted.nd.bam \
Remove_Duplicates=true  M=/work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.sorted_dup_metrics.txt ASSUME_SORTED=true


#samtools index /work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.sorted.nd.bam
#samtools index /work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.sorted.nd.bam
samtools index /work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.sorted.nd.bam

COMMAND=""

#for SCAF in HchaChr201001 HchaChr201002 HchaChr202001 HchaChr203003 HchaChr204001 HchaChr204002 HchaChr205001 HchaChr205002 HchaChr206001 HchaChr206002 HchaChr207001 HchaChr208001 HchaChr208002 HchaChr208003 HchaChr208004 HchaChr209001 HchaChr210001 HchaChr210002 HchaChr210003 HchaChr211001 HchaChr211002 HchaChr212001 HchaChr212002 HchaChr213001 HchaChr213002 HchaChr2140041 HchaChr2140042 HchaChr215003 HchaChr2160021 HchaChr2160022 HchaChr217001 HchaChr217002 HchaChr217003 HchaChr218003 HchaChr219001 HchaChr219002 HchaChr2200031 HchaChr2200032
#do

#~/work/Programs/samtools-0.1.19/samtools mpileup -q 20 -Q 20 -C 50 -u -r $SCAF -f $REF /work/rpapa/share/SVB_Maylin/PSMC/SVB003.filt.sorted.nd.bam |\
# ~/work/Programs/samtools-0.1.19/bcftools/bcftools view -cgI - |\
# ~/work/Programs/msmc-tools-master/bamCaller.py 30 /work/rpapa/share/SVB_Maylin/PSMC/out/SVB003_$SCAF.mask.bed.gz |\
# gzip -c > /work/rpapa/share/SVB_Maylin/PSMC/out/SVB003_$SCAF.vcf.gz

#~/work/Programs/msmc-tools-master/generate_multihetsep.py --mask=/work/rpapa/share/SVB_Maylin/PSMC/out/SVB003_$SCAF.mask.bed.gz /work/rpapa/share/SVB_Maylin/PSMC/out/SVB003_$SCAF.vcf.gz >\
# /work/rpapa/share/SVB_Maylin/PSMC/out/SVB003_$SCAF.txt

#done

#for SCAF in HchaChr201001 HchaChr201002 HchaChr202001 HchaChr203003 HchaChr204001 HchaChr204002 HchaChr205001 HchaChr205002 HchaChr206001 HchaChr206002 HchaChr207001 HchaChr208001 HchaChr208002 HchaChr208003 HchaChr208004 HchaChr209001 HchaChr210001 HchaChr210002 HchaChr210003 HchaChr211001 HchaChr211002 HchaChr212001 HchaChr212002 HchaChr213001 HchaChr213002 HchaChr2140041 HchaChr2140042 HchaChr215003 HchaChr2160021 HchaChr2160022 HchaChr217001 HchaChr217002 HchaChr217003 HchaChr218003 HchaChr219001 HchaChr219002 HchaChr2200031 HchaChr2200032
#do

#~/work/Programs/samtools-0.1.19/samtools mpileup -q 20 -Q 20 -C 50 -u -r $SCAF -f $REF /work/rpapa/share/SVB_Maylin/PSMC/char_001.filt.sorted.nd.bam |\
# ~/work/Programs/samtools-0.1.19/bcftools/bcftools view -cgI - |\
# ~/work/Programs/msmc-tools-master/bamCaller.py 30 /work/rpapa/share/SVB_Maylin/PSMC/out/char_001_$SCAF.mask.bed.gz |\
# gzip -c > /work/rpapa/share/SVB_Maylin/PSMC/out/char_001_$SCAF.vcf.gz

#~/work/Programs/msmc-tools-master/generate_multihetsep.py --mask=/work/rpapa/share/SVB_Maylin/PSMC/out/char_001_$SCAF.mask.bed.gz /work/rpapa/share/SVB_Maylin/PSMC/out/char_001_$SCAF.vcf.gz >\
# /work/rpapa/share/SVB_Maylin/PSMC/out/char_001_$SCAF.txt

#done

for SCAF in HchaChr201001 HchaChr201002 HchaChr202001 HchaChr203003 HchaChr204001 HchaChr204002 HchaChr205001 HchaChr205002 HchaChr206001 HchaChr206002 HchaChr207001 HchaChr208001 HchaChr208002 HchaChr208003 HchaChr208004 HchaChr209001 HchaChr210001 HchaChr210002 HchaChr210003 HchaChr211001 HchaChr211002 HchaChr212001 HchaChr212002 HchaChr213001 HchaChr213002 HchaChr2140041 HchaChr2140042 HchaChr215003 HchaChr2160021 HchaChr2160022 HchaChr217001 HchaChr217002 HchaChr217003 HchaChr218003 HchaChr219001 HchaChr219002 HchaChr2200031 HchaChr2200032
do

~/work/Programs/samtools-0.1.19/samtools mpileup -q 20 -Q 20 -C 50 -u -r $SCAF -f $REF /work/rpapa/share/SVB_Maylin/PSMC/char_002.filt.sorted.nd.bam |\
 ~/work/Programs/samtools-0.1.19/bcftools/bcftools view -cgI - |\
 ~/work/Programs/msmc-tools-master/bamCaller.py 30 /work/rpapa/share/SVB_Maylin/PSMC/out/char_002_$SCAF.mask.bed.gz |\
 gzip -c > /work/rpapa/share/SVB_Maylin/PSMC/out/char_002_$SCAF.vcf.gz

~/work/Programs/msmc-tools-master/generate_multihetsep.py --mask=/work/rpapa/share/SVB_Maylin/PSMC/out/char_002_$SCAF.mask.bed.gz /work/rpapa/share/SVB_Maylin/PSMC/out/char_002_$SCAF.vcf.gz >\
 /work/rpapa/share/SVB_Maylin/PSMC/out/char_002_$SCAF.txt

done


for SCAF in HchaChr201001 HchaChr201002 HchaChr202001 HchaChr203003 HchaChr204001 HchaChr204002 HchaChr205001 HchaChr205002 HchaChr206001 HchaChr206002 HchaChr207001 HchaChr208001 HchaChr208002 HchaChr208003 HchaChr208004 HchaChr209001 HchaChr210001 HchaChr210002 HchaChr210003 HchaChr211001 HchaChr211002 HchaChr212001 HchaChr212002 HchaChr213001 HchaChr213002 HchaChr2140041 HchaChr2140042 HchaChr215003 HchaChr2160021 HchaChr2160022 HchaChr217001 HchaChr217002 HchaChr217003 HchaChr218003 HchaChr219001 HchaChr219002 HchaChr2200031 HchaChr2200032
do

COMMAND="$COMMAND char_002_$SCAF.txt"

done

/rds/user/sv378/hpc-work/Programs/msmc_1.0.0_linux64bit -t 8 -o char_002 $COMMAND


#for FINAL in SVB003 
#do
#python ~/work/scripts/MSMC_plotInput.py -I /work/rpapa/share/SVB_Maylin/PSMC/out2/$FINAL.final.txt -u 2e-09 -g 0.25 > /work/rpapa/share/SVB_Maylin/PSMC/out2/$FINAL.final.Rin.txt
#done
