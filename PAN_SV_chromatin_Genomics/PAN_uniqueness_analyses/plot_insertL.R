#calculate chromosome coordinates
chrom.coords <- function(scafL,chromNames,gap = 2000000) {
  chromosome = vector()
  chromLengths  = vector()
  chromStarts = vector()
  chromEnds = vector()
  chromMid = vector()
  chrom = 1
  endLast = 0
  scafCurrent <- subset(scafL, chromosome == chromNames[1])
  chromosome[chrom] <- chrom
  chromLengths[chrom] <- sum(scafCurrent$length)
  chromStarts[chrom] <- endLast + 1
  chromEnds[chrom] <- endLast + chromLengths[chrom]
  chromMid[chrom] <- endLast + chromLengths[chrom]/2
  endLast = chromEnds[chrom]
  chrom = chrom + 1
  for (i in 2:length(chromNames)) {
    chromosome[chrom] <- chrom
    scafCurrent <- subset(scafL, chromosome == chromNames[i])
    chromLengths[chrom] <- sum(scafCurrent$length)
    chromStarts[chrom] <- endLast + gap + 1
    chromEnds[chrom] <- endLast + gap + chromLengths[chrom]
    chromMid[chrom] <- endLast + gap + chromLengths[chrom]/2
    endLast = chromEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(chromosome,chromLengths,chromStarts,chromEnds,chromMid))
  return(table)
}

#calculate scaffold coordinates
scaf.coords <- function(scafL,gap = 0) {
  scaffold = vector()
  scafStarts = vector()
  scafEnds = vector()
  chrom = 1
  endLast = 0
  for (e in 1:nrow(scafL)) {
    scaffold[chrom] <- levels(scafL$scaffold)[e]
    scafStarts[chrom] <- endLast + gap + 1
    scafEnds[chrom] <- endLast + gap + scafL$length[e]
    endLast = scafEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(scaffold,scafStarts,scafEnds))
  return(table)
}


###
scafL_dem <- read.table("inserts/scaffold_lengths_chrom_erato.txt", h=T)
chromNames <-c(1:21)
chrom_coords_dem <- chrom.coords(scafL_dem, chromNames)
scaf_coords_dem <- scaf.coords(scafL_dem)

scaf_coords2_dem <- merge(scafL_dem, chrom_coords_dem, by="chromosome", all.x=TRUE)

scafL_melp <- read.table("inserts/scaffold_lengths_chrom_melp.txt", h=T)
chromNames <-c(1:21)
chrom_coords_melp <- chrom.coords(scafL_melp, chromNames)
scaf_coords_melp <- scaf.coords(scafL_melp)

scaf_coords2_melp <- merge(scafL_melp, chrom_coords_melp, by="chromosome", all.x=TRUE)

scafL_cha <- read.table("inserts/scaffold_lengths_chrom_char.txt", h=T)
chromNames <-c(1:21)
chrom_coords_cha <- chrom.coords(scafL_cha, chromNames)
scaf_coords_cha <- scaf.coords(scafL_cha)

scaf_coords2_cha <- merge(scafL_cha, chrom_coords_melp, by="chromosome", all.x=TRUE)
######

chrom <- read.table("seqseqpan/chromosome_positions_pan.txt", h=T)
head(chrom)

chrom$length <- chrom$end -chrom$start

inserts_dem_blocks <- read.table("inserts/blocks_unique_1.txt")
inserts_melp_blocks <- read.table("inserts/blocks_unique_4.txt")
inserts_dem_missing <- read.table("inserts/missing_unique_1.txt")
inserts_melp_missing <- read.table("inserts/missing_unique_4.txt")

inserts_cha1_blocks <- read.table("inserts/blocks_unique_2.txt")
inserts_cha2_blocks <- read.table("inserts/blocks_unique_3.txt")

colnames(inserts_dem_blocks) <- c('scaf','start','end')
colnames(inserts_melp_blocks) <- c('scaf','start','end')
colnames(inserts_dem_missing) <- c('scaf','start','end')
colnames(inserts_melp_missing) <- c('scaf','start','end')

inserts_dem_blocks$insertL <- inserts_dem_blocks$end-inserts_dem_blocks$start
inserts_melp_blocks$insertL <- inserts_melp_blocks$end-inserts_melp_blocks$start
inserts_dem_missing$insertL <- inserts_dem_missing$end-inserts_dem_missing$start
inserts_melp_missing$insertL <- inserts_melp_missing$end-inserts_melp_missing$start


inserts_melpvsdem_blocks <- read.table("inserts/blocks_unique_4_vs_1.txt")
inserts_demvsmelp_blocks <- read.table("inserts/blocks_unique_1_vs_4.txt")

inserts_chavsera_blocks <- read.table("inserts/blocks_unique_2_vs_1.txt")
inserts_chavsmelp_blocks <- read.table("inserts/blocks_unique_2_vs_4.txt")

inserts_cha2vsera_blocks <- read.table("inserts/blocks_unique_3_vs_1.txt")
inserts_cha2vsmelp_blocks <- read.table("inserts/blocks_unique_3_vs_4.txt")

inserts_cha1vscha2_blocks <- read.table("inserts/blocks_unique_2_vs_3.txt")
inserts_cha2vscha1_blocks <- read.table("inserts/blocks_unique_3_vs_2.txt")

sum(inserts_cha1vscha2_blocks$V3-inserts_cha1vscha2_blocks$V2)
sum(inserts_cha2vscha1_blocks$V3-inserts_cha2vscha1_blocks$V2)

inserts_chavscha_blocks <- rbind(inserts_cha1vscha2_blocks,inserts_cha2vscha1_blocks)

inserts_melpvscha_blocks <- read.table("inserts/blocks_unique_4_vs_2_3.txt")
inserts_demvscha_blocks <- read.table("inserts/blocks_unique_1_vs_2_3.txt")
inserts_melpvschadem_blocks <- read.table("inserts/blocks_unique_4_vs_1_2_3.txt")

inserts_cha23vseramelp_blocks <- read.table("inserts/blocks_unique_2_3_combined.txt")
inserts_cha23vsmelp_blocks <- read.table("inserts/blocks_unique_2_3_vs_4.txt")
inserts_cha23vsera_blocks <- read.table("inserts/blocks_unique_2_3_vs_1.txt")

sum(inserts_cha23vseramelp_blocks$V3-inserts_cha23vseramelp_blocks$V2)
sum(inserts_cha23vsmelp_blocks$V3-inserts_cha23vsmelp_blocks$V2)
sum(inserts_cha23vsera_blocks$V3-inserts_cha23vsera_blocks$V2)

blocks1 <- read.table("inserts/1_blocks_intervals.corr.bed")
blocks2 <- read.table("inserts/2_blocks_intervals.corr.bed")
blocks3 <- read.table("inserts/3_blocks_intervals.corr.bed")
blocks4 <- read.table("inserts/4_blocks_intervals.corr.bed")

indels23_2 <- read.table("inserts/indels_2_3_uniqueEra_2.txt")
indels23_3 <- read.table("inserts/indels_2_3_uniqueEra_3.txt")

indels23 <- read.table("inserts/indels_2_3_uniqueEra.txt")

sum(indels23_2$V3-indels23_2$V2)
sum(indels23_3$V3-indels23_3$V2)

sum(indels23$V3-indels23$V2)

blocks23 <- read.table("inserts/blocks_merged_2_3.txt")


sum(blocks1$V3-blocks1$V2)
sum(blocks2$V3-blocks2$V2)
sum(blocks3$V3-blocks3$V2)
sum(blocks4$V3-blocks4$V2)

sum(blocks23$V3-blocks23$V2)

inserts_chavscha_blocks_s     <- subset(inserts_chavscha_blocks, inserts_chavscha_blocks$V2 <= 578665626)
inserts_chavsera_blocks_s     <- subset(inserts_cha23vsera_blocks, inserts_cha23vsera_blocks$V2 <= 578665626)
inserts_chavsmelp_blocks_s    <- subset(inserts_cha23vsmelp_blocks, inserts_cha23vsmelp_blocks$V2 <= 578665626)
inserts_melpvschadem_blocks_s <- subset(inserts_melpvschadem_blocks, inserts_melpvschadem_blocks$V2 <= 578665626)

inserts_dem_blocks_s <- subset(inserts_dem_blocks, inserts_dem_blocks$V2 <= 578665626)
inserts_cha23vseramelp_blocks_s <- subset(inserts_cha23vseramelp_blocks, inserts_cha23vseramelp_blocks$V2 <= 578665626)

par(mfrow=c(1,4))
hist(log10((as.numeric(as.character(inserts_chavscha_blocks_s$V3))-as.numeric(as.character(inserts_chavscha_blocks_s$V2)))), breaks=32, main='', xlab='', ylab='',         xlim=c(0,3),   ylim=c(0,800000))
hist(log10((as.numeric(as.character(inserts_chavsera_blocks_s$V3))-as.numeric(as.character(inserts_chavsera_blocks_s$V2)))), breaks=32, main='', xlab='', ylab='',         xlim=c(0,3),   ylim=c(0,800000))
hist(log10((as.numeric(as.character(inserts_chavsmelp_blocks_s$V3))-as.numeric(as.character(inserts_chavsmelp_blocks_s$V2)))), breaks=32, main='', xlab='', ylab='',       xlim=c(0,3), ylim=c(0,800000))
# hist(log10((as.numeric(as.character(inserts_demvsmelp_blocks$V3))-as.numeric(as.character(inserts_demvsmelp_blocks$V2)))), breaks=32, main='', xlab='', ylab='',       xlim=c(0,3), ylim=c(0,800000))
hist(log10((as.numeric(as.character(inserts_melpvschadem_blocks_s$V3))-as.numeric(as.character(inserts_melpvschadem_blocks_s$V2)))), breaks=32, main='', xlab='', ylab='', xlim=c(0,3), ylim=c(0,800000))



###
sum(inserts_chavscha_blocks$V3-inserts_chavscha_blocks$V2)
sum(inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2)
mean(inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2)
median(inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2)

inserts_chavscha_blocks_s_1 <- subset(inserts_chavscha_blocks_s, inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2 ==1)
inserts_chavscha_blocks_s_20 <- subset(inserts_chavscha_blocks_s, inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2 <20)
inserts_chavscha_blocks_s_1000 <- subset(inserts_chavscha_blocks_s, inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2 >1000)

sum1 <- sum(inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2)
nrow1 <- nrow(inserts_chavscha_blocks_s)

sum(inserts_chavscha_blocks_s_1$V3-inserts_chavscha_blocks_s_1$V2)  
sum(inserts_chavscha_blocks_s_1$V3-inserts_chavscha_blocks_s_1$V2)/sum1
sum(inserts_chavscha_blocks_s_1$V3-inserts_chavscha_blocks_s_1$V2)/403200000
nrow(inserts_chavscha_blocks_s_1)/nrow1

sum(inserts_chavscha_blocks_s_20$V3-inserts_chavscha_blocks_s_20$V2)
sum(inserts_chavscha_blocks_s_20$V3-inserts_chavscha_blocks_s_20$V2)/sum1
sum(inserts_chavscha_blocks_s_20$V3-inserts_chavscha_blocks_s_20$V2)/403200000
nrow(inserts_chavscha_blocks_s_20)/nrow1

sum(inserts_chavscha_blocks_s_1000$V3-inserts_chavscha_blocks_s_1000$V2)
sum(inserts_chavscha_blocks_s_1000$V3-inserts_chavscha_blocks_s_1000$V2)/sum1
sum(inserts_chavscha_blocks_s_1000$V3-inserts_chavscha_blocks_s_1000$V2)/403200000
nrow(inserts_chavscha_blocks_s_1000)/nrow1

max((inserts_chavscha_blocks_s$V3-inserts_chavscha_blocks_s$V2))

###
inserts_cha23vsera_blocks_s <- subset(inserts_cha23vsera_blocks, inserts_cha23vsera_blocks$V2 <= 578665626)

sum(inserts_cha23vsera_blocks$V3-inserts_cha23vsera_blocks$V2)
sum(inserts_cha23vsera_blocks_s$V3-inserts_cha23vsera_blocks_s$V2)

inserts_cha23vsera_blocks_s_1 <- subset(inserts_cha23vsera_blocks_s, inserts_cha23vsera_blocks_s$V3-inserts_cha23vsera_blocks_s$V2 ==1)
inserts_cha23vsera_blocks_s_20 <- subset(inserts_cha23vsera_blocks_s, inserts_cha23vsera_blocks_s$V3-inserts_cha23vsera_blocks_s$V2 <20)
inserts_cha23vsera_blocks_s_1000 <- subset(inserts_cha23vsera_blocks_s, inserts_cha23vsera_blocks_s$V3-inserts_cha23vsera_blocks_s$V2 >1000)

sum1 <- sum(inserts_cha23vsera_blocks_s$V3-inserts_cha23vsera_blocks_s$V2)
nrow1 <- nrow(inserts_cha23vsera_blocks_s)

sum(inserts_cha23vsera_blocks_s_1$V3-inserts_cha23vsera_blocks_s_1$V2)  
sum(inserts_cha23vsera_blocks_s_1$V3-inserts_cha23vsera_blocks_s_1$V2)/sum1
sum(inserts_cha23vsera_blocks_s_1$V3-inserts_cha23vsera_blocks_s_1$V2)/403200000
nrow(inserts_cha23vsera_blocks_s_1)/nrow1

sum(inserts_cha23vsera_blocks_s_20$V3-inserts_cha23vsera_blocks_s_20$V2)
sum(inserts_cha23vsera_blocks_s_20$V3-inserts_cha23vsera_blocks_s_20$V2)/sum1
sum(inserts_cha23vsera_blocks_s_20$V3-inserts_cha23vsera_blocks_s_20$V2)/403200000
nrow(inserts_cha23vsera_blocks_s_20)/nrow1

sum(inserts_cha23vsera_blocks_s_1000$V3-inserts_cha23vsera_blocks_s_1000$V2)
sum(inserts_cha23vsera_blocks_s_1000$V3-inserts_cha23vsera_blocks_s_1000$V2)/sum1
sum(inserts_cha23vsera_blocks_s_1000$V3-inserts_cha23vsera_blocks_s_1000$V2)/403200000
nrow(inserts_cha23vsera_blocks_s_1000)/nrow1


###
inserts_demvscha_blocks_s <- subset(inserts_demvscha_blocks, inserts_demvscha_blocks$V2 <= 578665626)

sum(inserts_demvscha_blocks$V3-inserts_demvscha_blocks$V2)
sum(inserts_demvscha_blocks_s$V3-inserts_demvscha_blocks_s$V2)

inserts_demvscha_blocks_s_1 <- subset(inserts_demvscha_blocks_s, inserts_demvscha_blocks_s$V3-inserts_demvscha_blocks_s$V2 ==1)
inserts_demvscha_blocks_s_20 <- subset(inserts_demvscha_blocks_s, inserts_demvscha_blocks_s$V3-inserts_demvscha_blocks_s$V2 <20)
inserts_demvscha_blocks_s_1000 <- subset(inserts_demvscha_blocks_s, inserts_demvscha_blocks_s$V3-inserts_demvscha_blocks_s$V2 >1000)

sum1 <- sum(inserts_demvscha_blocks_s$V3-inserts_demvscha_blocks_s$V2)
nrow1 <- nrow(inserts_demvscha_blocks_s)

sum(inserts_demvscha_blocks_s_1$V3-inserts_demvscha_blocks_s_1$V2)  
sum(inserts_demvscha_blocks_s_1$V3-inserts_demvscha_blocks_s_1$V2)/sum1
sum(inserts_demvscha_blocks_s_1$V3-inserts_demvscha_blocks_s_1$V2)/382800000
nrow(inserts_demvscha_blocks_s_1)/nrow1

sum(inserts_demvscha_blocks_s_20$V3-inserts_demvscha_blocks_s_20$V2)
sum(inserts_demvscha_blocks_s_20$V3-inserts_demvscha_blocks_s_20$V2)/sum1
sum(inserts_demvscha_blocks_s_20$V3-inserts_demvscha_blocks_s_20$V2)/382800000
nrow(inserts_demvscha_blocks_s_20)/nrow1

sum(inserts_demvscha_blocks_s_1000$V3-inserts_demvscha_blocks_s_1000$V2)
sum(inserts_demvscha_blocks_s_1000$V3-inserts_demvscha_blocks_s_1000$V2)/sum1
sum(inserts_demvscha_blocks_s_1000$V3-inserts_demvscha_blocks_s_1000$V2)/382800000
nrow(inserts_demvscha_blocks_s_1000)/nrow1

##
inserts_cha23vsmelp_blocks_s <- subset(inserts_cha23vsmelp_blocks, inserts_cha23vsmelp_blocks$V2 <= 578665626)

sum(inserts_cha23vsmelp_blocks$V3-inserts_cha23vsmelp_blocks$V2)
sum(inserts_cha23vsmelp_blocks_s$V3-inserts_cha23vsmelp_blocks_s$V2)

inserts_cha23vsmelp_blocks_s_1 <- subset(inserts_cha23vsmelp_blocks_s, inserts_cha23vsmelp_blocks_s$V3-inserts_cha23vsmelp_blocks_s$V2 ==1)
inserts_cha23vsmelp_blocks_s_20 <- subset(inserts_cha23vsmelp_blocks_s, inserts_cha23vsmelp_blocks_s$V3-inserts_cha23vsmelp_blocks_s$V2 <20)
inserts_cha23vsmelp_blocks_s_1000 <- subset(inserts_cha23vsmelp_blocks_s, inserts_cha23vsmelp_blocks_s$V3-inserts_cha23vsmelp_blocks_s$V2 >1000)

sum1 <- sum(inserts_cha23vsmelp_blocks_s$V3-inserts_cha23vsmelp_blocks_s$V2)
nrow1 <- nrow(inserts_cha23vsmelp_blocks_s)

sum(inserts_cha23vsmelp_blocks_s_1$V3-inserts_cha23vsmelp_blocks_s_1$V2)  
sum(inserts_cha23vsmelp_blocks_s_1$V3-inserts_cha23vsmelp_blocks_s_1$V2)/sum1
sum(inserts_cha23vsmelp_blocks_s_1$V3-inserts_cha23vsmelp_blocks_s_1$V2)/403200000
nrow(inserts_cha23vsmelp_blocks_s_1)/nrow1

sum(inserts_cha23vsmelp_blocks_s_20$V3-inserts_cha23vsmelp_blocks_s_20$V2)
sum(inserts_cha23vsmelp_blocks_s_20$V3-inserts_cha23vsmelp_blocks_s_20$V2)/sum1
sum(inserts_cha23vsmelp_blocks_s_20$V3-inserts_cha23vsmelp_blocks_s_20$V2)/403200000
nrow(inserts_cha23vsmelp_blocks_s_20)/nrow1

sum(inserts_cha23vsmelp_blocks_s_1000$V3-inserts_cha23vsmelp_blocks_s_1000$V2)
sum(inserts_cha23vsmelp_blocks_s_1000$V3-inserts_cha23vsmelp_blocks_s_1000$V2)/sum1
sum(inserts_cha23vsmelp_blocks_s_1000$V3-inserts_cha23vsmelp_blocks_s_1000$V2)/403200000
nrow(inserts_cha23vsmelp_blocks_s_1000)/nrow1


##
inserts_demvsmelp_blocks_s <- subset(inserts_demvsmelp_blocks, inserts_demvsmelp_blocks$V2 <= 578665626)

sum(inserts_demvsmelp_blocks$V3-inserts_demvsmelp_blocks$V2)
sum(inserts_demvsmelp_blocks_s$V3-inserts_demvsmelp_blocks_s$V2)

inserts_demvsmelp_blocks_s_1 <- subset(inserts_demvsmelp_blocks_s, inserts_demvsmelp_blocks_s$V3-inserts_demvsmelp_blocks_s$V2 ==1)
inserts_demvsmelp_blocks_s_20 <- subset(inserts_demvsmelp_blocks_s, inserts_demvsmelp_blocks_s$V3-inserts_demvsmelp_blocks_s$V2 <20)
inserts_demvsmelp_blocks_s_1000 <- subset(inserts_demvsmelp_blocks_s, inserts_demvsmelp_blocks_s$V3-inserts_demvsmelp_blocks_s$V2 >1000)

sum1 <- sum(inserts_demvsmelp_blocks_s$V3-inserts_demvsmelp_blocks_s$V2)
nrow1 <- nrow(inserts_demvsmelp_blocks_s)

sum(inserts_demvsmelp_blocks_s_1$V3-inserts_demvsmelp_blocks_s_1$V2)  
sum(inserts_demvsmelp_blocks_s_1$V3-inserts_demvsmelp_blocks_s_1$V2)/sum1
sum(inserts_demvsmelp_blocks_s_1$V3-inserts_demvsmelp_blocks_s_1$V2)/382800000
nrow(inserts_demvsmelp_blocks_s_1)/nrow1

sum(inserts_demvsmelp_blocks_s_20$V3-inserts_demvsmelp_blocks_s_20$V2)
sum(inserts_demvsmelp_blocks_s_20$V3-inserts_demvsmelp_blocks_s_20$V2)/sum1
sum(inserts_demvsmelp_blocks_s_20$V3-inserts_demvsmelp_blocks_s_20$V2)/382800000
nrow(inserts_demvsmelp_blocks_s_20)/nrow1

sum(inserts_demvsmelp_blocks_s_1000$V3-inserts_demvsmelp_blocks_s_1000$V2)
sum(inserts_demvsmelp_blocks_s_1000$V3-inserts_demvsmelp_blocks_s_1000$V2)/sum1
sum(inserts_demvsmelp_blocks_s_1000$V3-inserts_demvsmelp_blocks_s_1000$V2)/382800000
nrow(inserts_demvsmelp_blocks_s_1000)/nrow1

##
inserts_melpvscha_blocks_s <- subset(inserts_melpvscha_blocks, inserts_melpvscha_blocks$V2 <= 578665626)

sum(inserts_melpvscha_blocks$V3-inserts_melpvscha_blocks$V2)
sum(inserts_melpvscha_blocks_s$V3-inserts_melpvscha_blocks_s$V2)

inserts_melpvscha_blocks_s_1 <- subset(inserts_melpvscha_blocks_s, inserts_melpvscha_blocks_s$V3-inserts_melpvscha_blocks_s$V2 ==1)
inserts_melpvscha_blocks_s_20 <- subset(inserts_melpvscha_blocks_s, inserts_melpvscha_blocks_s$V3-inserts_melpvscha_blocks_s$V2 <20)
inserts_melpvscha_blocks_s_1000 <- subset(inserts_melpvscha_blocks_s, inserts_melpvscha_blocks_s$V3-inserts_melpvscha_blocks_s$V2 >1000)

sum1 <- sum(inserts_melpvscha_blocks_s$V3-inserts_melpvscha_blocks_s$V2)
nrow1 <- nrow(inserts_melpvscha_blocks_s)

sum(inserts_melpvscha_blocks_s_1$V3-inserts_melpvscha_blocks_s_1$V2)  
sum(inserts_melpvscha_blocks_s_1$V3-inserts_melpvscha_blocks_s_1$V2)/sum1
sum(inserts_melpvscha_blocks_s_1$V3-inserts_melpvscha_blocks_s_1$V2)/275200000
nrow(inserts_melpvscha_blocks_s_1)/nrow1

sum(inserts_melpvscha_blocks_s_20$V3-inserts_melpvscha_blocks_s_20$V2)
sum(inserts_melpvscha_blocks_s_20$V3-inserts_melpvscha_blocks_s_20$V2)/sum1
sum(inserts_melpvscha_blocks_s_20$V3-inserts_melpvscha_blocks_s_20$V2)/275200000
nrow(inserts_melpvscha_blocks_s_20)/nrow1

sum(inserts_melpvscha_blocks_s_1000$V3-inserts_melpvscha_blocks_s_1000$V2)
sum(inserts_melpvscha_blocks_s_1000$V3-inserts_melpvscha_blocks_s_1000$V2)/sum1
sum(inserts_melpvscha_blocks_s_1000$V3-inserts_melpvscha_blocks_s_1000$V2)/275200000
nrow(inserts_melpvscha_blocks_s_1000)/nrow1

###
inserts_melpvsdem_blocks_s <- subset(inserts_melpvsdem_blocks, inserts_melpvsdem_blocks$V2 <= 578665626)

sum(inserts_melpvsdem_blocks$V3-inserts_melpvsdem_blocks$V2)
sum(inserts_melpvsdem_blocks_s$V3-inserts_melpvsdem_blocks_s$V2)

inserts_melpvsdem_blocks_s_1 <- subset(inserts_melpvsdem_blocks_s, inserts_melpvsdem_blocks_s$V3-inserts_melpvsdem_blocks_s$V2 ==1)
inserts_melpvsdem_blocks_s_20 <- subset(inserts_melpvsdem_blocks_s, inserts_melpvsdem_blocks_s$V3-inserts_melpvsdem_blocks_s$V2 <20)
inserts_melpvsdem_blocks_s_1000 <- subset(inserts_melpvsdem_blocks_s, inserts_melpvsdem_blocks_s$V3-inserts_melpvsdem_blocks_s$V2 >1000)

sum1 <- sum(inserts_melpvsdem_blocks_s$V3-inserts_melpvsdem_blocks_s$V2)
nrow1 <- nrow(inserts_melpvsdem_blocks_s)

sum(inserts_melpvsdem_blocks_s_1$V3-inserts_melpvsdem_blocks_s_1$V2)  
sum(inserts_melpvsdem_blocks_s_1$V3-inserts_melpvsdem_blocks_s_1$V2)/sum1
sum(inserts_melpvsdem_blocks_s_1$V3-inserts_melpvsdem_blocks_s_1$V2)/275200000
nrow(inserts_melpvsdem_blocks_s_1)/nrow1

sum(inserts_melpvsdem_blocks_s_20$V3-inserts_melpvsdem_blocks_s_20$V2)
sum(inserts_melpvsdem_blocks_s_20$V3-inserts_melpvsdem_blocks_s_20$V2)/sum1
sum(inserts_melpvsdem_blocks_s_20$V3-inserts_melpvsdem_blocks_s_20$V2)/275200000
nrow(inserts_melpvsdem_blocks_s_20)/nrow1

sum(inserts_melpvsdem_blocks_s_1000$V3-inserts_melpvsdem_blocks_s_1000$V2)
sum(inserts_melpvsdem_blocks_s_1000$V3-inserts_melpvsdem_blocks_s_1000$V2)/sum1
sum(inserts_melpvsdem_blocks_s_1000$V3-inserts_melpvsdem_blocks_s_1000$V2)/275200000
nrow(inserts_melpvsdem_blocks_s_1000)/nrow1

###
sum(inserts_chavsera_blocks$V3-inserts_chavsera_blocks$V2)
sum(inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2)
mean(inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2)
median(inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2)

inserts_chavsera_blocks_s_1 <- subset(inserts_chavsera_blocks_s, inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2 ==1)
inserts_chavsera_blocks_s_20 <- subset(inserts_chavsera_blocks_s, inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2 <20)
inserts_chavsera_blocks_s_1000 <- subset(inserts_chavsera_blocks_s, inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2 >1000)

sum1 <- sum(inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2)

sum(inserts_chavsera_blocks_s_1$V3-inserts_chavsera_blocks_s_1$V2)/sum1
sum(inserts_chavsera_blocks_s_20$V3-inserts_chavsera_blocks_s_20$V2)/sum1
sum(inserts_chavsera_blocks_s_1000$V3-inserts_chavsera_blocks_s_1000$V2)/sum1

max((inserts_chavsera_blocks_s$V3-inserts_chavsera_blocks_s$V2))


###
sum(inserts_chavsmelp_blocks$V3-inserts_chavsmelp_blocks$V2)
sum(inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2)
mean(inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2)
median(inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2)

inserts_chavsmelp_blocks_s_1 <- subset(inserts_chavsmelp_blocks_s, inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2 ==1)
inserts_chavsmelp_blocks_s_20 <- subset(inserts_chavsmelp_blocks_s, inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2 <20)
inserts_chavsmelp_blocks_s_1000 <- subset(inserts_chavsmelp_blocks_s, inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2 >1000)

sum1 <- sum(inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2)

sum(inserts_chavsmelp_blocks_s_1$V3-inserts_chavsmelp_blocks_s_1$V2)/sum1
sum(inserts_chavsmelp_blocks_s_20$V3-inserts_chavsmelp_blocks_s_20$V2)/sum1
sum(inserts_chavsmelp_blocks_s_1000$V3-inserts_chavsmelp_blocks_s_1000$V2)/sum1

max((inserts_chavsmelp_blocks_s$V3-inserts_chavsmelp_blocks_s$V2))


##
sum(inserts_melpvschadem_blocks$V3-inserts_melpvschadem_blocks$V2)
sum(inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2)
mean(inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2)
median(inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2)

inserts_melpvschadem_blocks_s_1 <- subset(inserts_melpvschadem_blocks_s, inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2 ==1)
inserts_melpvschadem_blocks_s_20 <- subset(inserts_melpvschadem_blocks_s, inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2 <20)
inserts_melpvschadem_blocks_s_1000 <- subset(inserts_melpvschadem_blocks_s, inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2 >1000)

sum1 <- sum(inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2)

sum(inserts_melpvschadem_blocks_s_1$V3-inserts_melpvschadem_blocks_s_1$V2)/sum1
sum(inserts_melpvschadem_blocks_s_20$V3-inserts_melpvschadem_blocks_s_20$V2)/sum1
sum(inserts_melpvschadem_blocks_s_1000$V3-inserts_melpvschadem_blocks_s_1000$V2)/sum1

max((inserts_melpvschadem_blocks_s$V3-inserts_melpvschadem_blocks_s$V2))







inserts_chavscha_blocks$V4 <- max(inserts_cha2vscha1_blocks$V3-inserts_cha2vscha1_blocks$V2)
inserts_chavscha_blocks <- inserts_chavscha_blocks[order(inserts_chavscha_blocks$V4),]

sum((inserts_chavscha_blocks$V3-inserts_chavscha_blocks$V2))

sum1 <- sum(inserts_chavscha_blocks$V3-inserts_chavscha_blocks$V2)
inserts_chavscha_blocks_new <- subset(inserts_chavscha_blocks, inserts_chavscha_blocks$V3-inserts_chavscha_blocks$V2 >1000)

sum2 <- sum(inserts_chavscha_blocks_new$V3-inserts_chavscha_blocks_new$V2)

sum2/sum1



nrow(inserts_chavscha_blocks_new)/nrow(inserts_chavscha_blocks)

###############################
colnames(inserts_chavscha_blocks_s    ) <- c('scaffold', 'start', 'end')
colnames(inserts_chavsera_blocks_s    ) <- c('scaffold', 'start', 'end')
colnames(inserts_chavsmelp_blocks_s   ) <- c('scaffold', 'start', 'end')
colnames(inserts_melpvschadem_blocks_s) <- c('scaffold', 'start', 'end')

inserts_chavscha_blocks_s$insertL <- inserts_chavscha_blocks_s$end-inserts_chavscha_blocks_s$start
inserts_chavsera_blocks_s$insertL <- inserts_chavsera_blocks_s$end-inserts_chavsera_blocks_s$start
inserts_chavsmelp_blocks_s$insertL <- inserts_chavsmelp_blocks_s$end-inserts_chavsmelp_blocks_s$start
inserts_melpvschadem_blocks_s$insertL <- inserts_melpvschadem_blocks_s$end-inserts_melpvschadem_blocks_s$start

counts_chavscha <- c()
counts_chavsera <- c()
counts_chavsmelp <- c()
counts_melpvschadem <- c()

for(e in 1:21){
  table_sub <- subset(inserts_chavscha_blocks_s, inserts_chavscha_blocks_s$start >= chrom$start[e] & inserts_chavscha_blocks_s$end <= chrom$end[e])
  table_sub2 <- subset(table_sub, table_sub$insertL > 1 & table_sub$insertL < 100000)
  counts_chavscha <- rbind(counts_chavscha, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
}
for(e in 1:21){
  table_sub <- subset(inserts_chavsera_blocks_s, inserts_chavsera_blocks_s$start >= chrom$start[e] & inserts_chavsera_blocks_s$end <= chrom$end[e])
  table_sub2 <- subset(table_sub, table_sub$insertL > 1 & table_sub$insertL < 100000)
  counts_chavsera <- rbind(counts_chavsera, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
}
for(e in 1:21){
  table_sub <- subset(inserts_chavsmelp_blocks_s, inserts_chavsmelp_blocks_s$start >= chrom$start[e] & inserts_chavsmelp_blocks_s$end <= chrom$end[e])
  table_sub2 <- subset(table_sub, table_sub$insertL > 1 & table_sub$insertL < 100000)
  counts_chavsmelp <- rbind(counts_chavsmelp, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
}
for(e in 1:21){
  table_sub <- subset(inserts_melpvschadem_blocks_s, inserts_melpvschadem_blocks_s$start >= chrom$start[e] & inserts_melpvschadem_blocks_s$end <= chrom$end[e])
  table_sub2 <- subset(table_sub, table_sub$insertL > 1 & table_sub$insertL < 100000)
  counts_melpvschadem <- rbind(counts_melpvschadem, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
}

counts_chavscha <- as.data.frame(counts_chavscha)
counts_chavsera <- as.data.frame(counts_chavsera)
counts_chavsmelp <- as.data.frame(counts_chavsmelp)
counts_melpvschadem <- as.data.frame(counts_melpvschadem)

colnames(counts_chavscha) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')
colnames(counts_chavsera) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')
colnames(counts_chavsmelp) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')
colnames(counts_melpvschadem) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')

counts_chavscha <- cbind(counts_chavscha, chrom_coords_cha)
counts_chavsera <- cbind(counts_chavsera, chrom_coords_cha)
counts_chavsmelp <- cbind(counts_chavsmelp, chrom_coords_cha)
counts_melpvschadem <- cbind(counts_melpvschadem, chrom_coords_melp)


par(mfrow=c(2,2))

plot(counts_chavscha$chromLengths, counts_chavscha$insertL/counts_chavscha$chromLengths, pch = 19, cex=4, xlab = 'chromosome length', ylab = 'Proportion unique', ylim=c(0,1), main = 'indels charithona', xlim = c(8000000,22000000))
text(counts_chavscha$chromLengths, counts_chavscha$insertL/counts_chavscha$chromLengths, counts_chavscha$chromosome, col = 'white')
fit <- lm(counts_chavscha$insertL/counts_chavscha$chromLengths~counts_chavscha$chromLengths)
sumfit <- summary(fit)
abline(fit)
text(16000000, 0.9, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))

plot(counts_chavsera$chromLengths, counts_chavsera$insertL/counts_chavsera$chromLengths, pch = 19, cex=4, xlab = 'chromosome length', ylab = 'Proportion unique', ylim=c(0,1), main = 'unique cha vs erato', xlim = c(8000000,22000000))
text(counts_chavsera$chromLengths, counts_chavsera$insertL/counts_chavsera$chromLengths, counts_chavsera$chromosome, col = 'white')
fit <- lm(counts_chavsera$insertL/counts_chavsera$chromLengths~counts_chavsera$chromLengths)
sumfit <- summary(fit)
abline(fit)
text(16000000, 0.9, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))

plot(counts_chavsmelp$chromLengths, counts_chavsmelp$insertL/counts_chavsmelp$chromLengths, pch = 19, cex=4, xlab = 'chromosome length', ylab = 'Proportion unique', ylim=c(0,1), main = 'unique cha vs melp', xlim = c(8000000,22000000))
text(counts_chavsmelp$chromLengths, counts_chavsmelp$insertL/counts_chavsmelp$chromLengths, counts_chavsmelp$chromosome, col = 'white')
fit <- lm(counts_chavsmelp$insertL/counts_chavsmelp$chromLengths~counts_chavsmelp$chromLengths)
sumfit <- summary(fit)
abline(fit)
text(16000000, 0.9, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))

plot(counts_melpvschadem$chromLengths, counts_melpvschadem$insertL/counts_melpvschadem$chromLengths, pch = 19, cex=4, xlab = 'chromosome length', ylab = 'Proportion unique', ylim=c(0,1), main = 'unique melp vs cha', xlim = c(8000000,22000000))
text(counts_melpvschadem$chromLengths, counts_melpvschadem$insertL/counts_melpvschadem$chromLengths, counts_melpvschadem$chromosome, col = 'white')
fit <- lm(counts_melpvschadem$insertL/counts_melpvschadem$chromLengths~counts_melpvschadem$chromLengths)
sumfit <- summary(fit)
abline(fit)
text(16000000, 0.9, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))










# ## 
# test <- read.table("inserts/blocks_shared_2_3.txt", h=F)
# colnames(test) <- c('scaf','start','end')
# test$insertL <- test$end-test$start
# sum(test$insertL)



# ### test
# 
# angelo <- read.table("RM/H_e_dem_peaks_start_end_pan.rowsort.txt.final.bed", h=F)
# head(angelo)
# colnames(angelo) <- c('scaf','start','end')
# 
# counts <- c()
# for(e in 1:21){
#   table_sub <- subset(angelo, angelo$start >= chrom$start[e] & angelo$end <= chrom$end[e])
#   table_sub$insertL <- table_sub$end-table_sub$start
#   
#   table_sub2 <- subset(table_sub, table_sub$insertL > 10 & table_sub$insertL < 10000)
#   
#   counts <- rbind(counts, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
#   
# }
# counts
# 
# colnames(counts) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')
# 
# 
# counts <- as.data.frame(counts)
# colnames(counts) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')
# 
# counts <- cbind(counts, chrom_coords)
# 
# head(counts)
# 
# # plot(counts$chromlength, counts$count)
# # plot(counts$chromlength, counts$count/counts$chromlength)
# 
# plot(counts$chromLengths, counts$insertL/counts$chromLengths, pch = 19, cex=4, xlab = 'PAN chromosome length', ylab = 'Proportion unique')
# text(counts$chromLengths, counts$insertL/counts$chromLengths, counts$chromosome, col = 'white')
# 
# fit <- lm(counts$insertL/counts$chromLengths~counts$chromLengths)
# sumfit <- summary(fit)
# 
# abline(fit)
# 
# text(17000000, 0.35, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))

# dem <- read.table("inserts/1_blocks_intervals.txt", h=T)
# dem$len <- dem$end-dem$start
# sum(dem$len)
# 
# melp <- read.table("inserts/4_blocks_intervals.txt", h=T)
# melp$len <- melp$end-melp$start
# sum(melp$len)
# 
# cha1 <- read.table("inserts/2_blocks_intervals.txt", h=T)
# cha1$len <- cha1$end-cha1$start
# sum(cha1$len)
# 
# cha2 <- read.table("inserts/3_blocks_intervals.txt", h=T)
# cha2$len <- cha2$end-cha2$start
# sum(cha2$len)


### TEs
par(mfrow=c(1,3))
TE <- read.table("RM/H_e_dem_TE.bed", h=F)
colnames(TE) <- c('scaffold','start','end')
head(TE)

TE$insertL <- TE$end-TE$start
sum(TE$insertL)

TE2 <- merge(TE, scaf_coords2, by='scaffold')

head(TE2)

counts <- c()
for(e in 1:21){
  table_sub <- subset(TE2, TE2$chromosome == e)
  
  table_sub2 <- subset(table_sub, table_sub$insertL > 1 & table_sub$insertL < 100000)
  
  counts <- rbind(counts, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
  
}
counts <- as.data.frame(counts)
colnames(counts) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')

counts <- cbind(counts, chrom_coords)

head(counts)

plot(counts$chromLengths, counts$insertL/counts$chromLengths, pch = 19, cex=4, xlab = 'PAN chromosome length', ylab = 'Proportion unique', ylim=c(0.2,1))
text(counts$chromLengths, counts$insertL/counts$chromLengths, counts$chromosome, col = 'white')

fit <- lm(counts$insertL/counts$chromLengths~counts$chromLengths)
sumfit <- summary(fit)

abline(fit)
###

TE <- read.table("RM/H_e_dem_peaks_start_end_pan.rowsort.txt.final.bed", h=F)
colnames(TE) <- c('scaffold','start','end')
head(TE)

TE$insertL <- TE$end-TE$start
sum(TE$insertL)

# TE2 <- merge(TE, scaf_coords2, by='scaffold')

head(TE2)

counts <- c()
for(e in 1:21){
  table_sub <- subset(TE, TE$start >= chrom$start[e] & TE$end <= chrom$end[e])
  
  table_sub2 <- subset(table_sub, table_sub$insertL > 1 & table_sub$insertL < 100000)
  
  counts <- rbind(counts, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
  
}

counts <- as.data.frame(counts)
colnames(counts) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')

counts <- cbind(counts, chrom_coords)

head(counts)

plot(counts$chromLengths, counts$insertL/counts$chromLengths, pch = 19, cex=4, xlab = 'PAN chromosome length', ylab = 'Proportion unique', ylim=c(0.2,1))
text(counts$chromLengths, counts$insertL/counts$chromLengths, counts$chromosome, col = 'white')

fit <- lm(counts$insertL/counts$chromLengths~counts$chromLengths)
sumfit <- summary(fit)

abline(fit)

text(17000000, 0.4, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))
###



TE <- read.table("RM/H_e_dem_peaks_start_end_pan.rowsort.txt.final.bed", h=F)
colnames(TE) <- c('scaffold','start','end')
head(TE)

TE$insertL <- TE$end-TE$start
sum(TE$insertL)

# TE2 <- merge(TE, scaf_coords2, by='scaffold')

head(TE2)

counts <- c()
for(e in 1:21){
  table_sub <- subset(TE, TE$start >= chrom$start[e] & TE$end <= chrom$end[e])
  
  table_sub2 <- subset(table_sub, table_sub$insertL > 1 & table_sub$insertL < 100000)
  
  counts <- rbind(counts, c(chrom$chromosome[e], nrow(table_sub), chrom$length[e], sum(table_sub2$insertL)))
  
}

counts <- as.data.frame(counts)
colnames(counts) <- c('chromosome', 'count', 'chromlength_pan', 'insertL')

counts <- cbind(counts, chrom_coords)

head(counts)

plot(counts$chromlength_pan, counts$insertL/counts$chromlength_pan, pch = 19, cex=4, xlab = 'PAN chromosome length', ylab = 'Proportion unique', ylim=c(0.2,1))
text(counts$chromlength_pan, counts$insertL/counts$chromlength_pan, counts$chromosome, col = 'white')

fit <- lm(counts$insertL/counts$chromlength_pan~counts$chromlength_pan)
sumfit <- summary(fit)

abline(fit)

text(17000000, 0.4, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))
###








### dxy

dxyTable <- read.table('inserts/dxy_2_3.txt', h=T)
head(dxyTable)

dxyTable$dxy <- dxyTable$SNP/dxyTable$good_sites
mean(na.omit(dxyTable$dxy))

dxyMean <- c()
for(e in 1:21){
  
  table_sub <- subset(dxyTable, dxyTable$start >= chrom$start[e] & dxyTable$end <= chrom$end[e])
  table_sub2 <- subset(table_sub, table_sub$good_sites > 1000)
  
  dxyMean <- rbind(dxyMean, c(chrom$chromosome[e], chrom$length[e], mean(table_sub2$dxy)))
  
}

head(dxyMean)
mean(dxyMean$dxyM)

dxyMean <- as.data.frame(dxyMean)
colnames(dxyMean) <- c('chromosome', 'chromlength_pan', 'dxyM')

dxyMean <- cbind(dxyMean, chrom_coords_cha)


plot(dxyMean$chromlength_pan, dxyMean$dxyM, pch = 19, cex=4, xlab = 'chromosome length', ylab = 'dxy', ylim=c(0.08,0.16))
text(dxyMean$chromlength_pan, dxyMean$dxyM, counts_chavscha$chromosome, col = 'white')

fit <- lm(dxyMean$dxyM~dxyMean$chromlength_pan)
sumfit <- summary(fit)

abline(fit)

text(27000000, 0.15, paste('R2 = ', round(sumfit$r.squared,2), '; p = ', round(sumfit$coefficients[2,4], 6), sep = ''))

