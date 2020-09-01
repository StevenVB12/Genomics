###
# Functions
###

#calculate chromosome coordinates
chrom.coords <- function(scafL,chromNames,gap = 0) {
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
########
########




# calculate chromsome and scaf starts
scafL<-read.table("scaffold_lengths_chrom.txt", h=T)
chromStarts <- read.table('chrom_starts.txt', h = T)
chromNames <-c(1:21)
chrom_coords <- chrom.coords(scafL,chromNames)
scaf_coords <- scaf.coords(scafL)
scaf_coords2 <- merge(scafL,chrom_coords,by="chromosome", all.x=TRUE)
head(scaf_coords2)


# open gene density and merge to chromsome and scaf starts table
geneDens <- read.table('gene_density/Herato_gff_50k.txt', h=T)
geneDens <- merge(geneDens,scaf_coords2,by="scaffold", all.x=TRUE)
geneDens$GenomePos <- geneDens$position + geneDens$scafStart + geneDens$chromStarts -2


compsLD <- c("cyr", "dem", "emm", "era", "ety", "fav", "hydP", "hydFG", "lat", "not", "pet", "phy", "ven")

files_LDh <- list.files(path='LDhelmet', pattern="_ama", full.names = T)

LDh <- c()
for(f in 1:length(files_LDh)){
  table <- read.table(files_LDh[f], h=T)
  LDh <-rbind(LDh,table)
}
colnames(LDh) <- c('chromscaf','scaffold', 'position','Cposition','mean','p0.025','p0.500','p0.975')

# Merge LDhelmet table and chromosome positions
LDh_chrom <- merge(LDh,scaf_coords2, by="scaffold", all.x=TRUE)
LDh_chrom$GenomePos <- LDh_chrom$Cposition + LDh_chrom$scafStart + LDh_chrom$chromStarts -2


for(e in 1:length(compsLD)){
# read in LDhelmet files
files_LDh1 <- list.files(path='LDhelmet', pattern=paste('_',compsLD[e], sep=''), full.names = T)

LDh1 <- c()
for(f in 1:length(files_LDh1)){
  table <- read.table(files_LDh1[f], h=T)
  LDh1 <-rbind(LDh1,table)
}
colnames(LDh1) <- c('chromscaf','scaffold', 'position','Cposition','mean','p0.025','p0.500','p0.975')

# # read in LDhelmet files
# files_LDh2 <- list.files(path='LDhelmet', pattern='_fav', full.names = T)
# 
# LDh2 <- c()
# for(f in 1:length(files_LDh2)){
#   table <- read.table(files_LDh2[f], h=T)
#   LDh2 <-rbind(LDh2,table)
# }
# colnames(LDh2) <- c('chromscaf','scaffold', 'position','Cposition','mean','p0.025','p0.500','p0.975')
# 
# plot(LDh1$mean, LDh2$mean, xlim=c(0,0.6),ylim=c(0,0.6), pch=19, cex=0.5, col=adjustcolor("black", alpha.f = 0.1))
# 
# head(LDh1)
# head(geneDens)
# head(chromStarts)

# Merge LDhelmet table and chromosome positions
LDh_chrom1 <- merge(LDh1,scaf_coords2, by="scaffold", all.x=TRUE)
LDh_chrom1$GenomePos <- LDh_chrom1$Cposition + LDh_chrom1$scafStart + LDh_chrom1$chromStarts -2

# LDh_chrom2 <- merge(LDh2,scaf_coords2, by="scaffold", all.x=TRUE)
# LDh_chrom2$GenomePos <- LDh_chrom2$Cposition + LDh_chrom2$scafStart + LDh_chrom2$chromStarts -2

LDh_chrom <- merge(LDh_chrom,LDh_chrom1[,c(5,17)], by="GenomePos", all.x=TRUE)

# LDh_chrom$mean <- (LDh_chrom$mean.x +LDh_chrom$mean.y)/2
# 
# head(LDh_chrom)
}
# Read in cM table (calculated from linkage map)

meanT   <- LDh_chrom[ , grepl( "mean" , names( LDh_chrom ) ) ]
meanTm <- rowMeans(meanT, na.rm = F)

LDh_chrom$mean <- meanTm
head(LDh_chrom)

cM <- read.table('chrom_cM_length.txt', h=T)

# calculate scaler for cM/Mb transformation
for(e in 1:21){
  chr <- subset(LDh_chrom, LDh_chrom$chromosome == e)
  sumMb <- sum(chr$mean*50000)/1000000
  cM$scaler[e] <- cM$cM[e]/sumMb
}

LDh_chromS <- merge(LDh_chrom, cM, by="chromosome", all.x = TRUE)
LDh_chromS$meanS <- LDh_chromS$mean*LDh_chromS$scaler
head(LDh_chromS)

# outp <- LDh_chromS[,c(3,5,2,19,23)]
# write.table(outp, "LDhelmet_petiverana.txt", sep="\t", quote=F, row.names=F)


















# merge LDh table with gene density
LDh_genedens <- merge(geneDens,LDh_chromS,by="GenomePos", all.x=TRUE)
head(LDh_genedens)

# read in stats for hybrid zone
HZstats1 <- read.table("Fst_stats/himCYR_cyr.stats",h=T)
#HZstats <- read.table("Fst_stats/erato/emmFAV_favEMM.stats",h=T)
HZstats1 <- merge(HZstats1,scaf_coords2,by="scaffold", all.x=TRUE)
HZstats1$GenomePos <- HZstats1$position+HZstats1$scafStart+HZstats1$chromStart-3
head(HZstats1)

HZstats2 <- read.table("Fst_stats/emmHIM_himEM.stats",h=T)
HZstats2 <- merge(HZstats2,scaf_coords2,by="scaffold", all.x=TRUE)
HZstats2$GenomePos <- HZstats2$position+HZstats2$scafStart+HZstats2$chromStart-3
head(HZstats2)

HZstats3 <- read.table("Fst_stats/favEMM-himEM.stats",h=T)
HZstats3 <- merge(HZstats3,scaf_coords2,by="scaffold", all.x=TRUE)
HZstats3$GenomePos <- HZstats3$position+HZstats3$scafStart+HZstats3$chromStart-3
head(HZstats3)


LDh_stats1 <- merge(HZstats1,LDh_genedens,by="GenomePos", all.x=TRUE)
LDh_stats2 <- merge(HZstats2,LDh_genedens,by="GenomePos", all.x=TRUE)
LDh_stats3 <- merge(HZstats3,LDh_genedens,by="GenomePos", all.x=TRUE)

head(LDh_stats)

# read in D stats for hyrid zone
files1 <- list.files(path='ABBA_BABA', pattern='cyrN_cyrS_himN_her', full.names = T)
ABBA1 <- c()
for(f in 1:length(files1)){
  table <- read.csv(files1[f], h=T)
  ABBA1 <-rbind(ABBA1,table)
}
ABBA1 <- merge(ABBA1,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA1$GenomePos <- (ABBA1$start+ABBA1$end-1)/2 + ABBA1$scafStart + ABBA1$chromStart-2
LDh_final1 <- merge(LDh_stats1,ABBA1,by="GenomePos", all.x=TRUE)


files2 <- list.files(path='ABBA_BABA', pattern='emmE_emmW_himS_her', full.names = T)
ABBA2 <- c()
for(f in 1:length(files2)){
  table <- read.csv(files2[f], h=T)
  ABBA2 <-rbind(ABBA2,table)
}
ABBA2 <- merge(ABBA2,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA2$GenomePos <- (ABBA2$start+ABBA2$end-1)/2 + ABBA2$scafStart + ABBA2$chromStart-2
LDh_final2 <- merge(LDh_stats2,ABBA2,by="GenomePos", all.x=TRUE)

files3 <- list.files(path='ABBA_BABA', pattern='favE_favW_himS_her', full.names = T)
ABBA3 <- c()
for(f in 1:length(files3)){
  table <- read.csv(files3[f], h=T)
  ABBA3 <-rbind(ABBA3,table)
}
ABBA3 <- merge(ABBA3,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA3$GenomePos <- (ABBA3$start+ABBA3$end-1)/2 + ABBA3$scafStart + ABBA3$chromStart-2
LDh_final3 <- merge(LDh_stats3,ABBA3,by="GenomePos", all.x=TRUE)






par(mfcol=c(3,2),mai=c(0.2,0.2,0,0), oma=c(1,1,0.3,0.3)+0, by.row=F)

LDh_final1 <- subset(LDh_final1,LDh_final1$fd >=0)
LDh_final2 <- subset(LDh_final2,LDh_final2$fd >=0)
LDh_final3 <- subset(LDh_final3,LDh_final3$fd >=0)

plot(LDh_final1$exonP, LDh_final1$Fst, xlim=c(0,0.6), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt = 'none', xlab = "Proportion Coding sequence", ylab = "")
Axis(side=1, labels=FALSE)
# lo<-loess.smooth(LDh_final1$exonP, LDh_final1$Fst, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final1$Fst~LDh_final1$exonP ), lwd=2, col = 'green')
summary(lm(LDh_final1$Fst~LDh_final1$exonP))

plot(LDh_final2$exonP, LDh_final2$Fst, xlim=c(0,0.6), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt = 'none', xlab = "Proportion Coding sequence", ylab = "")
Axis(side=1, labels=FALSE)
# lo<-loess.smooth(LDh_final2$exonP, LDh_final2$Fst, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final2$Fst~LDh_final2$exonP ), lwd=2, col = 'green')
summary(lm(LDh_final2$Fst~LDh_final2$exonP))

plot(LDh_final3$exonP, LDh_final3$Fst, xlim=c(0,0.6), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xlab = "Proportion Coding sequence", ylab = "")
# lo<-loess.smooth(LDh_final3$exonP, LDh_final3$Fst, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final3$Fst~LDh_final3$exonP ), lwd=2, col = 'green')
summary(lm(LDh_final3$Fst~LDh_final3$exonP))

plot(LDh_final1$mean, LDh_final1$Fst, xlim=c(0,0.25), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt = 'none', yaxt = 'none', xlab = "Recombination rate (cM/Mb)", ylab = "")
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)
# lo<-loess.smooth(LDh_final1$meanS*2, LDh_final1$Fst, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final1$Fst~I(LDh_final1$mean) ), lwd=2, col = 'green')
summary(lm(LDh_final1$Fst~I(LDh_final1$mean) ))

plot(LDh_final2$mean, LDh_final2$Fst, xlim=c(0,0.25), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt = 'none', yaxt = 'none', xlab = "Recombination rate (cM/Mb)", ylab = "")
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)
# lo<-loess.smooth(LDh_final2$meanS*2, LDh_final2$Fst, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final2$Fst~I(LDh_final2$mean) ), lwd=2, col = 'green')
summary(lm(LDh_final2$Fst~I(LDh_final2$mean) ))

plot(LDh_final3$mean, LDh_final3$Fst, xlim=c(0,0.25), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, yaxt = 'none', xlab = "Recombination rate (cM/Mb)", ylab = "")
Axis(side=2, labels=FALSE)
# lo<-loess.smooth(LDh_final3$meanS*2, LDh_final3$Fst, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final3$Fst~I(LDh_final3$mean) ), lwd=2, col = 'green')
summary(lm(LDh_final3$Fst~I(LDh_final3$mean) ))




plot(LDh_final1$exonP, LDh_final1$fd, xlim=c(0,0.6), ylim=c(0,0.4), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xlab = "", ylab = "", xaxt='none')
Axis(side=1, labels=FALSE)
# lo<-loess.smooth(LDh_final1$exonP, LDh_final1$fd, span = 1/10, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final1$fd~LDh_final1$exonP ), lwd=2, col = 'green')
summary(lm(LDh_final1$fd~LDh_final1$exonP ))

plot(LDh_final2$exonP, LDh_final2$fd, xlim=c(0,0.6), ylim=c(0,0.4), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xlab = "", ylab = "", xaxt='none')
Axis(side=1, labels=FALSE)
# lo<-loess.smooth(LDh_final2$exonP, LDh_final2$fd, span = 1/10, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final2$fd~LDh_final2$exonP ), lwd=2, col = 'green')
summary(lm(LDh_final2$fd~LDh_final2$exonP ))

plot(LDh_final3$exonP, LDh_final3$fd, xlim=c(0,0.6), ylim=c(0,0.4), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xlab = "Proportion Coding sequence", ylab = "")
# lo<-loess.smooth(LDh_final3$exonP, LDh_final3$fd, span = 1/10, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final3$fd~LDh_final3$exonP ), lwd=2, col = 'green')
summary(lm(LDh_final3$fd~LDh_final3$exonP ))

plot(LDh_final1$mean, LDh_final1$fd, xlim=c(0,0.25), ylim=c(0,0.4), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt = 'none', yaxt = 'none', xlab = "", ylab = "fd")
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)
# lo<-loess.smooth(LDh_final1$mean, LDh_final1$fd, span = 1/10, degree = 1, 
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final1$fd~I(LDh_final1$mean) ), lwd=2, col = 'green')
summary(lm(LDh_final1$fd~I(LDh_final1$mean)))

plot(LDh_final2$mean, LDh_final2$fd, xlim=c(0,0.25), ylim=c(0,0.4), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt = 'none', yaxt = 'none', xlab = "", ylab = "fd")
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)
# lo<-loess.smooth(LDh_final2$mean, LDh_final2$fd, span = 1/10, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final2$fd~I(LDh_final2$mean) ), lwd=2, col = 'green')
summary(lm(LDh_final2$fd~I(LDh_final2$mean)))

plot(LDh_final3$mean, LDh_final3$fd, xlim=c(0,0.25), ylim=c(0,0.4), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, yaxt = 'none', xlab = "Recombination rate (cM/Mb)", ylab = "fd")
Axis(side=2, labels=FALSE)
# lo<-loess.smooth(LDh_final3$mean, LDh_final3$fd, span = 1/10, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final3$fd~I(LDh_final3$mean) ), lwd=2, col = 'green')
summary(lm(LDh_final3$fd~I(LDh_final3$mean)))





plot(LDh_final1$fd, LDh_final1$Fst, xlim=c(0,0.4), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt='none', xlab = "fd", ylab = "")
Axis(side=1, labels=FALSE)
# lo<-loess.smooth(LDh_final1$fd, LDh_final1$Fst, span = 1/1, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final1$Fst~LDh_final1$fd ), lwd=2, col = 'green')
summary(lm(LDh_final1$Fst~LDh_final1$fd ))

plot(LDh_final2$fd, LDh_final2$Fst, xlim=c(0,0.4), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xaxt='none', xlab = "fd", ylab = "")
Axis(side=1, labels=FALSE)
# lo<-loess.smooth(LDh_final2$fd, LDh_final2$Fst, span = 1/1, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final2$Fst~LDh_final2$fd ), lwd=2, col = 'green')
summary(lm(LDh_final2$Fst~LDh_final2$fd ))

plot(LDh_final3$fd, LDh_final3$Fst, xlim=c(0,0.4), ylim=c(0,1), pch=19, col=adjustcolor("black", alpha.f=0.2), cex=1, xlab = "fd", ylab = "")

# lo<-loess.smooth(LDh_final3$fd, LDh_final3$Fst, span = 1/1, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col = 'green', lwd =2)
abline(lm(LDh_final3$Fst~LDh_final3$fd ), lwd=2, col = 'green')
summary(lm(LDh_final3$Fst~LDh_final3$fd ))








