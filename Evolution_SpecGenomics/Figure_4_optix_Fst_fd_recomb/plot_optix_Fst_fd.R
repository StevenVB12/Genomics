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
scafL<-read.table("scaffold_lengths_chrom.txt",h=T)
chromNames <-c(1:21)
chrom_coords <- chrom.coords(scafL,chromNames)
scaf_coords <- scaf.coords(scafL)

scaf_coords2 <- merge(scafL,chrom_coords,by="chromosome", all.x=TRUE)
###
#Read data

EMM_FAV <- read.table("Fst_stats/emmFAV_favEMM.stats",h=T)
HIM_CYR <- read.table("Fst_stats/himCYR_cyr.stats",h=T)
EMM_HIM <- read.table("Fst_stats/emmHIM_himEM.stats",h=T)
FAV_HIM <- read.table("Fst_stats/favEMM-himEM.stats",h=T)

EMM_FAV <- merge(EMM_FAV,scaf_coords2,by="scaffold", all.x=TRUE)
HIM_CYR <- merge(HIM_CYR,scaf_coords2,by="scaffold", all.x=TRUE)
EMM_HIM <- merge(EMM_HIM,scaf_coords2,by="scaffold", all.x=TRUE)
FAV_HIM <- merge(FAV_HIM,scaf_coords2,by="scaffold", all.x=TRUE)

EMM_FAV <- cbind(EMM_FAV,chromPos=EMM_FAV$position+EMM_FAV$scafStart+EMM_FAV$chromStart-2)
HIM_CYR <- cbind(HIM_CYR,chromPos=HIM_CYR$position+HIM_CYR$scafStart+HIM_CYR$chromStart-2)
EMM_HIM <- cbind(EMM_HIM,chromPos=EMM_HIM$position+EMM_HIM$scafStart+EMM_HIM$chromStart-2)
FAV_HIM <- cbind(FAV_HIM,chromPos=FAV_HIM$position+FAV_HIM$scafStart+FAV_HIM$chromStart-2)


EMM_FAVsub <- subset(EMM_FAV, EMM_FAV$scaffold == 'Herato1801' & EMM_FAV$chromPos > 335236276+1239943-2000000 & EMM_FAV$chromPos < 335236276+1239943+2000000)
HIM_CYRsub <- subset(HIM_CYR, HIM_CYR$scaffold == 'Herato1801' & HIM_CYR$chromPos > 335236276+1239943-2000000 & HIM_CYR$chromPos < 335236276+1239943+2000000)
EMM_HIMsub <- subset(EMM_HIM, EMM_HIM$scaffold == 'Herato1801' & EMM_HIM$chromPos > 335236276+1239943-2000000 & EMM_HIM$chromPos < 335236276+1239943+2000000)
FAV_HIMsub <- subset(FAV_HIM, FAV_HIM$scaffold == 'Herato1801' & FAV_HIM$chromPos > 335236276+1239943-2000000 & FAV_HIM$chromPos < 335236276+1239943+2000000)

EMM_FAVsub <- EMM_FAVsub[order(EMM_FAVsub$chromPos),] 
HIM_CYRsub <- HIM_CYRsub[order(HIM_CYRsub$chromPos),] 
EMM_HIMsub <- EMM_HIMsub[order(EMM_HIMsub$chromPos),] 
FAV_HIMsub <- FAV_HIMsub[order(FAV_HIMsub$chromPos),] 


comp <- list(EMM_FAVsub,HIM_CYRsub,EMM_HIMsub,FAV_HIMsub)
cols <- c("blue", "red", "red", "red")
linet <- c(1,1,2,3)

layout(matrix(c(1:4), nrow=4, byrow=TRUE), heights=c(1,5,5,5))
layout.show(n=4)


par(mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

plot(NULL, xlim = c(0,1239943+2000000), ylim = c(0,1), axes=FALSE,ann=FALSE)

ANNOT <- read.table("GFF/Herato1801.gff",sep="\t")
names(ANNOT) <- c("contig", "HGC", "type", "con_start", "con_end", "dot", "str", "unk", "descr")
ANNOT <- subset(ANNOT,ANNOT$con_start > start & ANNOT$con_end < end)
for (g in 1:nrow(ANNOT)){
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "-") 
    rect(ANNOT$con_start[g], 0.75, ANNOT$con_end[g], 0.75, col = NULL, border = "black")
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "-") 
    rect(ANNOT$con_start[g], 0.5, ANNOT$con_end[g], 1, col = "black", border = "black", lwd=0.3)
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "+") 
    rect(ANNOT$con_start[g], 0.25, ANNOT$con_end[g], 0.25, col = NULL, border = "black")
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "+") 
    rect(ANNOT$con_start[g], 0, ANNOT$con_end[g], 0.5, col = "black", border = "black", lwd=0.3)
}
rect(1239943,0.5,1251211,1, col = "red", border='red')

plot(0, pch = "",xlim = c(1239943-1200000,1239943+2000000), ylim = c(0,1), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)


for(e in 1:length(comp)){
  par(new=T)
  plot(comp[[e]]$chromPos-335236276, comp[[e]]$Fst, type='l', lwd =2, col=cols[e], lty=linet[e],
       xlim = c(0,1239943+2000000), ylim = c(0,1),
       ylab = 'Fst', xlab = 'Position', xaxt='none')
}

###

files3_emm <- list.files(path='ABBA_BABA', pattern='emmE_emmW_himS_her', full.names = T)
files3_fav <- list.files(path='ABBA_BABA', pattern='favE_favW_himS_her', full.names = T)
files1_cyr <- list.files(path='ABBA_BABA', pattern='cyrN_cyrS_himN_her', full.names = T)
files1_emmfav <- list.files(path='ABBA_BABA', pattern='emmW_emmE_favE_her', full.names = T)


ABBA3_emm <- c()
ABBA3_fav <- c()
ABBA1_cyr <- c()
ABBA1_emmfav <- c()

for(f in 1:length(files3_emm)){
  table3 <- read.csv(files3_emm[f], h=T)
  ABBA3_emm <-rbind(ABBA3_emm,table3)
}
for(f in 1:length(files3_fav)){
  table3 <- read.csv(files3_fav[f], h=T)
  ABBA3_fav <-rbind(ABBA3_fav,table3)
}
for(f in 1:length(files1_cyr)){
  table1 <- read.csv(files1_cyr[f], h=T)
  ABBA1_cyr <-rbind(ABBA1_cyr,table1)
}
for(f in 1:length(files1_emmfav)){
  table1 <- read.csv(files1_emmfav[f], h=T)
  ABBA1_emmfav <-rbind(ABBA1_emmfav,table1)
}

ABBA3_emm <- merge(ABBA3_emm,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA3_emm$chromPos <- ABBA3_emm$mid+ABBA3_emm$scafStart+ABBA3_emm$chromStart-2

ABBA3_fav <- merge(ABBA3_fav,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA3_fav$chromPos <- ABBA3_fav$mid+ABBA3_fav$scafStart+ABBA3_fav$chromStart-2

ABBA1_cyr <- merge(ABBA1_cyr,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA1_cyr$chromPos <- ABBA1_cyr$mid+ABBA1_cyr$scafStart+ABBA1_cyr$chromStart-2

ABBA1_emmfav <- merge(ABBA1_emmfav,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA1_emmfav$chromPos <- ABBA1_emmfav$mid+ABBA1_emmfav$scafStart+ABBA1_emmfav$chromStart-2


ABBA3_emmsub <- subset(ABBA3_emm, ABBA3_emm$scaffold == 'Herato1801' & ABBA3_emm$chromPos > 335236276+1239943-2000000 & ABBA3_emm$chromPos < 335236276+1239943+2000000)
ABBA3_favsub <- subset(ABBA3_fav, ABBA3_fav$scaffold == 'Herato1801' & ABBA3_fav$chromPos > 335236276+1239943-2000000 & ABBA3_fav$chromPos < 335236276+1239943+2000000)
ABBA1_cyrsub <- subset(ABBA1_cyr, ABBA1_cyr$scaffold == 'Herato1801' & ABBA1_cyr$chromPos > 335236276+1239943-2000000 & ABBA1_cyr$chromPos < 335236276+1239943+2000000)
ABBA1_emmfavsub <- subset(ABBA1_emmfav, ABBA1_emmfav$scaffold == 'Herato1801' & ABBA1_emmfav$chromPos > 335236276+1239943-2000000 & ABBA1_emmfav$chromPos < 335236276+1239943+2000000)
# ABBA1_emmfavsub$fd <- -ABBA1_emmfavsub$fd
comp <- list(ABBA1_emmfavsub,ABBA3_emmsub,ABBA3_favsub,ABBA1_cyrsub)
# cols <- c("red", "red", "red")
# linet <- c(1,2,3)


plot(0, pch = "",xlim = c(1239943-1200000,1239943+2000000), ylim = c(0,1), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
# rect(335236276+1239943,0.9,335236276+1251211,1, col = "black")
for(e in 1:length(comp)){
  par(new=T)
  plot(comp[[e]]$chromPos-335236276, comp[[e]]$fd, type='l', lwd =2, col=cols[e], lty=linet[e],
       xlim = c(1239943-1200000,1239943+2000000), ylim = c(0,1),
       ylab = 'fd', xlab = 'Position', xaxt='none')
}











####


compsLD <- c("cyr", "dem", "emm", "era", "ety", "fav", "hydP", "hydFG", "lat", "not", "pet", "phy", "ven")

files_LDh <- list.files(path='LDhelmet', pattern="_ama", full.names = T)

LDh <- c()
for(f in 1:length(files_LDh)){
  table <- read.table(files_LDh[f], h=T)
  LDh <-rbind(LDh,table)
}
colnames(LDh) <- c('chromscaf','scaffold', 'position','Cposition',paste('mean_','ama',sep=''),'p0.025','p0.500','p0.975')

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
  colnames(LDh1) <- c('chromscaf','scaffold', 'position','Cposition',paste('mean_',compsLD[e],sep=''),'p0.025','p0.500','p0.975')
  
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
meanT   <- LDh_chrom[ , grepl( "mean" , names( LDh_chrom ) ) ]
meanTm <- rowMeans(meanT, na.rm = F)

LDh_chrom$mean <- meanTm

LDh_chrom <- subset(LDh_chrom, LDh_chrom$scaffold == 'Herato1801' & LDh_chrom$GenomePos > 335236276+1239943-2000000 & LDh_chrom$GenomePos < 335236276+1239943+2000000)

plot(LDh_chrom$GenomePos-335236276, LDh_chrom$mean, type='l', lwd =2, col='black', lty=1,
     xlim = c(1239943-1200000,1239943+2000000), ylim = c(0.0,0.1),
     ylab = 'Recomb', xlab = 'Position')

lo<-loess.smooth(LDh_chrom$GenomePos-335236276, LDh_chrom$mean, span = 1/2, degree = 1,
family = c("symmetric", "gaussian"), evaluation = 50)
lines(lo, col = "green", lwd =2)


# open gene density and merge to chromsome and scaf starts table
geneDens <- read.table('gene_density/Herato_gff_50k.txt', h=T)
geneDens <- merge(geneDens,scaf_coords2,by="scaffold", all.x=TRUE)
geneDens$GenomePos <- geneDens$position + geneDens$scafStart + geneDens$chromStarts -2
geneDens <- subset(geneDens, geneDens$scaffold == 'Herato1801' & geneDens$GenomePos > 335236276+1239943-2000000 & geneDens$GenomePos < 335236276+1239943+2000000)

par(new=T)
plot(geneDens$GenomePos-335236276, geneDens$exonP, type='l', lwd =2, col='black', lty=3,
     xlim = c(1239943-1200000,1239943+2000000), ylim = c(0,1),
     ylab = 'Recomb', xlab = 'Position', yaxt='none')
lo<-loess.smooth(geneDens$GenomePos-335236276, geneDens$exonP, span = 1/2, degree = 1,
                 family = c("symmetric", "gaussian"), evaluation = 50)
lines(lo, col = "blue", lwd =2)

axis(4, at=seq(0,1,by=0.2))
