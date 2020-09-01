library(tidyr)
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


####################################################################################
scafL<-read.table("scaffold_lengths_chrom.txt",h=T)
chromNames <-c(1:21)
chrom_coords <- chrom.coords(scafL,chromNames)
scaf_coords <- scaf.coords(scafL)

scaf_coords2 <- merge(scafL,chrom_coords,by="chromosome", all.x=TRUE)
####################################################################################
# 
# files1_emm <- list.files(path='ABBA_BABA_new', pattern='phy_himN_emmW_her', full.names = T)
# files2_emm <- list.files(path='ABBA_BABA_new', pattern='phy_himN_emmE_her', full.names = T)
files3_emm <- list.files(path='ABBA_BABA', pattern='emmE_emmW_himS_her', full.names = T)
# files3_emm <- list.files(path='ABBA_BABA_new', pattern='himS_himN_emmW_her', full.names = T)
# files4_emm <- list.files(path='ABBA_BABA_new', pattern='phy_himS_emmE_her', full.names = T)
# 
# files1_fav <- list.files(path='ABBA_BABA_new', pattern='phy_himN_favWs_her', full.names = T)
# files2_fav <- list.files(path='ABBA_BABA_new', pattern='phy_himN_favE_her', full.names = T)
files3_fav <- list.files(path='ABBA_BABA', pattern='favE_favW_himS_her', full.names = T)
# files3_fav <- list.files(path='ABBA_BABA_new', pattern='himS_himN_favW_her', full.names = T)
# files4_fav <- list.files(path='ABBA_BABA_new', pattern='phy_himS_favE_her', full.names = T)

files1_cyr <- list.files(path='ABBA_BABA', pattern='cyrN_cyrS_himN_her', full.names = T)
# files1_cyr <- list.files(path='ABBA_BABA_new', pattern='himS_himN_cyrS_her', full.names = T)
# files2_cyr <- list.files(path='ABBA_BABA_new', pattern='phy_himN_cyrN_her', full.names = T)
# files3_cyr <- list.files(path='ABBA_BABA_new', pattern='phy_himS_cyrS_her', full.names = T)
# files4_cyr <- list.files(path='ABBA_BABA_new', pattern='phy_himS_cyrN_her', full.names = T)
# 
# files_EW <- list.files(path='ABBA_BABA_new', pattern='ven_west_east_her', full.names = T)

# files1 <- list.files(path='ABBA_BABA_new', pattern='himS_himN_cyrS_her', full.names = T)
# files2 <- list.files(path='ABBA_BABA_new', pattern='himS_himN_cyrN_her', full.names = T)

# files1 <- list.files(path='ABBA_BABA_new', pattern='himS_himN_emmW_her', full.names = T)
# files2 <- list.files(path='ABBA_BABA_new', pattern='himS_himN_emmE_her', full.names = T)

# filesX1 <- list.files(path='ABBA_BABA_new', pattern='himS_himN_cyrS_her', full.names = T)
# filesX2 <- list.files(path='ABBA_BABA_new', pattern='himS_himN_emmW_her', full.names = T)
# filesX3 <- list.files(path='ABBA_BABA_new', pattern='himS_himN_favW_her', full.names = T)

ABBA1_emm <- c()
ABBA2_emm <- c()
ABBA3_emm <- c()
ABBA4_emm <- c()

ABBA1_fav <- c()
ABBA2_fav <- c()
ABBA3_fav <- c()
ABBA4_fav <- c()

ABBA1_cyr <- c()
ABBA2_cyr <- c()
ABBA3_cyr <- c()
ABBA4_cyr <- c()

ABBA_EW <- c()

for(f in 1:length(files1_emm)){
  table1 <- read.csv(files1_emm[f], h=T)
  ABBA1_emm <-rbind(ABBA1_emm,table1)
}

for(f in 1:length(files2_emm)){
  table2 <- read.csv(files2_emm[f], h=T)
  ABBA2_emm <-rbind(ABBA2_emm,table2)
}

for(f in 1:length(files3_emm)){
  table3 <- read.csv(files3_emm[f], h=T)
  ABBA3_emm <-rbind(ABBA3_emm,table3)
}

for(f in 1:length(files4_emm)){
  table4 <- read.csv(files4_emm[f], h=T)
  ABBA4_emm <-rbind(ABBA4_emm,table4)
}


for(f in 1:length(files1_fav)){
  table1 <- read.csv(files1_fav[f], h=T)
  ABBA1_fav <-rbind(ABBA1_fav,table1)
}

for(f in 1:length(files2_fav)){
  table2 <- read.csv(files2_fav[f], h=T)
  ABBA2_fav <-rbind(ABBA2_fav,table2)
}

for(f in 1:length(files3_fav)){
  table3 <- read.csv(files3_fav[f], h=T)
  ABBA3_fav <-rbind(ABBA3_fav,table3)
}

for(f in 1:length(files4_fav)){
  table4 <- read.csv(files4_fav[f], h=T)
  ABBA4_fav <-rbind(ABBA4_fav,table4)
}


for(f in 1:length(files1_cyr)){
  table1 <- read.csv(files1_cyr[f], h=T)
  ABBA1_cyr <-rbind(ABBA1_cyr,table1)
}

for(f in 1:length(files2_cyr)){
  table2 <- read.csv(files2_cyr[f], h=T)
  ABBA2_cyr <-rbind(ABBA2_cyr,table2)
}

for(f in 1:length(files3_cyr)){
  table3 <- read.csv(files3_cyr[f], h=T)
  ABBA3_cyr <-rbind(ABBA3_cyr,table3)
}

for(f in 1:length(files4_cyr)){
  table4 <- read.csv(files4_cyr[f], h=T)
  ABBA4_cyr <-rbind(ABBA4_cyr,table4)
}

for(f in 1:length(files_EW)){
  table <- read.csv(files_EW[f], h=T)
  ABBA_EW <-rbind(ABBA_EW,table)
}

ABBA1_emm <- merge(ABBA1_emm,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA2_emm <- merge(ABBA2_emm,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA3_emm <- merge(ABBA3_emm,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA4_emm <- merge(ABBA4_emm,scaf_coords2,by="scaffold", all.x=TRUE)

ABBA1_emm$chromPos <- ABBA1_emm$mid+ABBA1_emm$scafStart+ABBA1_emm$chromStart-2
ABBA2_emm$chromPos <- ABBA2_emm$mid+ABBA2_emm$scafStart+ABBA2_emm$chromStart-2
ABBA3_emm$chromPos <- ABBA3_emm$mid+ABBA3_emm$scafStart+ABBA3_emm$chromStart-2
ABBA4_emm$chromPos <- ABBA4_emm$mid+ABBA4_emm$scafStart+ABBA4_emm$chromStart-2

ABBA1_emm$windID <-  paste(ABBA1_emm$scaffold,"_",ABBA1_emm$start, "_", ABBA1_emm$end, sep="")
ABBA2_emm$windID <-  paste(ABBA2_emm$scaffold,"_",ABBA2_emm$start, "_", ABBA2_emm$end, sep="")
ABBA3_emm$windID <-  paste(ABBA3_emm$scaffold,"_",ABBA3_emm$start, "_", ABBA3_emm$end, sep="")
ABBA4_emm$windID <-  paste(ABBA4_emm$scaffold,"_",ABBA4_emm$start, "_", ABBA4_emm$end, sep="")


ABBA1_fav <- merge(ABBA1_fav,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA2_fav <- merge(ABBA2_fav,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA3_fav <- merge(ABBA3_fav,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA4_fav <- merge(ABBA4_fav,scaf_coords2,by="scaffold", all.x=TRUE)

ABBA1_fav$chromPos <- ABBA1_fav$mid+ABBA1_fav$scafStart+ABBA1_fav$chromStart-2
ABBA2_fav$chromPos <- ABBA2_fav$mid+ABBA2_fav$scafStart+ABBA2_fav$chromStart-2
ABBA3_fav$chromPos <- ABBA3_fav$mid+ABBA3_fav$scafStart+ABBA3_fav$chromStart-2
ABBA4_fav$chromPos <- ABBA4_fav$mid+ABBA4_fav$scafStart+ABBA4_fav$chromStart-2

ABBA1_fav$windID <-  paste(ABBA1_fav$scaffold,"_",ABBA1_fav$start, "_", ABBA1_fav$end, sep="")
ABBA2_fav$windID <-  paste(ABBA2_fav$scaffold,"_",ABBA2_fav$start, "_", ABBA2_fav$end, sep="")
ABBA3_fav$windID <-  paste(ABBA3_fav$scaffold,"_",ABBA3_fav$start, "_", ABBA3_fav$end, sep="")
ABBA4_fav$windID <-  paste(ABBA4_fav$scaffold,"_",ABBA4_fav$start, "_", ABBA4_fav$end, sep="")


ABBA1_cyr <- merge(ABBA1_cyr,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA2_cyr <- merge(ABBA2_cyr,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA3_cyr <- merge(ABBA3_cyr,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA4_cyr <- merge(ABBA4_cyr,scaf_coords2,by="scaffold", all.x=TRUE)

ABBA1_cyr$chromPos <- ABBA1_cyr$mid+ABBA1_cyr$scafStart+ABBA1_cyr$chromStart-2
ABBA2_cyr$chromPos <- ABBA2_cyr$mid+ABBA2_cyr$scafStart+ABBA2_cyr$chromStart-2
ABBA3_cyr$chromPos <- ABBA3_cyr$mid+ABBA3_cyr$scafStart+ABBA3_cyr$chromStart-2
ABBA4_cyr$chromPos <- ABBA4_cyr$mid+ABBA4_cyr$scafStart+ABBA4_cyr$chromStart-2

ABBA1_cyr$windID <-  paste(ABBA1_cyr$scaffold,"_",ABBA1_cyr$start, "_", ABBA1_cyr$end, sep="")
ABBA2_cyr$windID <-  paste(ABBA2_cyr$scaffold,"_",ABBA2_cyr$start, "_", ABBA2_cyr$end, sep="")
ABBA3_cyr$windID <-  paste(ABBA3_cyr$scaffold,"_",ABBA3_cyr$start, "_", ABBA3_cyr$end, sep="")
ABBA4_cyr$windID <-  paste(ABBA4_cyr$scaffold,"_",ABBA4_cyr$start, "_", ABBA4_cyr$end, sep="")

ABBA_EW <- merge(ABBA_EW,scaf_coords2,by="scaffold", all.x=TRUE)
ABBA_EW$chromPos <- ABBA_EW$mid+ABBA_EW$scafStart+ABBA_EW$chromStart-2
ABBA_EW$windID <-  paste(ABBA_EW$scaffold,"_",ABBA_EW$start, "_", ABBA_EW$end, sep="")

par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)
top = 1
bot = -1

begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[21]/1000000

plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray85","gray80"), lwd = 0, border = c("gray85","gray80"))
  
par(new=TRUE)
plot(ABBA_EW$chromPos/1000000,ABBA_EW$D, type="p",pch=16, cex=0.8,col=adjustcolor('black',alpha.f = 0.5),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  
axis(2,cex.axis = 1, line = -1.5)

mtext(side = 2, text = expression(paste("D")), cex=1,line = 0.8)

axis(1, at=chrom_coords[,5][1:21]/1000000, labels=(1:21),lwd=0, lwd.ticks=0)

#####
DFOIL_emm <- list.files(path='DFOIL/him_emm', pattern='DFOIL', full.names = T)
DFOIL_emm <- DFOIL_emm[!grepl('BC2639', DFOIL_emm)]
DF_emm <- read.table(DFOIL_emm[1], h = T, comment.char = "")

DF_emm$X.chrom <- gsub('/home/rpapa/sbelleghem/work/DFOIL_out/him_emm/', '', DF_emm$X.chrom)
DF_emm <- separate(DF_emm, X.chrom, into = c('wind','indiv'), sep="_himera")
DF_emm <- separate(DF_emm, wind, into = c('scaffold','startDF','endDF'), sep="_")

DF_emm <- DF_emm[order(DF_emm$scaffold, as.integer(DF_emm$start)),  ]
DF_emm$windID <-  paste(DF_emm$scaffold,"_",DF_emm$startDF, "_", DF_emm$endDF, sep="")

ABBA1DF_emm <- merge(ABBA1_emm, DF_emm, by = 'windID')
ABBA2DF_emm <- merge(ABBA2_emm, DF_emm, by = 'windID')
ABBA3DF_emm <- merge(ABBA3_emm, DF_emm, by = 'windID')
ABBA4DF_emm <- merge(ABBA4_emm, DF_emm, by = 'windID')

for(e in 2:length(DFOIL_emm)){
  DF_emm <- read.table(DFOIL_emm[e], h = T, comment.char = "")
  
  DF_emm$X.chrom <- gsub('/home/rpapa/sbelleghem/work/DFOIL_out/him_emm/', '', DF_emm$X.chrom)
  DF_emm <- separate(DF_emm, X.chrom, into = c('wind','indiv'), sep="_himera")
  DF_emm <- separate(DF_emm, wind, into = c('scaffold','startDF','endDF'), sep="_")
  
  DF_emm <- DF_emm[order(DF_emm$scaffold, as.integer(DF_emm$start)),  ]
  DF_emm$windID <-  paste(DF_emm$scaffold,"_",DF_emm$startDF, "_", DF_emm$endDF, sep="")
  
  ABBA1DF2_emm <- merge(ABBA1_emm, DF_emm, by = 'windID')
  ABBA2DF2_emm <- merge(ABBA2_emm, DF_emm, by = 'windID')
  ABBA3DF2_emm <- merge(ABBA3_emm, DF_emm, by = 'windID')
  ABBA4DF2_emm <- merge(ABBA4_emm, DF_emm, by = 'windID')
  
  ABBA1DF_emm$introg13 <- ABBA1DF_emm$introg13 + ABBA1DF2_emm$introg13
  ABBA1DF_emm$introg14 <- ABBA1DF_emm$introg14 + ABBA1DF2_emm$introg14
  ABBA1DF_emm$introg23 <- ABBA1DF_emm$introg23 + ABBA1DF2_emm$introg23
  ABBA1DF_emm$introg24 <- ABBA1DF_emm$introg24 + ABBA1DF2_emm$introg24
  ABBA1DF_emm$introg31 <- ABBA1DF_emm$introg31 + ABBA1DF2_emm$introg31
  ABBA1DF_emm$introg41 <- ABBA1DF_emm$introg41 + ABBA1DF2_emm$introg41
  ABBA1DF_emm$introg32 <- ABBA1DF_emm$introg32 + ABBA1DF2_emm$introg32
  ABBA1DF_emm$introg42 <- ABBA1DF_emm$introg42 + ABBA1DF2_emm$introg42
  ABBA1DF_emm$introg123 <- ABBA1DF_emm$introg123 + ABBA1DF2_emm$introg123
  ABBA1DF_emm$introg124 <- ABBA1DF_emm$introg124 + ABBA1DF2_emm$introg124
  ABBA1DF_emm$intrognone <- ABBA1DF_emm$intrognone + ABBA1DF2_emm$intrognone
  
  ABBA2DF_emm$introg13 <- ABBA2DF_emm$introg13 + ABBA2DF2_emm$introg13
  ABBA2DF_emm$introg14 <- ABBA2DF_emm$introg14 + ABBA2DF2_emm$introg14
  ABBA2DF_emm$introg23 <- ABBA2DF_emm$introg23 + ABBA2DF2_emm$introg23
  ABBA2DF_emm$introg24 <- ABBA2DF_emm$introg24 + ABBA2DF2_emm$introg24
  ABBA2DF_emm$introg31 <- ABBA2DF_emm$introg31 + ABBA2DF2_emm$introg31
  ABBA2DF_emm$introg41 <- ABBA2DF_emm$introg41 + ABBA2DF2_emm$introg41
  ABBA2DF_emm$introg32 <- ABBA2DF_emm$introg32 + ABBA2DF2_emm$introg32
  ABBA2DF_emm$introg42 <- ABBA2DF_emm$introg42 + ABBA2DF2_emm$introg42
  ABBA2DF_emm$introg123 <- ABBA2DF_emm$introg123 + ABBA2DF2_emm$introg123
  ABBA2DF_emm$introg124 <- ABBA2DF_emm$introg124 + ABBA2DF2_emm$introg124
  ABBA2DF_emm$intrognone <- ABBA2DF_emm$intrognone + ABBA2DF2_emm$intrognone
  
  ABBA3DF_emm$introg13 <- ABBA3DF_emm$introg13 + ABBA3DF2_emm$introg13
  ABBA3DF_emm$introg14 <- ABBA3DF_emm$introg14 + ABBA3DF2_emm$introg14
  ABBA3DF_emm$introg23 <- ABBA3DF_emm$introg23 + ABBA3DF2_emm$introg23
  ABBA3DF_emm$introg24 <- ABBA3DF_emm$introg24 + ABBA3DF2_emm$introg24
  ABBA3DF_emm$introg31 <- ABBA3DF_emm$introg31 + ABBA3DF2_emm$introg31
  ABBA3DF_emm$introg41 <- ABBA3DF_emm$introg41 + ABBA3DF2_emm$introg41
  ABBA3DF_emm$introg32 <- ABBA3DF_emm$introg32 + ABBA3DF2_emm$introg32
  ABBA3DF_emm$introg42 <- ABBA3DF_emm$introg42 + ABBA3DF2_emm$introg42
  ABBA3DF_emm$introg123 <- ABBA3DF_emm$introg123 + ABBA3DF2_emm$introg123
  ABBA3DF_emm$introg124 <- ABBA3DF_emm$introg124 + ABBA3DF2_emm$introg124
  ABBA3DF_emm$intrognone <- ABBA3DF_emm$intrognone + ABBA3DF2_emm$intrognone
  
  ABBA4DF_emm$introg13 <- ABBA4DF_emm$introg13 + ABBA4DF2_emm$introg13
  ABBA4DF_emm$introg14 <- ABBA4DF_emm$introg14 + ABBA4DF2_emm$introg14
  ABBA4DF_emm$introg23 <- ABBA4DF_emm$introg23 + ABBA4DF2_emm$introg23
  ABBA4DF_emm$introg24 <- ABBA4DF_emm$introg24 + ABBA4DF2_emm$introg24
  ABBA4DF_emm$introg31 <- ABBA4DF_emm$introg31 + ABBA4DF2_emm$introg31
  ABBA4DF_emm$introg41 <- ABBA4DF_emm$introg41 + ABBA4DF2_emm$introg41
  ABBA4DF_emm$introg32 <- ABBA4DF_emm$introg32 + ABBA4DF2_emm$introg32
  ABBA4DF_emm$introg42 <- ABBA4DF_emm$introg42 + ABBA4DF2_emm$introg42
  ABBA4DF_emm$introg123 <- ABBA4DF_emm$introg123 + ABBA4DF2_emm$introg123
  ABBA4DF_emm$introg124 <- ABBA4DF_emm$introg124 + ABBA4DF2_emm$introg124
  ABBA4DF_emm$intrognone <- ABBA4DF_emm$intrognone + ABBA4DF2_emm$intrognone
}

DFOIL_fav <- list.files(path='DFOIL/him_fav', pattern='DFOIL', full.names = T)
DFOIL_fav <- DFOIL_fav[!grepl('BC2639', DFOIL_fav)]
DF_fav <- read.table(DFOIL_fav[1], h = T, comment.char = "")

DF_fav$X.chrom <- gsub('/home/rpapa/sbelleghem/work/DFOIL_out/him_fav/', '', DF_fav$X.chrom)
DF_fav <- separate(DF_fav, X.chrom, into = c('wind','indiv'), sep="_himera")
DF_fav <- separate(DF_fav, wind, into = c('scaffold','startDF','endDF'), sep="_")

DF_fav <- DF_fav[order(DF_fav$scaffold, as.integer(DF_fav$start)),  ]
DF_fav$windID <-  paste(DF_fav$scaffold,"_",DF_fav$startDF, "_", DF_fav$endDF, sep="")

ABBA1DF_fav <- merge(ABBA1_fav, DF_fav, by = 'windID')
ABBA2DF_fav <- merge(ABBA2_fav, DF_fav, by = 'windID')
ABBA3DF_fav <- merge(ABBA3_fav, DF_fav, by = 'windID')
ABBA4DF_fav <- merge(ABBA4_fav, DF_fav, by = 'windID')

for(e in 2:length(DFOIL_fav)){
  DF_fav <- read.table(DFOIL_fav[e], h = T, comment.char = "")
  
  DF_fav$X.chrom <- gsub('/home/rpapa/sbelleghem/work/DFOIL_out/him_fav/', '', DF_fav$X.chrom)
  DF_fav <- separate(DF_fav, X.chrom, into = c('wind','indiv'), sep="_himera")
  DF_fav <- separate(DF_fav, wind, into = c('scaffold','startDF','endDF'), sep="_")
  
  DF_fav <- DF_fav[order(DF_fav$scaffold, as.integer(DF_fav$start)),  ]
  DF_fav$windID <-  paste(DF_fav$scaffold,"_",DF_fav$startDF, "_", DF_fav$endDF, sep="")
  
  ABBA1DF2_fav <- merge(ABBA1_fav, DF_fav, by = 'windID')
  ABBA2DF2_fav <- merge(ABBA2_fav, DF_fav, by = 'windID')
  ABBA3DF2_fav <- merge(ABBA3_fav, DF_fav, by = 'windID')
  ABBA4DF2_fav <- merge(ABBA4_fav, DF_fav, by = 'windID')
  
  ABBA1DF_fav$introg13 <- ABBA1DF_fav$introg13 + ABBA1DF2_fav$introg13
  ABBA1DF_fav$introg14 <- ABBA1DF_fav$introg14 + ABBA1DF2_fav$introg14
  ABBA1DF_fav$introg23 <- ABBA1DF_fav$introg23 + ABBA1DF2_fav$introg23
  ABBA1DF_fav$introg24 <- ABBA1DF_fav$introg24 + ABBA1DF2_fav$introg24
  ABBA1DF_fav$introg31 <- ABBA1DF_fav$introg31 + ABBA1DF2_fav$introg31
  ABBA1DF_fav$introg41 <- ABBA1DF_fav$introg41 + ABBA1DF2_fav$introg41
  ABBA1DF_fav$introg32 <- ABBA1DF_fav$introg32 + ABBA1DF2_fav$introg32
  ABBA1DF_fav$introg42 <- ABBA1DF_fav$introg42 + ABBA1DF2_fav$introg42
  ABBA1DF_fav$introg123 <- ABBA1DF_fav$introg123 + ABBA1DF2_fav$introg123
  ABBA1DF_fav$introg124 <- ABBA1DF_fav$introg124 + ABBA1DF2_fav$introg124
  ABBA1DF_fav$intrognone <- ABBA1DF_fav$intrognone + ABBA1DF2_fav$intrognone
  
  ABBA2DF_fav$introg13 <- ABBA2DF_fav$introg13 + ABBA2DF2_fav$introg13
  ABBA2DF_fav$introg14 <- ABBA2DF_fav$introg14 + ABBA2DF2_fav$introg14
  ABBA2DF_fav$introg23 <- ABBA2DF_fav$introg23 + ABBA2DF2_fav$introg23
  ABBA2DF_fav$introg24 <- ABBA2DF_fav$introg24 + ABBA2DF2_fav$introg24
  ABBA2DF_fav$introg31 <- ABBA2DF_fav$introg31 + ABBA2DF2_fav$introg31
  ABBA2DF_fav$introg41 <- ABBA2DF_fav$introg41 + ABBA2DF2_fav$introg41
  ABBA2DF_fav$introg32 <- ABBA2DF_fav$introg32 + ABBA2DF2_fav$introg32
  ABBA2DF_fav$introg42 <- ABBA2DF_fav$introg42 + ABBA2DF2_fav$introg42
  ABBA2DF_fav$introg123 <- ABBA2DF_fav$introg123 + ABBA2DF2_fav$introg123
  ABBA2DF_fav$introg124 <- ABBA2DF_fav$introg124 + ABBA2DF2_fav$introg124
  ABBA2DF_fav$intrognone <- ABBA2DF_fav$intrognone + ABBA2DF2_fav$intrognone
  
  ABBA3DF_fav$introg13 <- ABBA3DF_fav$introg13 + ABBA3DF2_fav$introg13
  ABBA3DF_fav$introg14 <- ABBA3DF_fav$introg14 + ABBA3DF2_fav$introg14
  ABBA3DF_fav$introg23 <- ABBA3DF_fav$introg23 + ABBA3DF2_fav$introg23
  ABBA3DF_fav$introg24 <- ABBA3DF_fav$introg24 + ABBA3DF2_fav$introg24
  ABBA3DF_fav$introg31 <- ABBA3DF_fav$introg31 + ABBA3DF2_fav$introg31
  ABBA3DF_fav$introg41 <- ABBA3DF_fav$introg41 + ABBA3DF2_fav$introg41
  ABBA3DF_fav$introg32 <- ABBA3DF_fav$introg32 + ABBA3DF2_fav$introg32
  ABBA3DF_fav$introg42 <- ABBA3DF_fav$introg42 + ABBA3DF2_fav$introg42
  ABBA3DF_fav$introg123 <- ABBA3DF_fav$introg123 + ABBA3DF2_fav$introg123
  ABBA3DF_fav$introg124 <- ABBA3DF_fav$introg124 + ABBA3DF2_fav$introg124
  ABBA3DF_fav$intrognone <- ABBA3DF_fav$intrognone + ABBA3DF2_fav$intrognone
  
  ABBA4DF_fav$introg13 <- ABBA4DF_fav$introg13 + ABBA4DF2_fav$introg13
  ABBA4DF_fav$introg14 <- ABBA4DF_fav$introg14 + ABBA4DF2_fav$introg14
  ABBA4DF_fav$introg23 <- ABBA4DF_fav$introg23 + ABBA4DF2_fav$introg23
  ABBA4DF_fav$introg24 <- ABBA4DF_fav$introg24 + ABBA4DF2_fav$introg24
  ABBA4DF_fav$introg31 <- ABBA4DF_fav$introg31 + ABBA4DF2_fav$introg31
  ABBA4DF_fav$introg41 <- ABBA4DF_fav$introg41 + ABBA4DF2_fav$introg41
  ABBA4DF_fav$introg32 <- ABBA4DF_fav$introg32 + ABBA4DF2_fav$introg32
  ABBA4DF_fav$introg42 <- ABBA4DF_fav$introg42 + ABBA4DF2_fav$introg42
  ABBA4DF_fav$introg123 <- ABBA4DF_fav$introg123 + ABBA4DF2_fav$introg123
  ABBA4DF_fav$introg124 <- ABBA4DF_fav$introg124 + ABBA4DF2_fav$introg124
  ABBA4DF_fav$intrognone <- ABBA4DF_fav$intrognone + ABBA4DF2_fav$intrognone
}

DFOIL_cyr <- list.files(path='DFOIL/him_cyr', pattern='DFOIL', full.names = T)
DFOIL_cyr <- DFOIL_cyr[!grepl('BC2639', DFOIL_cyr)]
DF_cyr <- read.table(DFOIL_cyr[1], h = T, comment.char = "")

DF_cyr$X.chrom <- gsub('/home/rpapa/sbelleghem/work/DFOIL_out/him_cyr/', '', DF_cyr$X.chrom)
DF_cyr <- separate(DF_cyr, X.chrom, into = c('wind','indiv'), sep="_himera")
DF_cyr <- separate(DF_cyr, wind, into = c('scaffold','startDF','endDF'), sep="_")

DF_cyr <- DF_cyr[order(DF_cyr$scaffold, as.integer(DF_cyr$start)),  ]
DF_cyr$windID <-  paste(DF_cyr$scaffold,"_",DF_cyr$startDF, "_", DF_cyr$endDF, sep="")

ABBA1DF_cyr <- merge(ABBA1_cyr, DF_cyr, by = 'windID')
ABBA2DF_cyr <- merge(ABBA2_cyr, DF_cyr, by = 'windID')
ABBA3DF_cyr <- merge(ABBA3_cyr, DF_cyr, by = 'windID')
ABBA4DF_cyr <- merge(ABBA4_cyr, DF_cyr, by = 'windID')

for(e in 2:length(DFOIL_cyr)){
  DF_cyr <- read.table(DFOIL_cyr[e], h = T, comment.char = "")
  
  DF_cyr$X.chrom <- gsub('/home/rpapa/sbelleghem/work/DFOIL_out/him_cyr/', '', DF_cyr$X.chrom)
  DF_cyr <- separate(DF_cyr, X.chrom, into = c('wind','indiv'), sep="_himera")
  DF_cyr <- separate(DF_cyr, wind, into = c('scaffold','startDF','endDF'), sep="_")
  
  DF_cyr <- DF_cyr[order(DF_cyr$scaffold, as.integer(DF_cyr$start)),  ]
  DF_cyr$windID <-  paste(DF_cyr$scaffold,"_",DF_cyr$startDF, "_", DF_cyr$endDF, sep="")
  
  ABBA1DF2_cyr <- merge(ABBA1_cyr, DF_cyr, by = 'windID')
  ABBA2DF2_cyr <- merge(ABBA2_cyr, DF_cyr, by = 'windID')
  ABBA3DF2_cyr <- merge(ABBA3_cyr, DF_cyr, by = 'windID')
  ABBA4DF2_cyr <- merge(ABBA4_cyr, DF_cyr, by = 'windID')
  
  ABBA1DF_cyr$introg13 <- ABBA1DF_cyr$introg13 + ABBA1DF2_cyr$introg13
  ABBA1DF_cyr$introg14 <- ABBA1DF_cyr$introg14 + ABBA1DF2_cyr$introg14
  ABBA1DF_cyr$introg23 <- ABBA1DF_cyr$introg23 + ABBA1DF2_cyr$introg23
  ABBA1DF_cyr$introg24 <- ABBA1DF_cyr$introg24 + ABBA1DF2_cyr$introg24
  ABBA1DF_cyr$introg31 <- ABBA1DF_cyr$introg31 + ABBA1DF2_cyr$introg31
  ABBA1DF_cyr$introg41 <- ABBA1DF_cyr$introg41 + ABBA1DF2_cyr$introg41
  ABBA1DF_cyr$introg32 <- ABBA1DF_cyr$introg32 + ABBA1DF2_cyr$introg32
  ABBA1DF_cyr$introg42 <- ABBA1DF_cyr$introg42 + ABBA1DF2_cyr$introg42
  ABBA1DF_cyr$introg123 <- ABBA1DF_cyr$introg123 + ABBA1DF2_cyr$introg123
  ABBA1DF_cyr$introg124 <- ABBA1DF_cyr$introg124 + ABBA1DF2_cyr$introg124
  ABBA1DF_cyr$intrognone <- ABBA1DF_cyr$intrognone + ABBA1DF2_cyr$intrognone
  
  ABBA2DF_cyr$introg13 <- ABBA2DF_cyr$introg13 + ABBA2DF2_cyr$introg13
  ABBA2DF_cyr$introg14 <- ABBA2DF_cyr$introg14 + ABBA2DF2_cyr$introg14
  ABBA2DF_cyr$introg23 <- ABBA2DF_cyr$introg23 + ABBA2DF2_cyr$introg23
  ABBA2DF_cyr$introg24 <- ABBA2DF_cyr$introg24 + ABBA2DF2_cyr$introg24
  ABBA2DF_cyr$introg31 <- ABBA2DF_cyr$introg31 + ABBA2DF2_cyr$introg31
  ABBA2DF_cyr$introg41 <- ABBA2DF_cyr$introg41 + ABBA2DF2_cyr$introg41
  ABBA2DF_cyr$introg32 <- ABBA2DF_cyr$introg32 + ABBA2DF2_cyr$introg32
  ABBA2DF_cyr$introg42 <- ABBA2DF_cyr$introg42 + ABBA2DF2_cyr$introg42
  ABBA2DF_cyr$introg123 <- ABBA2DF_cyr$introg123 + ABBA2DF2_cyr$introg123
  ABBA2DF_cyr$introg124 <- ABBA2DF_cyr$introg124 + ABBA2DF2_cyr$introg124
  ABBA2DF_cyr$intrognone <- ABBA2DF_cyr$intrognone + ABBA2DF2_cyr$intrognone
  
  ABBA3DF_cyr$introg13 <- ABBA3DF_cyr$introg13 + ABBA3DF2_cyr$introg13
  ABBA3DF_cyr$introg14 <- ABBA3DF_cyr$introg14 + ABBA3DF2_cyr$introg14
  ABBA3DF_cyr$introg23 <- ABBA3DF_cyr$introg23 + ABBA3DF2_cyr$introg23
  ABBA3DF_cyr$introg24 <- ABBA3DF_cyr$introg24 + ABBA3DF2_cyr$introg24
  ABBA3DF_cyr$introg31 <- ABBA3DF_cyr$introg31 + ABBA3DF2_cyr$introg31
  ABBA3DF_cyr$introg41 <- ABBA3DF_cyr$introg41 + ABBA3DF2_cyr$introg41
  ABBA3DF_cyr$introg32 <- ABBA3DF_cyr$introg32 + ABBA3DF2_cyr$introg32
  ABBA3DF_cyr$introg42 <- ABBA3DF_cyr$introg42 + ABBA3DF2_cyr$introg42
  ABBA3DF_cyr$introg123 <- ABBA3DF_cyr$introg123 + ABBA3DF2_cyr$introg123
  ABBA3DF_cyr$introg124 <- ABBA3DF_cyr$introg124 + ABBA3DF2_cyr$introg124
  ABBA3DF_cyr$intrognone <- ABBA3DF_cyr$intrognone + ABBA3DF2_cyr$intrognone

  ABBA4DF_cyr$introg13 <- ABBA4DF_cyr$introg13 + ABBA4DF2_cyr$introg13
  ABBA4DF_cyr$introg14 <- ABBA4DF_cyr$introg14 + ABBA4DF2_cyr$introg14
  ABBA4DF_cyr$introg23 <- ABBA4DF_cyr$introg23 + ABBA4DF2_cyr$introg23
  ABBA4DF_cyr$introg24 <- ABBA4DF_cyr$introg24 + ABBA4DF2_cyr$introg24
  ABBA4DF_cyr$introg31 <- ABBA4DF_cyr$introg31 + ABBA4DF2_cyr$introg31
  ABBA4DF_cyr$introg41 <- ABBA4DF_cyr$introg41 + ABBA4DF2_cyr$introg41
  ABBA4DF_cyr$introg32 <- ABBA4DF_cyr$introg32 + ABBA4DF2_cyr$introg32
  ABBA4DF_cyr$introg42 <- ABBA4DF_cyr$introg42 + ABBA4DF2_cyr$introg42
  ABBA4DF_cyr$introg123 <- ABBA4DF_cyr$introg123 + ABBA4DF2_cyr$introg123
  ABBA4DF_cyr$introg124 <- ABBA4DF_cyr$introg124 + ABBA4DF2_cyr$introg124
  ABBA4DF_cyr$intrognone <- ABBA4DF_cyr$intrognone + ABBA4DF2_cyr$intrognone
}



ABBAs <- list(ABBA1DF_emm,ABBA2DF_emm,ABBA3DF_emm,ABBA4DF_emm)
ABBAs <- list(ABBA1DF_fav,ABBA2DF_fav,ABBA3DF_fav,ABBA4DF_fav)
ABBAs <- list(ABBA1DF_cyr,ABBA2DF_cyr,ABBA3DF_cyr,ABBA4DF_cyr)

ABBAs <- list(ABBA1DF_cyr,ABBA3DF_emm,ABBA3DF_fav)

# ABBAs <- list(ABBA1DF,ABBA2DF)
for(t in 1:length(ABBAs)){
  ABBAs[[t]]$I12_34 <- ABBAs[[t]]$introg13 + ABBAs[[t]]$introg14 + ABBAs[[t]]$introg23 + ABBAs[[t]]$introg24
  ABBAs[[t]]$I34_12 <- ABBAs[[t]]$introg31 + ABBAs[[t]]$introg32 + ABBAs[[t]]$introg41 + ABBAs[[t]]$introg42
  ABBAs[[t]]$I12_A34 <- ABBAs[[t]]$introg123 + ABBAs[[t]]$introg124 
  ABBAs[[t]]$sum <- 10#ABBAs[[t]]$I12_34 + ABBAs[[t]]$I34_12 #+ ABBAs[[t]]$I12_A34 #+ ABBAs[[t]]$intrognone
  
  
  colfuncR <- colorRampPalette(c("white", "red"))
  colfuncG <- colorRampPalette(c("white", "blue"))
  colfuncB <- colorRampPalette(c("white", "blue"))
  
  
  ABBAs[[t]]$col <- "white"
  ABBAs[[t]]$val <- 1
  
  multiplyFactor <- c(0.6, 0.96, 1)
  
  for(e in 1:nrow(ABBAs[[t]])){
    # if(ABBAs[[t]]$I12_34[e] == 0 & ABBAs[[t]]$I34_12[e] == 0 & ABBAs[[t]]$I12_A34[e] == 0){
    #   ABBAs[[t]]$col[e] <- "gray"
    #   ABBAs[[t]]$val[e] <- adjustcolor(ABBAs[[t]]$col[e], alpha.f = 1)
    # }
    if(ABBAs[[t]]$I12_34[e] > ABBAs[[t]]$I34_12[e] ){
      ABBAs[[t]]$col[e] <- "red"
      colR <- colfuncR(100)
      ABBAs[[t]]$val[e] <- colR[((log(ABBAs[[t]]$I12_34[e]+1)/log(max(ABBAs[[1]]$I12_34,ABBAs[[1]]$I34_12,ABBAs[[2]]$I12_34,ABBAs[[2]]$I34_12,ABBAs[[3]]$I12_34,ABBAs[[3]]$I34_12))*100))*multiplyFactor[t]] #adjustcolor(ABBAs[[t]]$col[e], alpha.f = ABBAs[[t]]$I12_34[e]/sampleCombs)
    }
    if(ABBAs[[t]]$I12_34[e] < ABBAs[[t]]$I34_12[e] ){
      ABBAs[[t]]$col[e] <- "blue"
      colG <- colfuncG(100)
      ABBAs[[t]]$val[e] <- colG[((log(ABBAs[[t]]$I34_12[e]+1)/log(max(ABBAs[[1]]$I12_34,ABBAs[[1]]$I34_12,ABBAs[[2]]$I12_34,ABBAs[[2]]$I34_12,ABBAs[[3]]$I12_34,ABBAs[[3]]$I34_12))*100))*multiplyFactor[t]] #adjustcolor(ABBAs[[t]]$col[e], alpha.f = ABBAs[[t]]$I34_12[e]/sampleCombs)
    }
    # if(ABBAs[[t]]$I12_A34[e] > ABBAs[[t]]$I12_34[e] & ABBAs[[t]]$I12_A34[e] > ABBAs[[t]]$I34_12[e]){
    #   ABBAs[[t]]$col[e] <- "blue"
    #   colB <- colfuncB(800)
    #   # ABBAs[[t]]$val[e] <- colB[(log(ABBAs[[t]]$I12_A34[e]+1)/log(max(ABBAs[[t]]$I12_34,ABBAs[[t]]$I34_12,ABBAs[[t]]$I12_A34))*100)] #adjustcolor(ABBAs[[t]]$col[e], alpha.f = ABBAs[[t]]$I12_A34[e]/sampleCombs)
    #   ABBAs[[t]]$val[e] <- colB[ABBAs[[t]]$I12_A34[e]] #adjustcolor(ABBAs[[t]]$col[e], alpha.f = ABBAs[[t]]$I12_A34[e]/sampleCombs)
    # }
  }
  
  # ABBAs[[t]] <- ABBAs[[t]][order(ABBAs[[t]]$I12_34, ABBAs[[t]]$I34_12),]
  ABBAs[[t]] <- ABBAs[[t]][order(ABBAs[[t]]$scaffold.x,ABBAs[[t]]$start),]
  
  }

####
n <- 5
quant1 <- ABBAs[[1]][ABBAs[[1]]$fd > quantile(ABBAs[[1]]$fd,prob=1-n/100),]

nrow(quant1)
nrow(subset(quant1, quant1$col == 'white'))
nrow(subset(quant1, quant1$col == 'darkorange'))
nrow(subset(quant1, quant1$col == 'blue'))

quant2 <- ABBAs[[2]][ABBAs[[2]]$fd > quantile(ABBAs[[2]]$fd,prob=1-n/100),]

nrow(quant2)
nrow(subset(quant2, quant2$col == 'white'))
nrow(subset(quant2, quant2$col == 'darkorange'))
nrow(subset(quant2, quant2$col == 'blue'))

quant3 <- ABBAs[[3]][ABBAs[[3]]$fd > quantile(ABBAs[[3]]$fd,prob=1-n/100),]

nrow(quant3)
nrow(subset(quant3, quant3$col == 'white'))
nrow(subset(quant3, quant3$col == 'darkorange'))
nrow(subset(quant3, quant3$col == 'blue'))

boxplot(quant1$I12_34, quant1$I34_12, col = c("darkorange","blue"), ylim=c(0,60))
boxplot(quant2$I12_34, quant2$I34_12, col = c("darkorange","blue"), ylim=c(0,60))
boxplot(quant3$I12_34, quant3$I34_12, col = c("darkorange","blue"), ylim=c(0,60))

boxplot(ABBAs[[1]]$I12_34, ABBAs[[1]]$I34_12, col = c("darkorange","blue"), ylim=c(0,10))
boxplot(ABBAs[[2]]$I12_34, ABBAs[[2]]$I34_12, col = c("darkorange","blue"), ylim=c(0,10))
boxplot(ABBAs[[3]]$I12_34, ABBAs[[3]]$I34_12, col = c("darkorange","blue"), ylim=c(0,10))

boxplot(ABBAs[[1]]$I12_A34 ~ ABBAs[[1]]$chromosome, ylim=c(0,200))
boxplot(ABBAs[[2]]$I12_A34 ~ ABBAs[[2]]$chromosome, ylim=c(0,200))
boxplot(ABBAs[[3]]$I12_A34 ~ ABBAs[[3]]$chromosome, ylim=c(0,200))

boxplot(ABBAs[[1]]$I12_34+ABBAs[[1]]$I34_12 ~ ABBAs[[1]]$chromosome, ylim=c(0,20))
boxplot(ABBAs[[2]]$I12_34+ABBAs[[2]]$I34_12 ~ ABBAs[[2]]$chromosome, ylim=c(0,20))
boxplot(ABBAs[[3]]$I12_34+ABBAs[[3]]$I34_12 ~ ABBAs[[3]]$chromosome, ylim=c(0,20))

boxplot(quant1$I34_12 ~ quant1$chromosome, ylim=c(0,10))
boxplot(quant2$I34_12 ~ quant2$chromosome, ylim=c(0,10))
boxplot(quant3$I34_12 ~ quant3$chromosome, ylim=c(0,10))


ABBA_corr <- ABBAs[[1]]
ABBA_corr <- subset(ABBA_corr,ABBA_corr$fd >=0)

plot(ABBA_corr$fd, -ABBA_corr$I34_12/800, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="black")
par(new=TRUE)
plot(quant1$fd, -quant1$I34_12/800, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="blue")
par(new=TRUE)
plot(ABBA_corr$fd, ABBA_corr$I12_34/800, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="black")
par(new=TRUE)
plot(quant1$fd, quant1$I12_34/800, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="red")

lo<-loess.smooth(ABBA_corr$fd, ABBA_corr$I12_34/800, span = 1/5, degree = 1,
                 family = c("symmetric", "gaussian"), evaluation = 50)
lines(lo, col="green", lwd=2)
# lo<-loess.smooth(ABBAs[[1]]$fd, -ABBAs[[1]]$I34_12, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col="pink")

ABBA_corr2 <- ABBAs[[2]]
ABBA_corr2 <- subset(ABBA_corr2,ABBA_corr2$fd >=0)

plot(ABBA_corr2$fd, -ABBA_corr2$I34_12/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="black")
par(new=TRUE)
plot(quant2$fd, -quant2$I34_12/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="blue")
par(new=TRUE)
plot(ABBA_corr2$fd, ABBA_corr2$I12_34/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="black")
par(new=TRUE)
plot(quant2$fd, quant2$I12_34/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="red")

lo<-loess.smooth(ABBA_corr2$fd, ABBA_corr2$I12_34/500, span = 1/5, degree = 1,
                 family = c("symmetric", "gaussian"), evaluation = 50)
lines(lo, col="green", lwd=2)
# lo<-loess.smooth(ABBAs[[2]]$fd, -ABBAs[[2]]$I34_12/500, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col="pink")

ABBA_corr3 <- ABBAs[[3]]
ABBA_corr3 <- subset(ABBA_corr3,ABBA_corr3$fd >=0)

plot(ABBA_corr3$fd, -ABBA_corr3$I34_12/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="black")
par(new=TRUE)
plot(quant3$fd, -quant3$I34_12/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="blue")
par(new=TRUE)
plot(ABBA_corr3$fd, ABBA_corr3$I12_34/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="black")
par(new=TRUE)
plot(quant3$fd, quant3$I12_34/500, pch=19, cex=0.5, ylim=c(-0.2,0.2), xlim=c(0,0.4), col="red")

lo<-loess.smooth(ABBA_corr3$fd, ABBA_corr3$I12_34/500, span = 1/5, degree = 1,
                 family = c("symmetric", "gaussian"), evaluation = 50)
lines(lo, col="green", lwd=2)
# lo<-loess.smooth(ABBAs[[3]]$fd, -ABBAs[[3]]$I34_12, span = 1/5, degree = 1,
#                  family = c("symmetric", "gaussian"), evaluation = 50)
# lines(lo, col="pink")



par(new=TRUE)
plot(ABBAs[[1]]$fd, ABBAs[[1]]$I34_12, pch=19, cex=0.8, ylim=c(0,60), xlim=c(-0.4,0.4), col="blue")

# boxplot(log(ABBA1DF_cyr_5p$I12_A34+1)/log(800)*100, log(ABBA1DF_cyr_5p$I12_34+1)/log(800)*100, log(ABBA1DF_cyr_5p$I34_12)/log(800)*100)

#####

#####


comps <- c('((phy,himN),emmW),her','((phy,himN),emmE),her','((phy,himS),emmW),her','((phy,himS),emmE),her')
# comps <- c('((phy,himN),favW),her','((phy,himN),favE),her','((phy,himS),favW),her','((phy,himS),favE),her')
# comps <- c('((phy,himN),cyrS),her','((phy,himN),cyrN),her','((phy,himS),cyrS),her','((phy,himS),cyrN),her')

comps <- c('((himS,himN),cyrS),her','((himS,himN),cyrN),her','((himS,himN),cyrS),her','((himS,himN),cyrN),her')

par(mfrow=c(3,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)
top = 0.5
bot = 0

begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[21]/1000000

for(e in 1:length(ABBAs)){
  ABBAsub <- subset(ABBAs[[e]], ABBAs[[e]]$col == "blue" | ABBAs[[e]]$col == "red")
  
  ABBAsubA <- subset(ABBAs[[e]], ABBAs[[e]]$col == "white")
  ABBAsubB <- subset(ABBAs[[e]], ABBAs[[e]]$col == "blue")  
  ABBAsubC <- subset(ABBAs[[e]], ABBAs[[e]]$col == "blue")
  ABBAsubD <- subset(ABBAs[[e]], ABBAs[[e]]$col == "red")
  
  ABBAsubA <- subset(ABBAsubA, ABBAsubA$fd >=0)
  ABBAsubB <- subset(ABBAsubB, ABBAsubB$fd >=0)
  ABBAsubC <- subset(ABBAsubC, ABBAsubC$fd >=0)
  ABBAsubD <- subset(ABBAsubD, ABBAsubD$fd >=0)
  ABBAsub <- subset(ABBAsub, ABBAsub$fd >=0)

  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray45","gray30"), lwd = 0, border = c("gray65","gray50"))

  par(new=TRUE)
  plot(ABBAsubA$chromPos/1000000,ABBAsubA$fd, type="p",pch=16, cex=0.8,col=ABBAsubA$col,xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")

  par(new=TRUE)
  plot(ABBAsub$chromPos/1000000,ABBAsub$fd, type="p",pch=16, cex=0.8,col=ABBAsub$val,xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  
  x <- ABBAs[[e]]$chromPos
  y <- ABBAs[[e]]$fd
  lo<-loess(y ~x, span = 1/200, degree = 2, family = c( "gaussian"))
  smoothed <- predict(lo)

  lines(smoothed, x = ABBAs[[e]]$chromPos/1000000, col='green')

  # par(new=TRUE)
  # plot(ABBA1DF_cyr_5p$chromPos/1000000,ABBA1DF_cyr_5p$D, type="p",pch=16, cex=0.8,col="purple",xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  # 
  # ABBA1DF_cyr_5p
  # par(new=TRUE)
  # plot(ABBAsubB$chromPos/1000000,ABBAsubB$fd, type="p",pch=16, cex=0.8,col=ABBAsubB$val,xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")

  #
  # par(new=TRUE)
  # plot(ABBAsubD$chromPos/1000000,ABBAsubD$fd, type="p",pch=16, cex=0.8,col=ABBAsubD$val,xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  # par(new=TRUE)
  # plot(ABBAsubC$chromPos/1000000,ABBAsubC$fd, type="p",pch=16, cex=0.8,col=ABBAsubC$val,xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")

  axis(2,cex.axis = 1, line = -1.5)
  
  text(50,0.8, comps[e])
  mtext(side = 2, text = expression(paste("fd")), cex=1,line = 0.8)
}

axis(1, at=chrom_coords[,5][1:21]/1000000, labels=(1:21),lwd=0, lwd.ticks=0)
segments(c(390,390,400), c(250.07,250.0702,250.0702), c(400,390,400) ,c(250.07,250.0702,250.0702), lwd = 1)
text(395,300.086,labels = "10Mb", cex = 1)


###
# length admixture tracts

tableInter <- list()
for(i in 1:3){
  compar <- ABBAs[[i]]
  
  x <- compar$chromPos
  y <- compar$fd
  lo<-loess(y ~x, span = 1/200, degree = 2, family = c( "gaussian"))
  smoothed <- predict(lo)
  
  # lines(smoothed, x = compar$chromPos/1000000, col='green')
  
  lo_out <- predict(lo)
  
  
  e <- 1
  intervals_lo <- c()
  while(e <= nrow(compar)){
    if(lo_out[e] >= 0.05){
      start <- compar$chromPos[e]
      chromos <- compar$chromosome[e]
      
      while(lo_out[e+1] >=0.05){
        end <- compar$chromPos[e+1]
        e = e +1
      }
      
      intervals_lo <- rbind(intervals_lo, c(as.integer(chromos), start, end, end-start))
    }
    e = e +1
  }
  
  intervals_lo <- as.data.frame(intervals_lo)
  intervals_lo$V4[intervals_lo$V4<0] <- 50000
  
  colnames(intervals_lo) <- c('chrom', 'start', 'end', 'length')
  
  tableInter[[i]] <- intervals_lo
  
  print(c(sum(intervals_lo$length), mean(intervals_lo$length), sd(intervals_lo$length)))
}


cyr <- tableInter[[1]][c('start','end')]
emm <- tableInter[[2]][c('start','end')]
fav <- tableInter[[3]][c('start','end')]

nrow(cyr)
nrow(emm)
nrow(fav)

library(intervals)

for(n in 1:nrow(cyr)){
  if(cyr$end[n] < cyr$start[n]){
    cyr$end[n] <- cyr$start[n]+50000
  }
}

for(n in 1:nrow(emm)){
  if(emm$end[n] < emm$start[n]){
    emm$end[n] <- emm$start[n]+50000
  }
}

for(n in 1:nrow(fav)){
  if(fav$end[n] < fav$start[n]){
    fav$end[n] <- fav$start[n]+50000
  }
}
cyrI <- Intervals(as.matrix(cyr))
emmI <- Intervals(as.matrix(emm))
favI <- Intervals(as.matrix(fav))

cyremm <- as.data.frame(interval_intersection(cyrI,emmI))
cyrfav <- as.data.frame(interval_intersection(cyrI,favI))
emmfav <- as.data.frame(interval_intersection(emmI,favI))

cyremmD <- sum(cyremm$V2-cyremm$V1)
cyrfavD <- sum(cyrfav$V2-cyrfav$V1)
emmfavD <- sum(emmfav$V2-emmfav$V1)


###








### block jack knife

jack <- function(stat, w){
  
  start <- 0
  end <- w
  
  n <- 1
  
  D <- sum(stat)/length(stat)
  
  len <- as.integer(length(stat)/w)
  pseudo_stat <- vector(length = len)
  
  while(end <= length(stat)){
    
    jacked_stat <- sum(stat[-c(start:end)])/length(stat[-c(start:end)])
    
    pseudo_stat[n] <- D*len - jacked_stat*(len-1)
    
    start <- start + w
    end <- start + w
    n <- n + 1
  }
  return(pseudo_stat)
}

jacks <- jack(ABBAs[[3]]$D, 20)
jacksJ <- sum(jacks)/length(jacks)
jackssd <- sd(jacks)
jackserr <- sqrt(var(jacks)/length(jacks))

D_Z <- mean(ABBAs[[3]]$fd)/jackserr
D_p <- 2*pnorm(-abs(D_Z))

1.96*jackserr

abline(1.96*jackserr, 0)
