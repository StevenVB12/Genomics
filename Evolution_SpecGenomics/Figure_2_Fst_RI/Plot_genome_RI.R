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



LAT_NOT <- read.table("Fst_stats/lat_notLAT.stats",h=T)
ETY_NOT <- read.table("Fst_stats/ety_notETY.stats",h=T)
ERA_HYDfg <- read.table("Fst_stats/era_hydFG.stats",h=T)
EMM_FAVem <- read.table("Fst_stats/emmFAV_favEMM.stats",h=T)
DEM_HYDp <- read.table("Fst_stats/dem_hydP.stats",h=T)
AMA_ERA <- read.table("Fst_stats/ama_era.stats",h=T)
AMA_HYDfg <- read.table("Fst_stats/ama_hydFG.stats",h=T)
HIM_CYR <- read.table("Fst_stats/himCYR_cyr.stats",h=T)
EMM_HIM <- read.table("Fst_stats/emmHIM_himEM.stats",h=T)
FAV_HIM <- read.table("Fst_stats/favEMM-himEM.stats",h=T)
VEN_CHE <- read.table("Fst_stats/ven_che.stats",h=T)
PHY_FAVem <- read.table("Fst_stats/favEMM_phy.stats",h=T)
CYR_VEN <- read.table("Fst_stats/cyr_ven.stats",h=T)
CYRN_VEN <- read.table("Fst_stats/cyrN-ven.stats", h=T)
ETY_HIMem <- read.table("Fst_stats/ety-himEM.stats",h=T)
cyrN_emmFAV <- read.table("Fst_stats/cyr-emmHIM.stats",h=T)
head(CYR_VEN)

LAT_NOT <- merge(LAT_NOT,scaf_coords2,by="scaffold", all.x=TRUE)
ETY_NOT <- merge(ETY_NOT,scaf_coords2,by="scaffold", all.x=TRUE)
ERA_HYDfg <- merge(ERA_HYDfg,scaf_coords2,by="scaffold", all.x=TRUE)
EMM_FAVem <- merge(EMM_FAVem,scaf_coords2,by="scaffold", all.x=TRUE)
DEM_HYDp <- merge(DEM_HYDp,scaf_coords2,by="scaffold", all.x=TRUE)
AMA_ERA <- merge(AMA_ERA,scaf_coords2,by="scaffold", all.x=TRUE)
AMA_HYDfg <- merge(AMA_HYDfg,scaf_coords2,by="scaffold", all.x=TRUE)
HIM_CYR <- merge(HIM_CYR,scaf_coords2,by="scaffold", all.x=TRUE)
EMM_HIM <- merge(EMM_HIM,scaf_coords2,by="scaffold", all.x=TRUE)
FAV_HIM <- merge(FAV_HIM,scaf_coords2,by="scaffold", all.x=TRUE)
VEN_CHE <- merge(VEN_CHE,scaf_coords2,by="scaffold", all.x=TRUE)
PHY_FAVem <- merge(PHY_FAVem,scaf_coords2,by="scaffold", all.x=TRUE)
CYR_VEN <- merge(CYR_VEN,scaf_coords2,by="scaffold", all.x=TRUE)
CYRN_VEN <- merge(CYRN_VEN,scaf_coords2,by="scaffold", all.x=TRUE)
ETY_HIMem <- merge(ETY_HIMem,scaf_coords2,by="scaffold", all.x=TRUE)
cyrN_emmFAV <- merge(cyrN_emmFAV,scaf_coords2,by="scaffold", all.x=TRUE)


LAT_NOT <- cbind(LAT_NOT,chromPos=LAT_NOT$position+LAT_NOT$scafStart+LAT_NOT$chromStart-2)
ETY_NOT <- cbind(ETY_NOT,chromPos=ETY_NOT$position+ETY_NOT$scafStart+ETY_NOT$chromStart-2)
ERA_HYDfg <- cbind(ERA_HYDfg,chromPos=ERA_HYDfg$position+ERA_HYDfg$scafStart+ERA_HYDfg$chromStart-2)
EMM_FAVem <- cbind(EMM_FAVem,chromPos=EMM_FAVem$position+EMM_FAVem$scafStart+EMM_FAVem$chromStart-2)
DEM_HYDp <- cbind(DEM_HYDp,chromPos=DEM_HYDp$position+DEM_HYDp$scafStart+DEM_HYDp$chromStart-2)
AMA_ERA <- cbind(AMA_ERA,chromPos=AMA_ERA$position+AMA_ERA$scafStart+AMA_ERA$chromStart-2)
AMA_HYDfg <- cbind(AMA_HYDfg,chromPos=AMA_HYDfg$position+AMA_HYDfg$scafStart+AMA_HYDfg$chromStart-2)
HIM_CYR <- cbind(HIM_CYR,chromPos=HIM_CYR$position+HIM_CYR$scafStart+HIM_CYR$chromStart-2)
EMM_HIM <- cbind(EMM_HIM,chromPos=EMM_HIM$position+EMM_HIM$scafStart+EMM_HIM$chromStart-2)
FAV_HIM <- cbind(FAV_HIM,chromPos=FAV_HIM$position+FAV_HIM$scafStart+FAV_HIM$chromStart-2)
VEN_CHE <- cbind(VEN_CHE,chromPos=VEN_CHE$position+VEN_CHE$scafStart+VEN_CHE$chromStart-2)
PHY_FAVem <- cbind(PHY_FAVem,chromPos=PHY_FAVem$position+PHY_FAVem$scafStart+PHY_FAVem$chromStart-2)
CYR_VEN <- cbind(CYR_VEN,chromPos=CYR_VEN$position+CYR_VEN$scafStart+CYR_VEN$chromStart-2)
CYRN_VEN <- cbind(CYRN_VEN,chromPos=CYRN_VEN$position+CYRN_VEN$scafStart+CYRN_VEN$chromStart-2)
ETY_HIMem <- cbind(ETY_HIMem,chromPos=ETY_HIMem$position+ETY_HIMem$scafStart+ETY_HIMem$chromStart-2)
cyrN_emmFAV <- cbind(cyrN_emmFAV,chromPos=cyrN_emmFAV$position+cyrN_emmFAV$scafStart+cyrN_emmFAV$chromStart-2)


comp=list(LAT_NOT,ETY_NOT,ERA_HYDfg,AMA_HYDfg,AMA_ERA,EMM_FAVem,DEM_HYDp,HIM_CYR,EMM_HIM,FAV_HIM,VEN_CHE,PHY_FAVem,CYR_VEN,CYRN_VEN,ETY_HIMem,cyrN_emmFAV)
names=c("LAT_NOT","ETY_NOT","ERA_HYDfg","AMA_HYDfg","AMA_ERA","EMM_FAV","DEM_HYDp","HIM_CYR","EMM_HIM","FAV_HIM","VEN_CHE","PHY_FAVem","CYR_VEN","CYRN_VEN","ETY_HIMem"," cyrN_emmFAV")
col=c("purple","blue","orange","black","brown","red","green","darkgreen","mediumorchid4","mediumorchid3","mediumorchid2","cadetblue","cyan4","cyan4","cyan4","cyan4")


sum(LAT_NOT$Fst)/nrow(LAT_NOT)
sum(ETY_NOT$Fst)/nrow(ETY_NOT)
sum(ERA_HYDfg$Fst)/nrow(ERA_HYDfg)
sum(EMM_FAVem$Fst)/nrow(EMM_FAVem)
sum(DEM_HYDp$Fst)/nrow(DEM_HYDp)
sum(AMA_ERA$Fst)/nrow(AMA_ERA)
sum(AMA_HYDfg$Fst)/nrow(AMA_HYDfg)
sum(HIM_CYR$Fst)/nrow(HIM_CYR)
sum(EMM_HIM$Fst)/nrow(EMM_HIM)
sum(FAV_HIM$Fst)/nrow(FAV_HIM)
sum(VEN_CHE$Fst)/nrow(VEN_CHE)
sum(PHY_FAVem$Fst)/nrow(PHY_FAVem)
sum(CYR_VEN$Fst)/nrow(CYR_VEN)


comp=list(DEM_HYDp,ERA_HYDfg,ETY_NOT,LAT_NOT,EMM_FAVem,FAV_HIM,EMM_HIM,VEN_CHE,HIM_CYR)
names=c("DEM_HYDp","ERA_HYDfg","ETY_NOT","LAT_NOT","EMM_FAVem","FAV_HIM","EMM_HIM","VEN_CHE","HIM_CYR")
col=c("green","orange","blue","purple","red","darkgreen","mediumvioletred","cyan4","black")

# comp=list(EMM_FAVem)
# names=c("EMM_FAVem")
# col=c("red")

comp=list(DEM_HYDp,EMM_FAVem,FAV_HIM,EMM_HIM,HIM_CYR,VEN_CHE)
names=c("DEM_HYDp","EMM_FAVem","FAV_HIM","EMM_HIM","HIM_CYR","VEN_CHE")
col=c("black","black","black","black","black","black")



library("plotrix")


par(mfrow=c(16,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 0.05
bot = 0

begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[21]/1000000

for (i in 7-1:length(comp)){
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))
  #gradient.rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromEnds)), col = smoothColors("gray95",38,"gray90"), border = c("gray95","gray90"))
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$DxyNoN-(comp[[i]]$PiPop2+comp[[i]]$PiPop2)/2, type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=1),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("da")), cex=0.5,line = 0.8)
  mtext(side = 4, text = names[i], cex=0.5)
  print(mean(comp[[i]]$DxyNoN-(comp[[i]]$PiPop1+comp[[i]]$PiPop2)/2))
}

# x<-subset(AMA_HYD,AMA_HYD$Fst>0.3)
# subset(x,x$scaffold=='Herato0310')
axis(1, at=chrom_coords[,5][1:21]/1000000, labels=(1:21),lwd=0, lwd.ticks=0)
segments(c(390,390,400), c(250.07,250.0702,250.0702), c(400,390,400) ,c(250.07,250.0702,250.0702), lwd = 1)
text(395,300.086,labels = "10Mb", cex = 1)




par(mfrow=c(16,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 1
bot = 0

begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[21]/1000000

for (i in 1:length(comp)){
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))
  #gradient.rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromEnds)), col = smoothColors("gray95",38,"gray90"), border = c("gray95","gray90"))
  par(new=TRUE)
  tableS <- subset(comp[[i]], comp[[i]]$Fst >=0)
  plot(tableS$chromPos/1000000,tableS$Fst, type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=1),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Fst")), cex=0.5,line = 0.8)
  mtext(side = 4, text = names[i], cex=0.5)
  print(mean(comp[[i]]$DxyNoN-(comp[[i]]$PiPop1+comp[[i]]$PiPop2)/2))
}

# x<-subset(AMA_HYD,AMA_HYD$Fst>0.3)
# subset(x,x$scaffold=='Herato0310')
axis(1, at=chrom_coords[,5][1:21]/1000000, labels=(1:21),lwd=0, lwd.ticks=0)
segments(c(390,390,400), c(250.07,250.0702,250.0702), c(400,390,400) ,c(250.07,250.0702,250.0702), lwd = 1)
text(395,300.086,labels = "10Mb", cex = 1)







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

# block jackknife

comp=list(DEM_HYDp,ERA_HYDfg,AMA_ERA,AMA_HYDfg,ETY_NOT,LAT_NOT,EMM_FAVem,FAV_HIM,EMM_HIM,VEN_CHE,HIM_CYR,CYRN_VEN,ETY_HIMem)
names=c("DEM_HYDp","ERA_HYDfg","AMA_ERA","AMA_HYDfg","ETY_NOT","LAT_NOT","EMM_FAVem","FAV_HIM","EMM_HIM","VEN_CHE","HIM_CYR","CyrN_VEN","ETY_HIMem")

w <- 100

blockJ <-c()
for(e in 1:length(comp)){
  
auto <- subset(comp[[e]], comp[[e]]$scaffold != 'Herato2101')
sex <- subset(comp[[e]],comp[[e]]$scaffold == 'Herato2101')

FstAjacks <- jack(auto$Fst, w)
FstAJ <- sum(FstAjacks)/length(FstAjacks)
FstAJsd <- sd(FstAjacks)
FstAJerr <- sqrt(var(FstAjacks)/length(FstAjacks))

dXYjacks <- jack(auto$DxyNoN, w)
dXYJ <- sum(dXYjacks)/length(dXYjacks)
dXYJsd <- sd(dXYjacks)
dXYJerr <- sqrt(var(dXYjacks)/length(dXYjacks))

dajacks <- jack(auto$DxyNoN-(auto$PiPop2+auto$PiPop2)/2, w)
daJ <- sum(dajacks)/length(dajacks)
daJsd <- sd(dajacks)
daJerr <- sqrt(var(dajacks)/length(dajacks))


FstAjacksZ <- jack(sex$Fst, w)
FstAJZ <- sum(FstAjacksZ)/length(FstAjacksZ)
FstAJsdZ <- sd(FstAjacksZ)
FstAJerrZ <- sqrt(var(FstAjacksZ)/length(FstAjacksZ))

dXYjacksZ <- jack(sex$DxyNoN, w)
dXYJZ <- sum(dXYjacksZ)/length(dXYjacksZ)
dXYJsdZ <- sd(dXYjacksZ)
dXYJerrZ <- sqrt(var(dXYjacksZ)/length(dXYjacksZ))

dajacksZ <- jack(sex$DxyNoN-(sex$PiPop2+sex$PiPop2)/2, w)
daJZ <- sum(dajacksZ)/length(dajacksZ)
daJsdZ <- sd(dajacksZ)
daJerrZ <- sqrt(var(dajacksZ)/length(dajacksZ))

blockJ <- rbind(blockJ, c(names[e], FstAJ, FstAJerr,FstAJsd,dXYJ,dXYJerr,dXYJsd,daJ,daJerr,daJsd, FstAJZ, FstAJerrZ,FstAJsdZ,dXYJZ,dXYJerrZ,dXYJsdZ,daJZ,daJerrZ,daJsdZ))

}

blockJ <- as.data.frame(blockJ)

colnames(blockJ) <- c('pop','FstJ','FstJerr','Fstsd','dXYJ','dXYJerr','dXYsd','daJ','daJerr','dasd',
                      'FstJZ','FstJerrZ','FstsdZ','dXYJZ','dXYJerrZ','dXYsdZ','daJZ','daJerrZ','dasdZ')

blockJ <- subset(blockJ, blockJ$pop != 'CyrN_VEN')

# RI <- c(0.11, 0.22, 0.11, 0.22, 0.22, 0.22, 0.33, 0.5, 0.5, 0.89, 0.67, 0.44, 0.5)
RI <- c(0.0002, 0.51, 0.04, 0.04, 0.44, 0.44, 0.61, (0.672+0.804)/2, (0.5226+0.7041)/2, 0.9068, 0.804, (0.5226+0.7041)/2)

RI1 <- c(0.0002, 0.51, 0.04, 0.04, 0.44, 0.44, 0.61, (0.672), (0.5226), 0.9068, 0.804, (0.5226))
RI2 <- c(0.0002, 0.51, 0.04, 0.04, 0.44, 0.44, 0.61, (0.804), (0.7041), 0.9068, 0.804, (0.7041))

blockJ <- cbind(blockJ, RI, RI1, RI2)

par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(NULL, pch = 19, ylab = 'Fst', xlab = 'RI',bty='none', xlim=c(0,1), ylim=c(0,1), cex=3)


x <- as.numeric(as.character(blockJ$RI))
y <- as.numeric(as.character(blockJ$FstJ))
m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = 0, 
                                                                   b = 0.6, c = 10, d = 0.5), trace = TRUE)
lines(seq(0,1,length.out = 100), predict(m.s, newdata=data.frame(x=seq(0,1,length.out = 100))), lty = 2, lwd = 2, col = "gray")

x <- as.numeric(as.character(blockJ$RI))
y <- as.numeric(as.character(blockJ$FstJZ))
m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = 0, 
                                                                   b = 0.6, c = 10, d = 0.5), trace = TRUE)
lines(seq(0,1,length.out = 100), predict(m.s, newdata=data.frame(x=seq(0,1,length.out = 100))), lty = 2, lwd = 2, col = "gray")

points(blockJ$RI, as.numeric(as.character(blockJ$FstJ)), col='black', pch=19, cex=2)
points(blockJ$RI, as.numeric(as.character(blockJ$FstJZ)), col='pink', pch=21, cex=2)

points(blockJ$RI1, as.numeric(as.character(blockJ$FstJ)), col='black', pch=19, cex=1)
points(blockJ$RI1, as.numeric(as.character(blockJ$FstJZ)), col='pink', pch=21, cex=1)

points(blockJ$RI2, as.numeric(as.character(blockJ$FstJ)), col='black', pch=19, cex=1)
points(blockJ$RI2, as.numeric(as.character(blockJ$FstJZ)), col='pink', pch=21, cex=1)



par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(NULL, pch = 19, ylab = 'da', xlab = 'RI',bty='none', xlim=c(0,1), ylim=c(0,0.03), cex=3)


x <- as.numeric(as.character(blockJ$RI))
y <- as.numeric(as.character(blockJ$daJ))
m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = 0.001, 
                                                                   b = 0.018, c = 10, d = 0.7), trace = TRUE)
lines(seq(0,1,length.out = 100), predict(m.s, newdata=data.frame(x=seq(0,1,length.out = 100))), lty = 2, lwd = 2, col = "gray")

x <- as.numeric(as.character(blockJ$RI))
y <- as.numeric(as.character(blockJ$daJZ))
m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = 0, 
                                                                   b = 0.26, c = 10, d = 0.7), trace = TRUE)
lines(seq(0,1,length.out = 100), predict(m.s, newdata=data.frame(x=seq(0,1,length.out = 100))), lty = 2, lwd = 2, col = "gray")

points(blockJ$RI, as.numeric(as.character(blockJ$daJ)), col='black', pch=19, cex=2)
points(blockJ$RI, as.numeric(as.character(blockJ$daJZ)), col='pink', pch=21, cex=2)

points(blockJ$RI1, as.numeric(as.character(blockJ$daJ)), col='black', pch=19, cex=1)
points(blockJ$RI1, as.numeric(as.character(blockJ$daJZ)), col='pink', pch=21, cex=1)

points(blockJ$RI2, as.numeric(as.character(blockJ$daJ)), col='black', pch=19, cex=1)
points(blockJ$RI2, as.numeric(as.character(blockJ$daJZ)), col='pink', pch=21, cex=1)


text(as.numeric(as.character(blockJ$RI)), as.numeric(as.character(blockJ$daJ)), labels = blockJ$pop)
# a = 0
# b = 0.6
# c = 10
# d = 0.5
# curve((b-a)/(1 + exp(-c * (x-d))), 0, 1)
# 
# fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y))
# summary(fit)

text(blockJ$RI, as.numeric(as.character(blockJ$FstJ)), blockJ$pop)
segments(blockJ$RI,as.numeric(as.character(blockJ$FstJ))-as.numeric(as.character(blockJ$FstJerr)),blockJ$RI,as.numeric(as.character(blockJ$FstJ))+as.numeric(as.character(blockJ$FstJerr)))
plot(blockJ$RI, as.numeric(as.character(blockJ$dXYJ)), pch = 19, ylab = 'dXY', xlab = 'RI')
segments(blockJ$RI,as.numeric(as.character(blockJ$dXYJ))-as.numeric(as.character(blockJ$dXYJerr)),blockJ$RI,as.numeric(as.character(blockJ$dXYJ))+as.numeric(as.character(blockJ$dXYJerr)))
plot(blockJ$RI, as.numeric(as.character(blockJ$daJ)), pch = 19, ylab = 'da', xlab = 'RI')
segments(blockJ$RI,as.numeric(as.character(blockJ$daJ))-as.numeric(as.character(blockJ$daJerr)),blockJ$RI,as.numeric(as.character(blockJ$daJ))+as.numeric(as.character(blockJ$daJerr)))

