# files <- list.files(path='msms_run3', pattern='msms', full.names = T)
# 
# table <- c()
# for(f in 1:length(files)){
#   if(file.info(files[f])$size !=0){
#   table1 <- read.table(files[f], h=T)
#   table <-rbind(table,table1)
#   }
# }

table <- read.table("msms/msms_run220319.txt.gz", h=T)
head(table)


migration <- unique(table$m23)
mchange <- unique(table$mchange)
mchange <- sort(mchange)
recs <- unique(table$rec)
split <- unique(table$split)
selStart <- unique(table$selStart)

tableS <- subset(table, table$selStart == 5e-01 & table$split == 0.25)
tableS$pos <- (tableS$pos-0.5)*10000
# plot(tableS$pos, tableS$Fst)
# 
library(plyr)
tableSM <- ddply(tableS, .(pos, rec, m23, mchange), summarize, Fst=mean(Fst), fd=mean(fd))
# tableSM$pos <- (tableSM$pos-0.5)*10000
# 
tableSM2 <- tableSM
tableSM2$pos <- tableSM2$pos*-1

tableSMM <- rbind(tableSM,tableSM2)
tableSMM <- subset(tableSMM, tableSMM$rec >0.000625)
# tableSMM$pos

tableSMM <- tableSMM[order(tableSMM$pos),]

recomb <- unique(tableSMM$rec)
colfunc <- colorRampPalette(c("green","purple"))
colList <- colfunc(length(recomb))


layout(matrix(c(1:6), nrow=3, byrow=TRUE))
# layout.show(n=2)

for(s in 1:length(mchange)){
  for(m in 1:length(migration)){
    plot(0, pch = "",xlim = c(-2000000,2000000), ylim = c(0,0.8), ylab = "Fst", lwd = 0.5, xlab = "Position", bty = "n", main = paste("m =", migration[m], "mch =", mchange[s], sep=' '))
    for(e in 1:length(recomb)){
      tb <- subset(tableSMM, tableSMM$rec == recomb[e] & tableSMM$m23 == migration[m] & tableSMM$mchange == mchange[s])
      # lo<-loess.smooth(tb$pos, tb$Fst, span = 1/5, degree = 1,
                       # family = c("symmetric", "gaussian"), evaluation = 50)
      # lines(lo, col = colList[e], lwd =2)
      par(new=T)
      plot(tb$pos, tb$Fst, col = colList[e], type='l', lwd =2, xlim = c(-2000000,2000000), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
    }
    plot(0, pch = "",xlim = c(-2000000,2000000), ylim = c(0,0.8), ylab = "fd", lwd = 0.5, xlab = "Position", bty = "n", main = paste("m =", migration[m], "mch =", mchange[s], sep=' '))
    for(e in 1:length(recomb)){
      tb <- subset(tableSMM, tableSMM$rec == recomb[e] & tableSMM$m23 == migration[m] & tableSMM$mchange == mchange[s])
      # lo<-loess.smooth(tb$pos, tb$fd, span = 1/5, degree = 1,
                       # family = c("symmetric", "gaussian"), evaluation = 50)
      # lines(lo, col = colList[e], lwd =2)
      par(new=T)
      plot(tb$pos, tb$fd, col = colList[e], type='l', lwd =2, xlim = c(-2000000,2000000), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
    }
  }
}

##################
library(plotrix)

table <- subset(table, table$rec > 0.000625 & table$mchange>=0.00025)

migration <- unique(table$m23)
mchange <- unique(table$mchange)
mchange <- sort(mchange)
recs <- unique(table$rec)
split <- unique(table$split)
selStart <- unique(table$selStart)

selStartF = 5e-1
m23F = 40
mchangeF = 0.25
splitsF = 0.25

tableF <- subset(table, table$selStart == selStartF & table$m23 == m23F & table$mchange == mchangeF & table$split == splitsF)
tableF$pos <- (tableF$pos-0.5)*10000

tableFM <- ddply(tableF, .(pos, rec), summarize, Fstm=mean(Fst,na.rm=TRUE), fdm=mean(fd,na.rm=TRUE), Fstsd=sd(Fst,na.rm=TRUE), fdsd=sd(fd,na.rm=TRUE))
head(tableFM)

tableFM2 <- tableFM
tableFM2$pos <- tableFM2$pos*-1

tableFMM <- rbind(tableFM,tableFM2)


tableFMM <- tableFMM[order(tableFMM$pos),]

recomb <- unique(tableFMM$rec)
colfunc <- colorRampPalette(c("green","purple"))
colList <- colfunc(length(recomb))

layout(matrix(c(1:15), nrow=5, byrow=TRUE))
par(mai=c(0.5,0.5,0.2,0.2), oma=c(2,2,0,0)+0)


plot(0, pch = "",xlim = c(-2000000,2000000), ylim = c(0,0.8), ylab = "Fst", lwd = 0.5, xlab = "Position (bp)", bty = "n", main = "")
segments(-500000,0,-500000,0.8, col = "red", lty = 2, lwd =2)
for(e in 1:length(recomb)){
  tb <- subset(tableFMM, tableFMM$rec == recomb[e])
  # lo<-loess.smooth(tb$pos, tb$Fst, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$pos, tb$Fstm, col = colList[e], type='l', lwd =2, xlim = c(-2000000,2000000), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
}
plot(0, pch = "",xlim = c(-2000000,2000000), ylim = c(0,0.5), ylab = "fd", lwd = 0.5, xlab = "Position (bp)", bty = "n", main = "")
segments(-500000,0,-500000,0.5, col = "red", lty = 2, lwd =2)
for(e in 1:length(recomb)){
  tb <- subset(tableFMM, tableFMM$rec == recomb[e] )
  # lo<-loess.smooth(tb$pos, tb$fd, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$pos, tb$fdm, col = colList[e], type='l', lwd =2, xlim = c(-2000000,2000000), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
}

plot(0, pch = "",xlim = c(0,10), ylim = c(0,10), ylab = "", lwd = 0.5, xlab = "", bty = "n", main = "", axes=F)
legend(1, 10, legend= recs, col=colList, lty=1, cex=0.8, box.lty=0)



tableF <- subset(table, table$selStart == selStartF & table$mchange == mchangeF & table$split == splitsF & table$pos == -49.5)
tableF$pos <- (tableF$pos-0.5)*10000

tableFM <- ddply(tableF, .(m23, rec), summarize, Fstm=mean(Fst,na.rm=TRUE), fdm=mean(fd,na.rm=TRUE), Fstsd=sd(Fst,na.rm=TRUE), fdsd=sd(fd,na.rm=TRUE))
head(tableFM)

migration <- sort(unique(tableF$m23), decreasing=T)
colfunc <- colorRampPalette(c("green","purple"))
colList <- colfunc(length(migration))

plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.8), ylab = "Fst", lwd = 0.5, xlab = "r", bty = "n", main = "migration")
gradient.rect(0,0,0.08,-.02,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(migration)){
  tb <- subset(tableFM, tableFM$m23 == migration[e])
  # lo<-loess.smooth(tb$pos, tb$Fst, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$Fstm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  if(tb$m23 == m23F){
    par(new=T)
    plot(tb$rec, tb$Fstm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.5), ylab = "fd", lwd = 0.5, xlab = "r", bty = "n", main = "migration")
gradient.rect(0,0,0.08,-.0125,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(migration)){
  tb <- subset(tableFM, tableFM$m23 == migration[e] )
  # lo<-loess.smooth(tb$pos, tb$fd, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$fdm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  if(tb$m23 == m23F){
    par(new=T)
    plot(tb$rec, tb$fdm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,10), ylim = c(0,10), ylab = "", lwd = 0.5, xlab = "", bty = "n", main = "", axes=F)
legend(1, 10, legend= migration, col="black", lty=1:3, cex=0.8, box.lty=0)


tableF <- subset(table, table$m23 == m23F & table$mchange == mchangeF & table$split == splitsF & table$pos == -49.5)
tableF$pos <- (tableF$pos-0.5)*10000

tableFM <- ddply(tableF, .(selStart, rec), summarize, Fstm=mean(Fst,na.rm=TRUE), fdm=mean(fd,na.rm=TRUE), Fstsd=sd(Fst,na.rm=TRUE), fdsd=sd(fd,na.rm=TRUE))
head(tableFM)

selStarts <- sort(unique(tableFM$selStart), decreasing=T)
colfunc <- colorRampPalette(c("green","purple"))
colList <- colfunc(length(selStart))

plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.8), ylab = "Fst", lwd = 0.5, xlab = "r", bty = "n", main = "selection T")
gradient.rect(0,0,0.08,-.02,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(selStarts)){
  tb <- subset(tableFM, tableFM$selStart == selStarts[e])
  # lo<-loess.smooth(tb$pos, tb$Fst, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$Fstm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  if(tb$selStart == selStartF){
    par(new=T)
    plot(tb$rec, tb$Fstm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.5), ylab = "fd", lwd = 0.5, xlab = "r", bty = "n", main = "selection T")
gradient.rect(0,0,0.08,-.0125,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(selStart)){
  tb <- subset(tableFM, tableFM$selStart == selStarts[e])
  print(tb)
  # lo<-loess.smooth(tb$pos, tb$fd, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$fdm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  if(tb$selStart == selStartF){
    par(new=T)
    plot(tb$rec, tb$fdm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,10), ylim = c(0,10), ylab = "", lwd = 0.5, xlab = "", bty = "n", main = "", axes=F)
legend(1, 10, legend= selStarts*4000000, col="black", lty=1:5, cex=0.8, box.lty=0)



tableF <- subset(table, table$m23 == m23F & table$selStart == selStartF & table$split == splitsF & table$pos == -49.5)
tableF$pos <- (tableF$pos-0.5)*10000

tableFM <- ddply(tableF, .(mchange, rec), summarize, Fstm=mean(Fst,na.rm=TRUE), fdm=mean(fd,na.rm=TRUE), Fstsd=sd(Fst,na.rm=TRUE), fdsd=sd(fd,na.rm=TRUE))
head(tableFM)

mchanges <- sort(unique(tableFM$mchange), decreasing=T)
colfunc <- colorRampPalette(c("green","purple"))
colList <- colfunc(length(mchanges))

plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.8), ylab = "Fst", lwd = 0.5, xlab = "r", bty = "n", main = "migration T")
gradient.rect(0,0,0.08,-.02,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(mchanges)){
  tb <- subset(tableFM, tableFM$mchange == mchanges[e])
  # lo<-loess.smooth(tb$pos, tb$Fst, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$Fstm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  if(tb$mchange == mchangeF){
    par(new=T)
    plot(tb$rec, tb$Fstm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.5), ylab = "fd", lwd = 0.5, xlab = "r", bty = "n", main = "migration T")
gradient.rect(0,0,0.08,-.0125,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(mchanges)){
  tb <- subset(tableFM, tableFM$mchange == mchanges[e])
  print(tb)
  # lo<-loess.smooth(tb$pos, tb$fd, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$fdm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  if(tb$mchange == mchangeF){
    par(new=T)
    plot(tb$rec, tb$fdm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,10), ylim = c(0,10), ylab = "", lwd = 0.5, xlab = "", bty = "n", main = "", axes=F)
legend(1, 10, legend= mchanges*4000000, col="black", lty=1:6, cex=0.8, box.lty=0)



tableF <- subset(table, table$m23 == m23F & table$selStart == selStartF & table$mchange == mchangeF & table$pos == -49.5)
tableF$pos <- (tableF$pos-0.5)*10000

tableFM <- ddply(tableF, .(split, rec), summarize, Fstm=mean(Fst,na.rm=TRUE), fdm=mean(fd,na.rm=TRUE), Fstsd=sd(Fst,na.rm=TRUE), fdsd=sd(fd,na.rm=TRUE))
head(tableFM)

splits <- sort(unique(tableFM$split), decreasing=T)
colfunc <- colorRampPalette(c("green","purple"))
colList <- colfunc(length(splits))

plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.8), ylab = "Fst", lwd = 0.5, xlab = "r", bty = "n", main = "split T")
gradient.rect(0,0,0.08,-.02,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(splits)){
  tb <- subset(tableFM, tableFM$split == splits[e])
  # lo<-loess.smooth(tb$pos, tb$Fst, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$Fstm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  if(tb$split == splitsF){
    par(new=T)
    plot(tb$rec, tb$Fstm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.8), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,0.08), ylim = c(0,0.5), ylab = "fd", lwd = 0.5, xlab = "r", bty = "n", main = "split T")
gradient.rect(0,0,0.08,-.0125,col=smoothColors("green",38,"purple"),border=NA)
for(e in 1:length(splits)){
  tb <- subset(tableFM, tableFM$split == splits[e])
  print(tb)
  # lo<-loess.smooth(tb$pos, tb$fd, span = 1/5, degree = 1,
  # family = c("symmetric", "gaussian"), evaluation = 50)
  # lines(lo, col = colList[e], lwd =2)
  par(new=T)
  plot(tb$rec, tb$fdm, col = "black", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  if(tb$split == splitsF){
    par(new=T)
    plot(tb$rec, tb$fdm, col = "red", type='l', lty = e, lwd =2, xlim = c(0,0.08), ylim = c(0,0.5), ylab = "", xlab = "", axes=F)
  }
}
plot(0, pch = "",xlim = c(0,10), ylim = c(0,10), ylab = "", lwd = 0.5, xlab = "", bty = "n", main = "", axes=F)
legend(1, 10, legend= splits*4000000, col="black", lty=1:3, cex=0.8, box.lty=0)
# 
# 
# library(lattice)
# library(gridExtra)
# library(viridis)
# 
# 
# 
# trellis.par.set("axis.line", list(col=NA,lty=1,lwd=1))
# wireframe(Fst ~ rec*pos, data=tableSMM, 
#                          xlab = "Recombination rate",
#                          ylab = "Position",
#                          main = "",
#                          drape = TRUE,
#                          colorkey = TRUE,
#                          scales = list(arrows=FALSE,cex=.5, tick.number = 10, z = list(arrows=F), distance =c(1.5, 1.5, 1.5)),
#                          light.source = c(10,0,10),
#                          col.regions = inferno(100, alpha= 0.8),
#                          screen = list(z = -60, x = -60)
# )
# 
# trellis.par.set("axis.line", list(col=NA,lty=1,lwd=1))
# wireframe(fd ~ -rec*pos, data=tableSMM, 
#           xlab = "Recombination rate",
#           ylab = "Position",
#           main = "",
#           drape = TRUE,
#           colorkey = TRUE,
#           scales = list(arrows=FALSE,cex=.5, tick.number = 10, z = list(arrows=F), distance =c(1.5, 1.5, 1.5)),
#           light.source = c(10,0,10),
#           col.regions = inferno(100, alpha = 0.8),
#           screen = list(z = -60, x = -60)
# )
