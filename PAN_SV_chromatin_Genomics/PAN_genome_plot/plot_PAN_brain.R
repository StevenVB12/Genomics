library(rtracklayer)


# interval
start = 568111292
end = 460265248

start =  568111292
end-start

startP = 1000
endP = 201000


minLemgth = 100
minLemgthTE = 100
maxLemgthTE = 10000

genomes <- c(6,5,3,1,4,2)
cols <- c('pink','pink','cornflowerblue','orange')
names <- c('H. e. cha1', 'H. e. cha2', 'H. e. dem', 'H. m. ros')

###
# ATAC data
###

range <- paste('pan:',568111292+startP,'-',568111292+endP,sep='')

BRAIN_cha1 <- import.bw("bw/brain_H_cha1_normalized_mean_start_end_pan.rowsort.clean.bw", which = GRanges(range))
BRAIN_cha2 <- import.bw("bw/brain_H_cha2_normalized_mean_start_end_pan.rowsort.clean.bw", which = GRanges(range))
BRAIN_dem <- import.bw("bw/brain_H_e_dem_normalized_mean_start_end_pan.rowsort.clean.bw", which = GRanges(range))
BRAIN_ros <- import.bw("bw/brain_H_m_ros_normalized_mean_start_end_pan.rowsort.clean.bw", which = GRanges(range))

bw_list <- list(BRAIN_cha1, BRAIN_cha2, BRAIN_dem, BRAIN_ros)

genomes <- c(2,3,1,4)

unique_dem <- read.table('blocks_unique_1.txt')
unique_cha1_poly <- read.table('blocks_unique_2and3.txt')
unique_cha2_poly <- read.table('blocks_unique_2and3.txt')
unique_ros <- read.table('blocks_unique_4.txt')
uniqueList <- list(unique_cha1_poly, unique_cha2_poly, unique_dem, unique_ros)

TEdem <- read.table('eratovsmelp.TE.unique.bed')
TEmelp <- read.table('melp.TE.unique.bed')
TEcha1 <- read.table('cha1.TE.unique.bed')
TEcha2 <- read.table('cha2.TE.unique.bed')

colnames(TEdem) <- c('pan','start','end')
colnames(TEmelp) <- c('pan','start','end')
colnames(TEcha1) <- c('pan','start','end')
colnames(TEcha2) <- c('pan','start','end')

TE_list <- list(TEcha1, TEcha2, TEdem, TEmelp)

# png('Hmel218003_PAN_brain.png', width=5000, height=5000)
pdf('Hmel218003_PAN_brain.pdf', width=50, height=50)

layout(matrix(c(1:16), nrow=16, byrow=TRUE), height = c(0.2,0.1,0.2,1,0.1,0.2,1,0.1,0.2,1,0.1,0.2,1,0.1,0.2,1))
# layout.show(n=25)

par(mar = c(0,4,0,4), oma=c(2,0,2,0))

###
# Annotations
###
# 
# # grik2
# 
# ANNOT_start <- read.table("grik2_start_pan.txt", sep='\t', h=T)
# ANNOT_end <- read.table("grik2_end_pan.txt", sep='\t', h=T)
# 
# ANNOT <- as.data.frame(cbind(ANNOT_start[,2], ANNOT_end[,2]))
# names(ANNOT) <- c("con_start", "con_end")
# 
# plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)
# 
# for (g in 1:nrow(ANNOT)){
#   rect(ANNOT$con_start[g], 0.2, ANNOT$con_end[g], 0.8, col = "black", border = "black", lwd=5)
# }
# rect(min(ANNOT$con_start), 0.5, max(ANNOT$con_end), 0.5, col = "black", border = "black", lwd=5)
# 
# # regulcalcin
# 
# ANNOT_start <- read.table("regulcalcin_start_pan.txt", sep='\t', h=T)
# ANNOT_end <- read.table("regulcalcin_end_pan.txt", sep='\t', h=T)
# 
# ANNOT <- as.data.frame(cbind(ANNOT_start[,2], ANNOT_end[,2]))
# names(ANNOT) <- c("con_start", "con_end")
# 
# # plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)
# 
# for (g in 1:nrow(ANNOT)){
#   rect(ANNOT$con_start[g], 0.2, ANNOT$con_end[g], 0.8, col = "black", border = "black", lwd=5)
# }
# rect(min(ANNOT$con_start), 0.5, max(ANNOT$con_end), 0.5, col = "black", border = "black", lwd=5)
# 

plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)

genesPAN <- read.table('PAN_H_charithonia_genes_eratoTransfer.gff', h=F, sep= '\t')

genesPAN <- subset(genesPAN, (genesPAN$V5-genesPAN$V4 >= minLemgth) & (genesPAN$V4 > 568111292+startP) & (genesPAN$V5 < 568111292+endP))

colnames(genesPAN) <- c("contig", "HGC", "type", "con_start", "con_end", "dot", "str", "unk", "descr")

# for (g in 1:nrow(genesPAN)){
#   rect(genesPAN$V4[g], 0.4, genesPAN$V5[g], 0.6, col = "black", border = "black", lwd=5)
# }




for (g in 1:nrow(genesPAN)){
  if (genesPAN$type[g] == "gene" && genesPAN$str[g] == "-") 
    rect(genesPAN$con_start[g], 0.75, genesPAN$con_end[g], 0.75, col = NULL, border = "black")
  if (genesPAN$type[g] == "exon" && genesPAN$str[g] == "-") 
    rect(genesPAN$con_start[g], 0.5, genesPAN$con_end[g], 1, col = "black", border = "black", lwd=0.3)
  if (genesPAN$type[g] == "gene" && genesPAN$str[g] == "+") 
    rect(genesPAN$con_start[g], 0.25, genesPAN$con_end[g], 0.25, col = NULL, border = "black")
  if (genesPAN$type[g] == "exon" && genesPAN$str[g] == "+") 
    rect(genesPAN$con_start[g], 0, genesPAN$con_end[g], 0.5, col = "black", border = "black", lwd=0.3)
}
###
# TEs
###

plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)

TEs <- read.table("total.TE.bed")

colnames(TEs) <- c('scaf', 'start', 'end')
TEs <- subset(TEs, (TEs$end-TEs$start >= minLemgthTE) & (TEs$end-TEs$start <= maxLemgthTE) & (TEs$start > 568111292+startP) & (TEs$end < 568111292+endP))

for(i in 1:nrow(TEs)){
  rect(TEs$start[i], 0, TEs$end[i], 1, col = 'darkred', border = NA)
}

###
# PAN genome
###
plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)

top = 1
bot = 0

rect(568111292+startP, 0, 568111292+endP, 1, col = "black", border = NA)

for(e in 1:length(genomes)){
  
  unique <- uniqueList[[e]]
  colnames(unique) <- c('scaf', 'start', 'end')
  unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > 568111292+startP) & (unique$end < 568111292+endP))
  
  if(nrow(unique) > 0){ 
    for(i in 1:nrow(unique)){
      rect(unique$start[i], bot, unique$end[i], top, col = cols[e], border = NA)
    }
  }
  
  if(e == 1){
    unique <- read.table('blocks_unique_2_vs_3.txt')
    colnames(unique) <- c('scaf', 'start', 'end')
    unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > 568111292+startP) & (unique$end < 568111292+endP))
    
    if(nrow(unique) > 0){ 
      for(i in 1:nrow(unique)){
        rect(unique$start[i], bot2, unique$end[i], top2[e], col = adjustcolor('pink', alpha = 1), border = NA)
      }
    }
  }
  
  if(e == 2){
    unique <- read.table('blocks_unique_3_vs_2.txt')
    colnames(unique) <- c('scaf', 'start', 'end')
    unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > 568111292+startP) & (unique$end < 568111292+endP))
    
    if(nrow(unique) > 0){ 
      for(i in 1:nrow(unique)){
        rect(unique$start[i], bot2, unique$end[i], top2[e], col = adjustcolor('pink', alpha = 1), border = NA)
      }
    }
  }
}

###
# ATAC
###

bot2 = 0
top2 = c(150,150,150,150)

for(e in 1:length(bw_list)){

  unique <- uniqueList[[e]]
  colnames(unique) <- c('scaf', 'start', 'end')
  unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > 568111292+startP) & (unique$end < 568111292+endP))

  plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(bot2,top2[e]), axes=FALSE, ann=FALSE)
  
  if(nrow(unique) > 0){ 
    for(i in 1:nrow(unique)){
      rect(unique$start[i], bot2, unique$end[i], top2[e], col = adjustcolor(cols[e], alpha = 0.5), border = NA)
    }
  }
  
  if(e == 1){
    unique <- read.table('blocks_unique_2_vs_3.txt')
    colnames(unique) <- c('scaf', 'start', 'end')
    unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > 568111292+startP) & (unique$end < 568111292+endP))
    
    if(nrow(unique) > 0){ 
      for(i in 1:nrow(unique)){
        rect(unique$start[i], bot2, unique$end[i], top2[e], col = adjustcolor('gray', alpha = 0.5), border = NA)
      }
    }
  }
  
  if(e == 2){
    unique <- read.table('blocks_unique_3_vs_2.txt')
    colnames(unique) <- c('scaf', 'start', 'end')
    unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > 568111292+startP) & (unique$end < 568111292+endP))
    
    if(nrow(unique) > 0){ 
      for(i in 1:nrow(unique)){
        rect(unique$start[i], bot2, unique$end[i], top2[e], col = adjustcolor('gray', alpha = 0.5), border = NA)
      }
    }
  }
  
  par(new=T)
  plot(0.5*(start(bw_list[[e]]) + end(bw_list[[e]])), bw_list[[e]]$score, type='l', xlim = c(568111292+startP,568111292+endP), ylim = c(bot2,top2[e]), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none')
  
  # mtext(names[e], side = 1, cex=5, padj = -2, las = 1, adj=0)
  
  ##
  plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)
  
  TE_spec <- TE_list[[e]]
  
  TE_spec <- subset(TE_spec, (TE_spec$end-TE_spec$start >= minLemgthTE) & (TE_spec$end-TE_spec$start <= maxLemgthTE) & (TE_spec$start > 568111292+startP) & (TE_spec$end < 568111292+endP))
  
  if(nrow(TE_spec) > 0){ 
    for(i in 1:nrow(TE_spec)){
      rect(TE_spec$start[i], 0, TE_spec$end[i], 1, col = 'darkred', border = NA)
    }
  }
  ##
  
  blocks <- read.table(paste(genomes[e],"_blocks_intervals.corr.bed",sep=""), h = F)
  colnames(blocks) <- c('scaf', 'start', 'end')
  blocks <- subset(blocks, (blocks$end-blocks$start >= 20) & (blocks$start > 568111292+startP) & (blocks$end < 568111292+endP))
  
  plot(NULL, xlim=c(568111292+startP,568111292+endP), ylim = c(bot,top), axes=FALSE, ann=FALSE)
  rect(startP+568111292, 0, endP+568111292, 1, col = "grey90", border = NA)
  
  if(nrow(blocks) > 0){ 
    for(i in 1:nrow(blocks)){
      rect(blocks$start[i], bot, blocks$end[i], top, col = adjustcolor(cols[e], alpha = 1), border = NA)
    }
  }
}

plot(NULL, xlim=c(startP,endP), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1, at = seq(startP,endP, by=10000), labels = NA, line =-70, lwd = 3, lwd.ticks = 3, tck = -0.05)
# axis(1, at = seq(startP,endP, by=100000), labels = seq(startP/1000000,endP/1000000, by=0.1), line =-3, cex.axis = 6, lwd=0)

dev.off()

