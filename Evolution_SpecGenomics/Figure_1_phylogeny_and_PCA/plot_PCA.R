pca <- read.table("Herato_040518.HeratoALL.thindist1000.PCA.comp.out", h=F)

pca

pca$col <- "gray"

for(e in 1:nrow(pca)){
  if(pca$V12[e] == "cyrS" | pca$V12[e] == "cyrN"){
    pca$col [e] <- "orange"
  }
  if(pca$V12[e] == "himS"| pca$V12[e] == "himN"){
    pca$col [e] <- "red"
  }
  if(pca$V12[e] == "himS"){
    pca$col [e] <- "brown"
  }
  if(pca$V12[e] == "emmW"| pca$V12[e] == "emmE"){
    pca$col [e] <- "green"
  }
  if(pca$V12[e] == "favW"| pca$V12[e] == "favE"){
    pca$col [e] <- "blue"
  }
  if(pca$V12[e] == "pet"){
    pca$col [e] <- "black"
  }  
  if(pca$V12[e] == "che"){
    pca$col [e] <- "pink"
  }
}

plot(pca$V2, pca$V3, col = pca$col, cex = 3, pch=19)
