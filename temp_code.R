
# droprate preprose(serialization)
step = 4
map = c(0, c(1:step) * (1/step))
droprate = sapply(1:dim(droprate)[2], function(i){
  temp = droprate[, i]
  for(j in 2:step){
    temp[which(temp <= map[j]) && temp > map[j-1]] = map[j]
  }
  return(temp)
})


### droprate based on group
map = unique(gt)
group_id = lapply(1:length(map), function(i){
  res = which(gt == map[i])
  return(res)
})
group_droprate = sapply(group_id, function(i){
  clust_droprate = colSums(droprate[i, ] / length(i))
  return(clust_droprate)
})


# visualize
meandrop = colSums(droprate) / N
sort_geneid = order(colSums(droprate), decreasing=TRUE)
plot(meandrop[sort_geneid])
dev.new()
plot(mu_g[sort_geneid])
sort_geneid = order(group_droprate[, 3], decreasing=TRUE)
for(i in 1:length(map)){
  png(paste0("../MGGE/result/figure/cluster_meandrop/cluster_meandrop of ", map[i], ".png"), width=600, height=600)
  #dev.new()
  plot(group_droprate[sort_geneid, i], main=paste0("gene\'s mean droprate of cluster ", map[i]),
       xlab = paste0("sorted gene index of cluster ", map[3]), ylab = paste0("gene\'s mean droprate"))
  dev.off()
}



sort_geneid = order(colSums(droprate), decreasing=TRUE)
### cell's droprate based on group
Nsp = 1
group = 2
sample_id = sample(group_id[[group]], Nsp, replace=FALSE)
#sort_geneid = order(group_droprate[, 2], decreasing=TRUE)
for(i in 1:Nsp){
  if(i == 1){
    #png(paste0("result/figure/cell drop by cluster/mean droprate of cluster ", map[group], ".png"), width=600, height=600)
    dev.new()
    plot(group_droprate[sort_geneid, group], main=paste0("mean droprate of cluster ", map[group], " -- index based on cluster ", map[1]),
         xlab = paste0("sorted gene index of cluster ", map[1]), ylab = paste0("mean droprate of cluster ", map[group]))
    #dev.off()
  }
  
  #png(paste0("result/figure/cell drop by cluster/cell droprate of cell ", i, " from cluster ", map[group], ".png"), width=600, height=600)
  dev.new()
  plot(droprate[sample_id[i], sort_geneid], main=paste0("droprate of cell ", i, " from cluster ", map[group]),
       xlab = paste0("sorted gene index of cluster ", map[1]), ylab = paste0("droprate of cell ", i, " from cluster ", map[group]))
  #dev.off()
}

dev.new()
plot((colSums(droprate)/dim(droprate)[1])[sort_geneid])
var_droprate = apply(droprate, 2, function(i){
  res = sqrt(var(i))
  return(res)
})
dev.new()
plot(var_droprate[sort_geneid])



### reading weighted res
nmi = matrix(0, 3, 3)
rownames(nmi) = as.character(c(.5, 1, 3))
colnames(nmi) = as.character(c(.5, 1, 3))

for(sigmac in c(.5, 1, 3)){
  for(sigmag in c(.5, 1, 3)){
    try({
      res = readRDS(paste0("../PDGE/result/mouse1c_", sigmac, '_g_', sigmag, ".rds"))
      res = res$nmi[301, 1]
      nmi[as.character(sigmac), as.character(sigmag)] = res
    })
  }
}


### visualize
dev.new()
plot(res$lambda1)
dev.new()
plot(res$lambda2)
dev.new()
plot(res$J_HE)
dev.new()
plot(res$J_LE)
dev.new()
plot(res$nmi[,1])
title("nmi_H")
dev.new()
plot(res$nmi[,2])
title("nmi_Lw")
dev.new()
plot(res$nmi[,3])
title("nmi_spW")
dev.new()
plot(res$J_DR[,1])
title("cost_low drop")
dev.new()
plot(res$J_DR[,2])
title("cost_mid drop")
dev.new()
plot(res$J_DR[,3])
title("cost_high drop")
dev.new()
plot(res$weight[,1])
title("weight_low drop")
dev.new()
plot(res$weight[,2])
title("weight_mid drop")
dev.new()
plot(res$weight[,3])
title("weight_high drop")


# weigth of DR
data = "mouse1"
png(paste0("result/figure/weight and lost of PDGE/NMI of W dist ", data, ".png"), width=600, height=600)
plot(res$nmi[,2], main="NMI of W dist", ylab = "NMI", xlab = "iteration")
dev.off()
png(paste0("result/figure/weight and lost of PDGE/Low drop lost of ", data, ".png"), width=600, height=600)
plot(res$J_DR[, 1], main = "Lost of Low drop", ylab = "lost", xlab = "iteration")
dev.off()
png(paste0("result/figure/weight and lost of PDGE/Mid drop lost of ", data, ".png"), width=600, height=600)
plot(res$J_DR[, 2], main = "Lost of Mid drop", ylab = "lost", xlab = "iteration")
dev.off()
png(paste0("result/figure/weight and lost of PDGE/High drop lost of ", data, ".png"), width=600, height=600)
plot(res$J_DR[, 3], main = "Lost of High drop", ylab = "lost", xlab = "iteration")
dev.off()
png(paste0("result/figure/weight and lost of PDGE/Low drop weight of ", data, ".png"), width=600, height=600)
plot(res$weight[, 1], main = "Weight of Low drop", ylab = "weight", xlab = "iteration")
dev.off()
png(paste0("result/figure/weight and lost of PDGE/Mid drop weight of ", data, ".png"), width=600, height=600)
plot(res$weight[, 2], main = "Weight of Mid drop", ylab = "weight", xlab = "iteration")
dev.off()
png(paste0("result/figure/weight and lost of PDGE/High drop weight of ", data, ".png"), width=600, height=600)
plot(res$weight[, 3], main = "Weight of High drop", ylab = "weight", xlab = "iteration")
dev.off()