setwd("E:/Project/Paper_debug/Clustering algorithm/PDGE/")


source("npc_cal.R")
#source("proj_sim.R")
source("constructW.R")
source("PDGE_update.R")
#source("NormalizeUV.R")
source("data_normalize.R")
source("imputing_update.R")
source("clustering_update.R")
#source("construct_locsim.R")
#source("localsim_integrated_update.R")
#source("informative_gene_selection.R")


library(Matrix)
library(parallel)
library(aricode)
library(Seurat)
library(R.matlab)
library(ggplot2)


# dataset info
datas = c("Biase/Biase.rds", "Deng_GSE45719/mouse_embryo.rds", "Zeisel/Zeisel.rds", 
          "mouse1/mouse1.rds", "mouse2/mouse2.rds", "human3/human3.rds", "Hrvatin/Hrvatin.rds",
          "human1/human1.rds", "human2/human2.rds", "human4/human4.rds")
data_path = "E:/Project/dataset/Bioinformatics/scRNA/Selected data/"
for(i in c(4)){
  for(Lmud in c(0.4)){
    for(Mmud in c(0.6)){
      for(Hmud in c(0.85)){
      sigmac = 0.5
      sigmag = 0.5

    try({
    X = readRDS(paste0(data_path, datas[i]))
    gt = X$gt
    X = X$expr     # for Baron dataset
    #gt = X$cluster
    #X = t(X@assays@data$rawcounts)
    N = nrow(X)
    P = ncol(X)
    message(paste0("## Loading raw data of ",N , "*", P, " and labels...\n"))
    
    
    ## gene filter & normalization
    message("## Data nomalization and log-transformation...\n")
    lg_X = data_normalize(X, N, P, gt, mu_probs=0.5, cv_probs=0.2)   #* 0.2 by default with gene select  0.1*round(log10(N))
    gt = lg_X$gt

    
    mu_g = lg_X$mu_g
    sd_g = lg_X$sd_g
    lg_X = as.matrix(lg_X$count_hv)
    
    
    ##* cluster number estimation
    #K = cluster_number_estimation(lg_X)
    K = length(unique(gt))
    message("## Estimated cluster number:", K, "\n")
    
    
    #* construct W & neighbors
    # NN = 5 for 10^2, else 20 for 1o^3
    res = constructW(lg_X, 20, K)
    S = res$S
    neighbors = res$neighbors
    rm(res)
    L = diag(colSums(S)) - S
    H = eigen(L)$vectors
    H = H[, (N-K):(N-1)]
    cluster = kmeans(H, K, nstart=20)$cluster
    
    
    #* 14\ no gene selection
    if(FALSE){
      #lg_X = informative_gene_selection(lg_X, delta=0.7)
      data_object = CreateSeuratObject(counts = Matrix(t(lg_X), sparse=TRUE),
                                       project = "MGGE", min.cells = 3)
      data_object = NormalizeData(data_object, normalization.method = "LogNormalize", scale.factor = 10000)
      ## gene selection
      #** 2\select about thousand genes
      data_object = FindVariableFeatures(data_object, selection.method = "vst", nfeatures = 3000)
      gene_id = VariableFeatures(data_object)
      if(all(gene_id %in% colnames(lg_X))){
      }else{
        gene_id = gene_id[-which(!(gene_id %in% colnames(lg_X)))]
      }
      
      id = sapply(gene_id, function(i){
        res = which(colnames(lg_X) == i)
        return(res)
      })
      message("## Finally choose ", length(gene_id), " genes for MGGE...")
      
      lg_X = lg_X[, gene_id]
      mu_g = mu_g[id]
      sd_g = sd_g[id]
      rm(X, data_object)
    }
    
    
    
    ## dimension reduction(npc) calculation
    npc = npc_cal(as.matrix(lg_X), K, var_thre=0.8)
    message("## Caluculated dim reduction npc:", npc, "\n")
    
    
    ## cluster update
    # clustering update, 200 iteration by default
    #* stop update when Jbefore < Jafter
    files = list.files("./scImpute/", pattern="*.R")
    sapply(paste0("./scImpute/", files), source)
    

    # main function
    res <- var_update(lg_X, K, npc, S, neighbors, cluster, gt, mu_g, sd_g, sigmac=sigmac, sigmag=sigmag, lambda1=2, lambda2=2, 
                      iteration=1, clust_iteration=300, imp_iteration=3000, 
                      res_save=FALSE, Lmud=Lmud, Mmud=Mmud, Hmud=Hmud, data_name = strsplit(datas[i], split="/")[[1]][1])
    
    #clust = res$cluster
    #nmi = NMI(clust, gt)
    #ari = ARI(clust, gt)
    #res$nmi = nmi
    #res$ari = ari
    #message("## NMI: ", nmi, "   ARI: ", ari, "\n")
    
    
    # save res
    output = strsplit(datas[i], split="/")[[1]][1]
    
    saveRDS(res, paste0("result/test/", output, 'Lmud_', Lmud, '_Mmud_', Mmud, '_Hmud_', Hmud, '.rds'))
    })
  }
}
}
}


