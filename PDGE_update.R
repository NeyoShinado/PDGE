# clustering vars updating
library(Matrix)
library(bignmf)


var_update <- function(lg_X, K, npc, S, neighbors, guide_cluster, gt, sigmac=3, sigmag=2, lambda1=0.1, lambda2=0.1, 
                       thre_J=1e-4, drop_thre=0.5, iteration=1, 
                       clust_iteration=100, imp_iteration=100, 
                       output="result/", res_save=TRUE){
  # lg_X without pseudo count
  cat("## Initializing for program...\n")
  if(!dir.exists(paste0(output, "clust_res"))){
    dir.create(paste0(output, "clust_res"))
  }
  if(!dir.exists(paste0(output, "imp_res"))){
    dir.create(paste0(output, "imp_res"))
  }
  
  N = dim(lg_X)[1]
  P = dim(lg_X)[2]

  
  imp_iter = 1
  droprate = 0
  imp_X = 0
  
  
### imputing_update
  cat(paste0("## carrying ", imp_iter, "th cluster imputation...\n\n"))
  #*clust = clust_set[[i]]
  imp_res = imputing_update(log10(10^(lg_X) + 0.01), guide_cluster, neighbors, 
                            imp_iteration=imp_iteration, 
                            out_file=paste0(output, "imp_res/", imp_iter, "th_impres.rds"), 
                            res_save=res_save)
  
  #imp_X = log10(10^(imp_X) - 0.01)   # log-tran 1.01 for imputing, 1.0 for clustering
  
  ## arrange imp_output
  #local_sim = imp_res$local_sim
  imp_X = imp_res$imp_X
  droprate = imp_res$droprate
  meandrop = apply(droprate, 2, mean)
  # L\M\H means high\middle\low droprate
  # 0.6\0.9\1 percent of drop
  Lid = which(meandrop < 0.4)
  Mid = which(meandrop < 0.65 & meandrop >= 0.4)
  Hid = which(meandrop >= 0.65)
  drop_id = c(Lid, Mid, Hid)
  
  
  # divide gene
  split_data = list()
  split_data[[1]] = imp_X[, Lid]
  split_data[[2]] = imp_X[, Mid]
  split_data[[3]] = droprate[, Hid]
  
  
  #* 12\weighted cell
  Dc = rowSums(droprate)
  Dg = colSums(droprate)[Mid]
  Dc = Matrix(diag(exp(-sigmac * Dc / P)), sparse=TRUE)
  Dg = Matrix(diag(exp(-sigmag * Dg / N)), sparse=TRUE)    # scale factor 2 for small dataset and 3 for the big one
  
  
  #** 3\test of best local_sim / imp_local_sim
  if(FALSE){
    local_sim[[i]] = imp_res$local_sim
    local_sim[[i]] = t(sapply(1:N, function(j){
      res = rep(0, N)
      id = order(local_sim[[i]][j, ], decreasing=TRUE)[1:NN]
      res[id] = local_sim[[i]][j, ][id]
      return(res)
    }))
    local_sim[[i]]  = (local_sim[[i]] + t(local_sim[[i]])) / 2
  }
  
  
  ### Clustering update
  cat(paste0("## running ", imp_iter, "th vars update via clustering and imputation...\n"))

  
  # init vars
  #H = matrix(runif(dim(lg_X)[1]*K), N, K)
  H = eigen(diag(colSums(S)) - S)$vector
  H = H[, (N-K):(N-1)]
  H[H<=0] = 1e-10
  W = matrix(runif(N * npc), N, npc) 
  V = lapply(1:3, function(i){
    pi = dim(split_data[[i]])[2]
    Vi = matrix(runif(npc * pi), npc, pi)
    return(Vi)
  })
  
  
  
  #* 7\ using split imp_X & dropout
  clust_res = clustering_update(split_data, K, npc, lambda1, lambda2, W=W, V=V, drop_id=drop_id,
                                H=H, Dc=Dc, Dg = Dg, S=S, iteration=clust_iteration,
                                out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), 
                                res_save=res_save)
  #clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, Dc=Dc, L=L, iteration=clust_iteration,
  #                             out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), res_save=res_save)
  
  H = clust_res$H
  V = clust_res$V
  W = clust_res$W
  dw = clust_res$dw
  alpha = clust_res$alpha
  cluster = clust_res$cluster
  nmi = clust_res$nmi
  res_J = clust_res$J
  lambda1 = clust_res$lambda1
  lambda2 = clust_res$lambda2
  J_DR = clust_res$J_DR
  J_HE = clust_res$J_HE
  J_LE = clust_res$J_LE
  #neighbors = clust_res$neighbors
  
  
  res = list(cluster=cluster, imp_X = imp_X, guide_cluster = guide_cluster, lambda1 = lambda1, lambda2 = lambda2,
             droprate=droprate, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, J = res_J, nmi = nmi,
             W = W, V = V, H = H, Dc = Dc, Dg = Dg, S = S, dw = dw, weight = alpha)  
  
  #  res = list(cluster=cluster, neighbors=neighbors, 
  #             imp_X=imp_X, droprate=droprate, local_sim = local_sim, 
  #             J = res_J, W = W, V = V, H = H)
  cat("# MGGE iteration complete!\n")
  #saveRDS(res, paste0(output, "MGGE_res.rds"))
  return(res)
}


imputing_update <- function(lg_X, cluster, neighbors, drop_thre=0.5, out_file="result/imp_res/res.rds", ncores=1, imp_iteration=100, res_save=TRUE){
  # return log-tran X
  # lg 1.01
  # not parallel version on imputation
  imp_res = imputation_wlabel_model(lg_X, cluster, neighbors=neighbors, point=log10(1.01),
                                    drop_thre=drop_thre, ncores=ncores, imp_iteration=imp_iteration)
  imp_X = imp_res$count_imp
  droprate = imp_res$droprate
  local_sim = imp_res$local_sim
  if(res_save){
    saveRDS(imp_res, out_file)
  }
  rm(imp_res)
  
  return(list(imp_X = imp_X, local_sim = local_sim, droprate = droprate))
}

