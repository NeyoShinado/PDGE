clustering_update <- function(X, K, npc, lambda1, lambda2, W=NULL, V=NULL, drop_id=NULL,
                              H=NULL, Dc=NULL, Dg=NULL, S=NULL, iteration=50, 
                              thre=1e-4, out_file="result/clust_res/res.rds", 
                              res_save=TRUE){
  cat("## clustering vars updating...\n")
  
  # initialization
  M = length(X)
  N = dim(W)[1]
  J_set = NULL
  J_DR = matrix(0, iteration, M)
  J_HE = NULL
  J_LE = NULL
  par1 = NULL
  par2 = NULL
  nmi = matrix(0, 1, 3)
  colnames(nmi) = c("nmiH", "nmiW_L", "nmiW_sp")
  Xe = c(sqrt(eigen(X[[1]] %*% t(X[[1]]))$values)[1], sqrt(eigen(X[[2]] %*% t(X[[2]]))$values)[1], 
         sqrt(eigen(X[[3]] %*% t(X[[3]]))$values)[1])
  Se = eigen(S)$values[1]
  Ds = diag(colSums(S))
  L = Ds - S
  alpha = matrix(1/3, iteration, 3)
  Lid = drop_id[1]
  Mid = drop_id[2]
  Hid = drop_id[3]
  #drop_H = droprate[, Hid]
  rm(droprate)
  
  
  for(iter in 1:iteration){
    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))
    cat("### updating W...\n")
    #
    LH = H %*% t(H)
    pr = alpha[iter, 1] * V[[1]] %*% t(V[[1]]) + alpha[iter, 3] * V[[3]] %*% t(V[[3]])
    sec_ord_g = lambda1 * LH %*% W + W %*% pr + alpha[iter, 2] * Dc^2 %*% W %*% V[[2]] %*% Dg^2 %*% t(V[[2]])
    nW = alpha[iter, 1] * X[[1]] %*% t(V[[1]]) + alpha[iter, 2] * Dc^2 %*% X[[2]] %*% Dg^2 %*% t(V[[2]]) + 
      alpha[iter, 3] * X[[3]] %*% t(V[[3]])
    W = W * ( nW/ sec_ord_g)
    rm(pr, nW)
    #g_W = -Dc^2 %*% X %*% t(V) + sec_ord_g#- 0.36s
    ## the N*k * N*k Hessian matrix is corresponse to Wij which expand by col
    # Hessian = kronecker(Matrix(diag(1, npc), sparse=TRUE), Matrix(lambda2 * LH, sparse=TRUE)) + 
    #   kronecker(V %*% t(V), Matrix(Dc^2, sparse=TRUE))  #- sparse 0.66s
    #      Hessian = Matrix(Hessian, sparse=TRUE)   #- 10s
    # sec_ord_g = Hessian %*% matrix(W, ncol=1)   
    
    
    cat("### updating V...\n")
    V[[1]] = V[[1]] * ((t(W) %*% X[[1]]) / (t(W) %*% W %*% V[[1]]))
    V[[2]] = V[[2]] * ((t(W) %*% Dc^2 %*% X[[2]] %*% Dg^2) / (t(W) %*% Dc^2 %*% W %*% V[[2]] %*% Dg^2))
    V[[3]] = V[[3]] * ((t(W) %*% X[[3]]) / (t(W) %*% W %*% V[[3]]))
    
    #* 5\ can not normalize W & V
    #norms = NormalizeUV(W, V)
    #W = norms$U
    #V = norms$V
    #rm(norms)
    
    
    cat("### updating H...\n")
    #d = -2 * W %*% t(W)
    d = rowSums(W^2)
    d = matrix(d, nrow=N, ncol=N)
    d = d + t(d) - 2*W %*% t(W)
    
    #- H = eigen(0.5*lambda1 * d + lambda2 * S)$vectors[, (N-K+1):N]
    H = H * ((lambda2 * S %*% H) / ((lambda2 * Ds + lambda1 / 2 * d) %*% H))
    
    norms = rowSums(H)    #* 5\row based H normalize
    norms[norms == 0] = 1e-10
    H = H / matrix(norms, N, K)
    
    
    cat('### updating alpha...\n')
    J_DR[iter, ] = t(sapply(1:M, function(i){
      if(i == 2){
        cost = norm(Dc^2 %*% (X[[i]] - W %*% V[[i]]) %*% Dg^2, 'F')^2
      }else{
        cost = norm(X[[i]] - W %*% V[[i]], 'F')^2
      }
      
      return(cost)
    }))
    
    alpha[iter, ] = t(sapply(1:M, function(i){
      res = 0.5 / sqrt(J_DR[iter, i])
      return(res)
    }))
    alpha[iter, ] = alpha[iter, ] / sum(alpha[iter, ])
    
    # integration update
    if(FALSE){
      cat('### updating S...\n')
      d_H = matrix(rowSums(H^2), N, N) + t(matrix(rowSums(H^2), N, N)) + 2* H %*% t(H)
      A = lapply(1:length(local_sim), function(i){
        local_sim[[i]] * alpha[i]
      })
      A = Reduce('+', A)
      S = t(sapply(1:N, function(i){
        res = rep(0, N)
        id = which(A[i, ] > 0)
        #* 10\update L with H
        ad = (A[i, id]) / sum(alpha)   #* not H constraint
        #ad = (A[i, id] - 0.25*lambda2*d_H[i, id]) / sum(alpha)
        res[id] = proj_sim(ad)
        return(res)
      }))
      S = (S + t(S)) / 2
      
    }
    
    
    #** 6\update parameters    option one
    if(iter %% 10 == 0){
      He = sqrt(eigen(H %*% t(H))$values)[1]
      Ve = c(sqrt((eigen(V[[1]] %*% t(V[[1]]))$values)[1]), sqrt((eigen(V[[2]] %*% t(V[[2]]))$values)[1]),
             sqrt((eigen(V[[3]] %*% t(V[[3]]))$values)[1]))
      lambda1 = sum(alpha[iter, ] * Xe) / sum(He + alpha[iter, ] * Ve)
      
      de = eigen(d)$values[1]
      lambda2 = (lambda1 * de) / (2 * Se)
      #par1[iter] = lambda1
      #par2[iter] = lambda2
      par1[(iter-9):iter] = lambda1
      par2[(iter-9):iter] = lambda2
    }
    
    
    # cost calculation
    # if verbose
    J_LE[iter] = sum(diag(t(H) %*% L %*% H))
    J_HE[iter] = sum(diag(t(W) %*% LH %*% W))
    J_set[iter] = alpha[iter, ] * J_DR[iter, ] + lambda1 * J_HE[iter] + lambda2 * J_LE[iter]
    cat("### Current cost:", J_set[iter], "\n")
    
    
    if(iter > 1){
      var_J = abs(J_set[iter] - J_set[iter-1])
      
      # convergence check
      if(var_J <= thre){
        cluster = sapply(1:N, function(i){
          which.max(H[i, ])
        })
        
        res = list(W = W, V = V, H = H, cluster = cluster,
                   J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, 
                   lambda1 = par1, lambda2 = par2, alpha = alpha, S = S)
        if(res_save){
          saveRDS(res, out_file)
        }
        return(res)
      }
    }
    
    
    # recording nmi
    #cluster = specc(as.matrix(d), K)@.Data
    #nmiW_sp = NMI(cluster, gt)
    nmiW_sp = 0
    
    F = eigen(d)$vectors
    F = F[, (N-K+1):(N)]
    clust = kmeans(F, K)$cluster
    nmiW_L = NMI(clust, gt)
    
    cluster = sapply(1:N, function(i){
      which.max(H[i, ])
    })
    nmiH = NMI(cluster, gt)
    
    nmi = rbind(nmi, c(nmiH, nmiW_L, nmiW_sp))
    
  }
  
  
  res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc), Dg = diag(Dg),
             J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, nmi = nmi,
             lambda1 = par1, lambda2 = par2, alpha = alpha, dw = d)   #* lambda save
  if(res_save){
    saveRDS(res, out_file)
  }
  return(res)
  
}