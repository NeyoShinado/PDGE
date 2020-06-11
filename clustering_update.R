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
  #par_sigma = NULL
  #par_beta = NULL
  nmi = c()
  #nmi = matrix(0, iteration, 3)
  # eigen val
  Xe = c(sqrt(eigen(X[[1]] %*% t(X[[1]]))$values)[1], sqrt(eigen(X[[2]] %*% t(X[[2]]))$values)[1], 
         sqrt(eigen(X[[3]] %*% t(X[[3]]))$values)[1])
  #Xn1 = sapply(1:3, function(i){
  #  res = norm(X[[i]], '1')
  #  return(res)
  #})
  
  Se = eigen(S)$values[1]
  # var init
  #E = W
  #theta = matrix(0, N, npc)
  #sigma = 5e+1
  #beta = 1    #* par of sparse norm
  Ds = diag(colSums(S))
  L = Ds - S
  alpha = matrix(1/3, iteration, 3)
  Lid = drop_id[1]
  Mid = drop_id[2]
  Hid = drop_id[3]
  ignore = sapply(1:M, function(i){
    if(length(drop_id[[i]]) < 50){
      return(i)
    }else{
      return(NA)
    }
  })
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
    W = W * (nW / sec_ord_g)
    #W = W * ((nW + sigma*E + theta) / (sec_ord_g + sigma*W))
    rm(pr, nW)
    
    # 5\ unitize
    #norms = rowSums(W)
    #norms[norms == 0] = 1e-10
    #W = W / matrix(norms, N, dim(W)[2])
    
    #g_W = -Dc^2 %*% X %*% t(V) + sec_ord_g#- 0.36s
    ## the N*k * N*k Hessian matrix is corresponse to Wij which expand by col
    # Hessian = kronecker(Matrix(diag(1, npc), sparse=TRUE), Matrix(lambda2 * LH, sparse=TRUE)) + 
    #   kronecker(V %*% t(V), Matrix(Dc^2, sparse=TRUE))  #- sparse 0.66s
    #      Hessian = Matrix(Hessian, sparse=TRUE)   #- 10s
    # sec_ord_g = Hessian %*% matrix(W, ncol=1)
    
    ##* sparse regularize
    #cat("### updating E and theta...\n")
    #E = soft(W - theta/sigma, beta/sigma)
    #sigma = 1.618*sigma
    #theta = theta + sigma * (E - W)
    
    
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
    
    H = eigen(0.5*lambda1 * d + lambda2 * L)$vectors[, (N-K+1):N]
    #H = H * ((lambda2 * S %*% H) / ((lambda2 * Ds + lambda1 / 2 * d) %*% H))
    
    #norms = rowSums(H)    #* 5\row based H normalize
    #norms[norms == 0] = 1e-10
    #H = H / matrix(norms, N, K)
    
    
    cat('### updating alpha...\n')
    J_DR[iter, ] = t(sapply(1:M, function(i){
      if(i == 2){
        cost = norm(Dc^2 %*% (X[[i]] - W %*% V[[i]]) %*% Dg^2, 'F')^2
      }else{
        cost = norm(X[[i]] - W %*% V[[i]], 'F')^2
      }
      
      return(cost)
    }))
    
    # veiw check
    alpha[iter, ] = t(sapply(1:M, function(i){
      if(i %in% ignore){
        return(0)
      }else{
        res = 0.5 / sqrt(J_DR[iter, i])
        return(res)
      }
    }))
    # 4\unitize weight
    #alpha[iter, ] = alpha[iter, ] / sum(alpha[iter, ])
    
    
    #** 6\update parameters    option one
    if(iter %% 10 == 0){
      He = sqrt(eigen(H %*% t(H))$values)[1]
      Ve = c(sqrt((eigen(V[[1]] %*% t(V[[1]]))$values)[1]), sqrt((eigen(V[[2]] %*% t(V[[2]]))$values)[1]),
             sqrt((eigen(V[[3]] %*% t(V[[3]]))$values)[1]))
      lambda1 = sum(alpha[iter, ] * Xe) / sum(He + alpha[iter, ] * Ve)
      
      de = eigen(d)$values[1]
      lambda2 = (2 * Se) / (lambda1 * de)
      par1[iter] = lambda1
      par2[iter] = lambda2
      
      #* sparse norm
      #beta = sum(alpha[iter, ] * Xn1) / sum(alpha[iter, ] * Ve)
      #par_beta[iter] = beta
      #par_sigma[iter] = sigma
      
      #par1[(iter-9):iter] = lambda1
      #par2[(iter-9):iter] = lambda2
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
        
        res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc), Dg = diag(Dg),
                   J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, nmi = nmi,
                   lambda1 = par1, lambda2 = par2, alpha = alpha, dw = d)   #* lambda save
        if(res_save){
          saveRDS(res, out_file)
        }
        return(res)
      }
    }
    
    
    # recording nmi
    if(FALSE){
      cluster = specc(as.matrix(d), K)@.Data
      nmiW_sp = NMI(cluster, gt)
      
      F = eigen(d)$vectors
      F = F[, (N-K):(N-1)]
      clust = kmeans(F, K)$cluster
      nmiW_L = NMI(clust, gt)
      
      cluster = sapply(1:N, function(i){
        which.max(H[i, ])
      })
      nmiH = NMI(cluster, gt)
    }
    
    cluster = kmeans(H, K, nstart = 25, iter.max = 100)$cluster
    nmiH = NMI(cluster, gt)
    nmi = c(nmi, nmiH)
    #nmi[iter, ] = c(nmiH, nmiW_sp, nmiW_L)
    
  }
  
  
  res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc), Dg = diag(Dg),
             J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, nmi = nmi,
             lambda1 = par1, lambda2 = par2, alpha = alpha, dw = d)   #* lambda save
  if(res_save){
    saveRDS(res, out_file)
  }
  return(res)
}


soft <- function(x, theta){
  if(max(theta == 0)){
    res = x
    return(x)
  }else{
    res = abs(x) - theta
    res[which(res<0)] = 0
    res = sign(x) * res
    return(res)
  }
}