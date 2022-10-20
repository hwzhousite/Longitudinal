############## packages ###########
require(mixtools)
require(Matrix)
library(splines)
#load MDSP_gamma function

###################################
#x: matrix of beta
#m.iter: maximum number of iterations
#eps: tolerance level for stopping rule
#sp : spline type
#center : number of center (default 2)
###############################

MDSP_temp_gamma <- function (x, m, m.iter, eps, id = NULL, sp = "ns", center = 2){
  
  
  #################################
  # p: p_x
  # n: number of identity
  # g: matrix for gamma
  # sp_mat: matrix stores regression and classification information
  ###################################
  
  p <- nrow(x)
  n <- ncol(x)/m
  g <- matrix(0, nrow = m, ncol = center * p) 
  sp_mat <- data.frame(t = rep(seq(1,m),each =n), cl = rep(0, n*m))
  model <- list()
  
  
  if(is.null(id)){
    
    fix_id <- sample(1:n,1)
    
  }else{
    
    fix_id <- as.numeric(id)
    
  }
  
  fix_tag <- 1
  
  # Initialization of cluster labels 
  for (i in 1:m) {
    
    ind <- which(sp_mat[,"t"] == i)
    
    temp <- as.data.frame(t(x[,ind, drop = FALSE]))
    
    cl_temp <- kmeans(temp, center= center)
    
    cl <- cl_temp$cluster
    cent <- cl_temp$centers
    
    cl_group <- as.numeric(names(sort(table(cl[fix_id]), decreasing = TRUE))[1])
    
    if( cl_group != fix_tag ){
      
      ex_id <- which(cl == fix_tag)
      
      cl[which(cl!=fix_tag)] <- fix_tag
      cl[ex_id] <- fix_tag
      cent <- cent[c(2,1),]
      
    }
    
    sp_mat[ind,"cl"] <- cl
    g[ i , ]  <- as.vector(t(cent))
    
  }
  
  
  # Initialization of regression centers
  for (k in 1:center){
    
    temp_model <- list()
    for (j in 1:p) {
      
      data <- data.frame(y = g[  , j + (k-1)*p ], t = seq(1,m))
      
      sp_fit <- lm(y ~ ns(t, df = 5), data = data)
      
      g[,j + (k-1) * p ] <- sp_fit$fitted.values
      
      temp_model[[j]] <- sp_fit
      
    }
    
    model[[k]] <- temp_model
  }
  
  
  iter = 0
  g.a = g ; g.b = g + 1
  
  start_time <- Sys.time()
  
  while (sum((g.a - g.b)^2) > eps & iter < m.iter) {
    
    g.b <- g.a
    iter <- iter + 1
    # assign the class label 
    
    for (i in 1:m) {
      
      temp_id <- which(sp_mat[,"t"] == i)
      
      dis <- matrix(NA, ncol = center, nrow = length(temp_id))
      
      temp_g <- g.a[i,, drop = FALSE]
      
      for (k in 1:center) {
        
        cent_mat <- matrix(temp_g[seq(1,p) + p * (k-1)], ncol = length(temp_id), nrow = p)
        
        dis[,k] <- apply((x[,temp_id] - cent_mat)^2,2,sum)
        
      }
      
      cl <- apply(dis, 1, which.min)
      cl_group <- as.numeric(names(sort(table(cl[fix_id]), decreasing = TRUE))[1])
      
      if( cl_group != fix_tag ){
        
        ex_id <- which(cl == fix_tag)
        
        cl[which(cl!=fix_tag)] <- fix_tag
        cl[ex_id] <- fix_tag
        
      }
      
      sp_mat[temp_id ,"cl"] <- cl
      
    }
    
    for (k in 1:center){
      
      temp_model <- list()
      
      for (j in 1:p) {
        
        data <- data.frame(y = g.a[  , j + (k-1)*p ], t = seq(1,m))
        
        sp_fit <- lm(y ~ ns(t, df = 5), data = data)
        
        g.a[,j + (k-1) * p ] <- sp_fit$fitted.values
        
        temp_model[[j]] <- sp_fit
        
      }
      
      model[[k]] <- temp_model
    }
    
  }
  
  end_time <- Sys.time()
  #print(end_time - start_time )
  
  g <- g.b
  
  return(list(info=sp_mat, gamma=g, model = model))
}

MDSP_temp_clfix <- function (x, m, m.iter, eps, id = NULL, sp = "ns", center = 2,cl.t ){
  
  ##############################
  # p: p_x
  # n: number of identity
  # g: matrix for gamma
  # sp_mat: matrix stores regression and classification information
  #################################
  
  p <- nrow(x)
  n <- ncol(x)/m
  g <- matrix(0, nrow = m, ncol = center * p) 
  sp_mat <- data.frame(t = rep(seq(1,m),each =n), cl = rep(0, n *m))
  sp_mat[,"cl"] <- rep(cl.t, m)
  model <- list()
  
  for (i in 1:m) {
    
    for (k in 1:center) {
      
      temp_id <- which(sp_mat[,"t"] == i & sp_mat[,"cl"] == k)
      
      g[i,seq(1,p) + (k-1) * p ] <- apply(x[ ,temp_id, drop = FALSE], 1, mean)
      
    }
  }
  
  for (k in 1:center){
    
    temp_model <- list()
    
    for (j in 1:p) {
      
      data <- data.frame(y = g[  , j + (k-1)*p ], t = seq(1,m))
      
      sp_fit <- lm(y ~ ns(t, df = 5), data = data)
      
      g[,j + (k-1) * p ] <- sp_fit$fitted.values
      
      temp_model[[j]] <- sp_fit
      
    }
    
    model[[k]] <- temp_model
  }
  
  
  return(list(info=sp_mat, gamma=g, model = model))
  
}


#########################################

MDSP_temporal = function(Y,X,Z,m, l1,kp,id = NULL, k=2, cl.t = NULL, beta.int=NULL, beta.true , label ){
  
  #############################################
  #############  ADMMM Algorithm ##############
  #############################################
  
  
  n=length(Y)/m
  p.b = dim(X)[2]
  p.a = dim(Z)[2]
  #X.td= bdiag(list(t(KhatriRao(diag(n * m), t(X)))))
  
  ### set tunning hyperparameters and initial values
  ## hyper para
  kpp=kp # ADMM step
  lm.s=l1 # L1 regularization
  
  ## initials
  if (is.null(beta.int)){nu.int=rep(0,n*m*p.b)}else{nu.int=beta.int}
  lm.int=rep(0,n*m*p.b)
  
  lm.h.new = lm.int
  nu.h.new = nu.int
  a.h.new = rep(0,p.a)
  b.h.new = rep(0,n*p.b)
  theta.h.new = matrix(0,ncol= m * n, nrow = p.b)
  theta.h = matrix(1,ncol= m * n, nrow = p.b)
  
  iter = 1
  max.iter = 500
  eps = 0.0001
  
  error_rate = matrix(NA, nrow = m, ncol = max.iter-1)
  error_beta = rep(NA,max.iter-1)
  error = rep(NA,max.iter-1)
  
  beta_track = matrix(NA, nrow =  max.iter, ncol =n)
  nu_track = matrix(NA, nrow =  max.iter, ncol = n)
  lm_track = matrix(NA, nrow =  max.iter, ncol = n)
  
  ##### Main loop
  while(iter < max.iter & sum((theta.h-theta.h.new)^2)/(n*p.b*m) > eps){
    
    error[iter] <- sum((theta.h - theta.h.new)^2/(n*p.b*m))
    
    error_beta[iter] <- sqrt(sum((theta.h.new - beta.true)^2)/(n * p.x * m))
    
    ###
    lm.h=lm.h.new
    nu.h=nu.h.new
    b.h=b.h.new
    theta.h=theta.h.new
    
    ############## theta updates #################
    #theta.h.new = solve(t(X.td) %*% X.td + 
    #               bdiag(list(kpp*diag(n * m * p.b)))) %*% 
    #               (t(X.td) %*% Y + kpp*nu.h-lm.h)
    
    start_time <- Sys.time()
    
    t <- 1
    
    for (i in 1:m) {
      
      for (j in 1:n) {
        
        bindex <- seq(1, p.b ) + (t-1) * p.b
        
        xx <- matrix(X[t, ], nrow = 1)
        theta.h.new[,t] <- solve(t(xx) %*% xx + kpp * diag(p.b)) %*% 
          (t(xx) * Y[t] + kpp*nu.h[bindex] - lm.h[bindex])
        
        t <- t+1
        
      }
    }
    
   
    
    end_time <- Sys.time()
    #print( end_time - start_time )
    
    b.h.new = theta.h.new
    
    ######### nu,gamma updates ###############
    center = k
    
    B.h.new = matrix(b.h.new,p.b,n * m)
    Nu.h.new = matrix(nu.h,p.b,n * m)
    g.h=matrix(0,p.b, m * center)
    g.h.new=matrix(0,p.b, m * center)
    cl.h=rep(0,n * m)
    Lm=matrix(lm.h,p.b,n * m)
    
    if( is.null( cl.t ) ){
      
      sp.0 <- MDSP_temp_gamma(B.h.new, m = m, 50,0.001,id = id)
      
    }else{
      
      sp.0 <- MDSP_temp_clfix(B.h.new, m = m, 50,0.001,id = id,cl.t = cl.t )
      
    }
    
    
    g.h = sp.0$gamma - 0.01
    
    g.h.new = g.h + 1
    
    cl.h = sp.0$info
    
    iter.nu=1
    
    ### find gamma, nu
    
    while (iter.nu <50 & sum((g.h-g.h.new)^2) >0.01){
      
      b.t.new = B.h.new + kpp^{-1} * Lm
      g.h =  g.h.new
      
      for (j in 1:m) {
        
        for (i in 1:n) {
          
          temp_id <- i + (j-1) * n
          temp_b <- b.t.new[ , temp_id] 
          temp_g <-  g.h.new[j, seq(1,p.b) + (cl.h[i,"cl"]-1) * p.b]
          
          nu = temp_g + 
            sign(temp_b - temp_g) * 
            as.numeric(abs(temp_b - temp_g)>lm.s/kpp) * 
            (abs(temp_b - temp_g)-lm.s/kpp)
          
          Nu.h.new[,temp_id] <- nu
          
        }
        
      }
      
      if( is.null( cl.t ) ){
        
        sp.0 <- MDSP_temp_gamma( Nu.h.new, m = m, 50,0.001,id = id)
        
      }else{
        
        sp.0 <- MDSP_temp_clfix( Nu.h.new, m = m, 50,0.001,id = id,cl.t = cl.t )
        
      }
      
      g.h.new = sp.0$gamma
      
      cl.h = sp.0$info
      
      
      
      iter.nu <- iter.nu + 1
      
      if(iter.nu==50) {print("WARNING: maxium iteration (nu) achieves")}
      
    }
    
    
    group_alg <- cl.h
    
    for (i in 1:m) {
    
      error_rate[i,iter] <- mean( group_alg[which(  group_alg[,"t"] == i), "cl" ] == label)
    
    }
    
    ### lambda updates
    nu.h.new = as.vector(Nu.h.new)
    lm.h.new = lm.h + kpp*(b.h.new-nu.h.new)
    
    
    beta_track[iter, ] <- b.h.new
    nu_track[iter, ] <- nu.h.new
    lm_track[iter, ] <- lm.h.new
    
    iter = iter+1
    
    #print(iter)
    end_time <- Sys.time()
    #print("Whole:")
    #print(end_time - start_time)
    #print("====")
    if(iter==500) {print("WARNING: maxium iteration (nu) achieves")}
    
  }
  
  
  if( is.null( cl.t ) ){
    
    sp.h <- MDSP_temp_gamma(B.h.new, m = m, 50,0.001,id = id)
    
  }else{
    
    sp.h <- MDSP_temp_clfix(B.h.new, m = m, 50,0.001,id = id,cl.t = cl.t )
    
  }
  
  
  g.h.new = sp.h$gamma
  
  cl.h = sp.h$info
  
  group_alg <- cl.h
  
  
  sp_model = sp.h  
  
  model = sp.h$model
  
  return(list( beta.int = beta.int, Beta=B.h.new, gamma=g.h.new, group=cl.h, model = model, it.errors=error[1:(iter-1)],
               error.rate = error_rate[,1:(iter-1)], error.beta = error_beta[1:(iter-1)], 
               beta_track = beta_track[1:(iter-1),], nu_track =nu_track[1:(iter-1),],
               lm_track = lm_track[1:(iter-1),]))
  
  
}



MDSP_step = function(Y,X,Z,m, l1,kp,id = NULL, k=2, cl.t = NULL, beta.int=NULL, beta.true , label ){
  
  #############################################
  #############  ADMMM Algorithm ##############
  #############################################
  
  
  n=length(Y)/m
  p.b = dim(X)[2]
  p.a = dim(Z)[2]
  #X.td= bdiag(list(t(KhatriRao(diag(n * m), t(X)))))
  
  ### set tunning hyperparameters and initial values
  ## hyper para
  kpp=kp # ADMM step
  lm.s=l1 # L1 regularization
  
  ## initials
  if (is.null(beta.int)){nu.int=rep(0,n*m*p.b)}else{nu.int=beta.int}
  lm.int=rep(0,n*m*p.b)
  
  lm.h.new = lm.int
  nu.h.new = nu.int
  a.h.new = rep(0,p.a)
  b.h.new = rep(0,n*p.b)
  theta.h.new = matrix(0,ncol= m * n, nrow = p.b)
  theta.h = matrix(1,ncol= m * n, nrow = p.b)
  
  iter = 1
  max.iter = 500
  eps = 0.0001
  
  error_rate = matrix(NA, nrow = m, ncol = max.iter-1)
  error_beta = rep(NA,max.iter-1)
  error = rep(NA,max.iter-1)
  
  beta_track = matrix(NA, nrow =  max.iter, ncol =n)
  nu_track = matrix(NA, nrow =  max.iter, ncol = n)
  lm_track = matrix(NA, nrow =  max.iter, ncol = n)
  
  ##### Main loop
  while(iter < max.iter & sum((theta.h-theta.h.new)^2)/(n*p.b*m) > eps){
    
    error[iter] <- sum((theta.h - theta.h.new)^2/(n*p.b*m))
    
    error_beta[iter] <- sqrt(sum((theta.h.new - beta.true)^2)/(n * p.x * m))
    
    ###
    lm.h=lm.h.new
    nu.h=nu.h.new
    b.h=b.h.new
    theta.h=theta.h.new
    
    ############## theta updates #################
    #theta.h.new = solve(t(X.td) %*% X.td + 
    #               bdiag(list(kpp*diag(n * m * p.b)))) %*% 
    #               (t(X.td) %*% Y + kpp*nu.h-lm.h)
    
    start_time <- Sys.time()
    
    t <- 1
    
    for (i in 1:m) {
      
      for (j in 1:n) {
        
        bindex <- seq(1, p.b ) + (t-1) * p.b
        
        xx <- matrix(X[t, ], nrow = 1)
        theta.h.new[,t] <- solve(t(xx) %*% xx + kpp * diag(p.b)) %*% 
          (t(xx) * Y[t] + kpp*nu.h[bindex] - lm.h[bindex])
        
        t <- t+1
        
      }
    }
    
    
    
    end_time <- Sys.time()
    #print( end_time - start_time )
    
    b.h.new = theta.h.new
    
    ######### nu,gamma updates ###############
    center = k
    
    B.h.new = matrix(b.h.new,p.b,n * m)
    Nu.h.new = matrix(nu.h,p.b,n * m)
    g.h=matrix(0,p.b, m * center)
    g.h.new=matrix(0,p.b, m * center)
    cl.h=rep(0,n * m)
    Lm=matrix(lm.h,p.b,n * m)
    
    if( is.null( cl.t ) ){
      
      sp.0 <- MDSP_temp_gamma(B.h.new, m = m, 50,0.001,id = id)
      
    }else{
      
      sp.0 <- MDSP_temp_clfix(B.h.new, m = m, 50,0.001,id = id,cl.t = cl.t )
      
    }
    
    
    g.h = sp.0$gamma - 0.01
    
    g.h.new = g.h + 0.01
    
    cl.h = sp.0$info
    
    iter.nu=1
    
    ### find gamma, nu
    
    while (iter.nu <50 & sum((g.h-g.h.new)^2) >0.01){
      
      b.t.new = B.h.new + kpp^{-1} * Lm
      g.h =  g.h.new
      
      for (t in 1:m) {
        
        for (i in 1:n) {
          
          temp_id <- i + (t-1) * n
          temp_b <- b.t.new[ , temp_id] 
          temp_g <-  g.h.new[t, seq(1,p.b) + (cl.h[i,"cl"]-1) * p.b]
          
          nu = temp_g + 
            sign(temp_b - temp_g) * 
            as.numeric(abs(temp_b - temp_g)>lm.s/kpp) * 
            (abs(temp_b - temp_g)-lm.s/kpp)
          
          Nu.h.new[,temp_id] <- nu
          
        }
        
      }
      
      if( is.null( cl.t ) ){
        
        sp.0 <- MDSP_temp_gamma( Nu.h.new, m = m, 50,0.001,id = id)
        
      }else{
        
        sp.0 <- MDSP_temp_clfix( Nu.h.new, m = m, 50,0.001,id = id,cl.t = cl.t )
        
      }
      
      g.h.new = sp.0$gamma
      
      cl.h = sp.0$info
      
      
      
      iter.nu <- iter.nu + 1
      
      if(iter.nu==50) {print("WARNING: maxium iteration (nu) achieves")}
      
    }
    
    group_alg <- cl.h
    
    for (i in 1:m) {
      
      error_rate[i,iter] <- mean( group_alg[which(  group_alg[,"t"] == i), "cl" ] == label)
      
    }
    
    ### lambda updates
    nu.h.new = as.vector(Nu.h.new)
    
    step = kpp
    lm.h.new = lm.h + step*(b.h.new-nu.h.new)
    
    
    beta_track[iter, ] <- b.h.new
    nu_track[iter, ] <- nu.h.new
    lm_track[iter, ] <- lm.h.new
    
    iter = iter+1
    
    #print(iter)
    end_time <- Sys.time()
    #print("Whole:")
    #print(end_time - start_time)
    #print("====")
    if(iter==500) {print("WARNING: maxium iteration (nu) achieves")}
    
  }
  
  
  if( is.null( cl.t ) ){
    
    sp.h <- MDSP_temp_gamma(Nu.h.new, m = m, 50,0.001,id = id)
    
  }else{
    
    sp.h <- MDSP_temp_clfix(Nu.h.new, m = m, 50,0.001,id = id,cl.t = cl.t )
    
  }
  
  
  g.h.new = sp.h$gamma
  
  cl.h = sp.h$info
  
  group_alg <- cl.h
  
  
  sp_model = sp.h  
  
  model = sp.h$model
  
  return(list( beta.int = beta.int, Beta=B.h.new, gamma=g.h.new, group=cl.h, model = model, it.errors=error[1:(iter-1)],
               error.rate = error_rate[,1:(iter-1)], error.beta = error_beta[1:(iter-1)], 
               beta_track = beta_track[1:(iter-1),], nu_track =nu_track[1:(iter-1),],
               lm_track = lm_track[1:(iter-1),]))
  
  
}
