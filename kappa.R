
MDSP_kappa = function(Y,X,Z,m, l1,kp,id = NULL, k=2, cl.t = NULL, beta.int=NULL, label ){
  
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
  kappa_track = rep(NA,max.iter-1)
  
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
    
    s = -kpp * (nu.h.new - nu.h)
    r = (b.h.new-nu.h.new)
    
    kappa_track[iter] = kpp
      
    ### penalty parameter updating
    if(sqrt(sum(r^2)) > 10 * sqrt(sum(s^2)) ){
      
      kpp = kpp * 2
      
    }
    
    if(sqrt(sum(s^2)) > 10 * sqrt(sum(r^2))){
      
      kpp = kpp / 2
      
    }
    
    
    lm.h.new = lm.h + kpp*r
    
    
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
               lm_track = lm_track[1:(iter-1),]),
         kappa_track = kappa_track)
  
  
}
