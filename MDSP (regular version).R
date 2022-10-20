require(mixtools)
############## packages
require(Matrix)

#load MDSP_gamma function
MDSP_gamma <- function (x, m.iter, eps){
  
  #x: numeric vector
  #m.iter: maximum number of iterations
  #eps: tolerance level for stopping rule
  
  #g.int=3*mean(x)
  
  if(length(unique(x))<2)
  {
    class.g=rep(1,length(x)); g.b=mean(x); print("homogeneous!!! in MDSP Gamma")
    
  }else{
    
         k.int=kmeans(x, centers=2)
         g.int=k.int$centers[which.max(abs(k.int$centers))]*0.95
         iter=1
         g.a=g.int; g.b=g.int+1
         
         while (abs(g.b-g.a)>eps & iter < m.iter)
         {
           
           class.g=as.numeric(abs(x-g.a)< abs(x))
           g.b=sum(class.g * x)/sum(class.g)
           iter=iter+1
           g.a=g.b
           
         }
         
         if(iter==m.iter) {print("WARNING: maxium iteration achieves")}
  }
  
  return(list(l.group=class.g, gamma=g.b))
  
}

#l.group: group indicator (non-zero subgroup =1) 
# gamma: sub-homogeneous effect of non-zero subgroup



### Modified BIC function


#########################################

MDSP <- function(Y,X,Z,l1,kp,beta.int=NULL){
  
    #############################################
    #############  ADMMM Algorithm ##############
    #############################################
    
    ## Input data X (heterogeneous, n by p.b), 
    ## Z(homogeneous including intercept, n by p.a),  and 
    ## Y (n by 1)
  
    n=length(Y)
    p.b=dim(X)[2]
    p.a=dim(Z)[2]
    X.td=cbind(Z, KhatriRao(diag(n), t(X)))
    
    ### set tunning hyperparameters and initial values
    ## hyper para
    ## kpp for $\| \beta - \mu \|^2$
    ## l1 for $lambda_{N,m}$
    
    kpp=kp #ADMM step
    lm.s=l1 #L1 regularization
    
    ## initials
    if (is.null(beta.int)){
      
      nu.int=rep(0,n)
      
    }else{
      
      nu.int=beta.int
        
    }
    
    #nu.int=beta.int
    #nu.int=rep(0,n)
    lm.int=rep(0,n)
    
    lm.h.new=lm.int
    nu.h.new=nu.int
    
    a.h.new=rep(0,p.a)
    b.h.new=rep(0,n*p.b)
    
    theta.h=rep(0,n*p.b+p.a)
    theta.h.new=rep(1,n*p.b+p.a)
    
    iter=1
    max.iter=100
    eps=0.001
    
    error=rep(1,max.iter-1)
    
    
    #NOTE: It's better to have a warm start
    
    ptm <- proc.time()
    
    ##### Main loop
    while(iter < max.iter & sqrt(mean((theta.h-theta.h.new)^2)) > eps)
    {
      error[iter]=sqrt(mean((theta.h-theta.h.new)^2))
      
      ###
      lm.h = lm.h.new
      nu.h = nu.h.new
      
      a.h = a.h.new
      b.h = b.h.new
      theta.h=theta.h.new
      
      ############## theta updates##############
      theta.h.new = solve(t(X.td) %*% X.td + bdiag(list(matrix(0,p.a,p.a),kpp*diag(n)))) %*% 
        (t(X.td) %*% Y + c(rep(0,p.a),kpp*nu.h-lm.h))
      a.h.new=theta.h.new[1:p.a]
      b.h.new=theta.h.new[-(1:p.a)]
      
      ######### nu,gamma updates ###############
      sp.0 = MDSP_gamma(b.h.new, 50,0.001)
      g.h = sp.0$gamma - 0.01
      g.h.new = g.h + 0.01
      cl.h=sp.0$l.group
      iter.nu=1
      
      ### find gamma, nu
      while (iter.nu <50 & abs(g.h-g.h.new)>0.001)
      {
        g.h = g.h.new  
        b.t.new = b.h.new+kpp^{-1}*lm.h
        #lm.s.ad=lm.s*min(iter/5,1)
        lm.s.ad = lm.s
        
        nu.0 = sign( b.t.new ) * as.numeric( abs( b.t.new ) > lm.s.ad/kpp ) * ( abs( b.t.new ) - lm.s.ad/kpp )
        nu.g = g.h + sign( b.t.new - g.h ) * as.numeric(abs( b.t.new - g.h ) > lm.s.ad/kpp ) * ( abs(b.t.new-g.h) - lm.s.ad/kpp )
        nu.h.new = nu.g*cl.h+nu.0*(1-cl.h)
        
        sp.h = MDSP_gamma(nu.h.new, 50,0.001)
        g.h.new = sp.h$gamma
        cl.h = sp.h$l.group
        iter.nu = iter.nu+1
        #if(iter.nu==50) {print("WARNING: maxium iteration (nu) achieves")}
      }
      
      
      ### lambda updates
      lm.h.new=lm.h+kpp*(b.h.new-nu.h.new)
      
      iter=iter+1
      if(iter==max.iter) {print("WARNING: maxium iteration achieves")}
    }
    
    return(list(alpha.h=a.h.new, beta.h=b.h.new, nu.h=nu.h.new,gamma.h=g.h.new, group=cl.h, it.errors=error[1:(iter-1)]))
    
    print(proc.time()-ptm)

}


#################################################################################

# 
# # #############Data generation
 n=100
 p.x=1
 p.z=2
 
 set.seed(1122)
 ### parameter
 gamma=4
 id.beta=rep(c(0,1), each=n/2)
 id.beta=rep(1,n)
 beta=gamma*id.beta
 alpha=rep(1,p.z)
 
# 
# ### data
 X=matrix(rnorm(n*p.x),n,p.x)
 Z=cbind(1, matrix(rnorm(n*(p.z-1)),n,(p.z-1)))
 #X.td=cbind(Z, KhatriRao(diag(n), t(X)))
 X.td=cbind(Z, bdiag(as.list(X)))
 Y=X.td%*%c(alpha,beta)+rnorm(n)
# 
 fit.m=MDSP(Y,X,Z,l1=0.2,kp=1)
 fit.m1=MDSP(Y,X,Z,l1=0.2,kp=1,beta.int=round(fit.m$beta.h,2))
 fit.m2=MDSP(Y,X,Z,l1=1,kp=1,beta.int=round(fit.m1$beta.h,2))
 par(mfrow=c(1,3))
 plot(fit.m$it.errors)
 plot(fit.m1$it.errors)
 plot(fit.m2$it.errors)
# 
# 
 plot(fit.m$nu.h)
 plot(fit.m1$nu.h)
 plot(fit.m2$nu.h)
# 
# 
# 
#  plot(fit.m$beta.h)
#  points(beta,col="red", pch=17)
#  abline(h=fit.m$gamma.h,col="blue")
#  fit.m$gamma.h
#  mean(fit.m$group==as.numeric(beta!=0))
# 
#  plot(fit.m1$beta.h)
#  points(beta,col="red", pch=17)
#  abline(h=fit.m1$gamma.h,col="blue")
#  fit.m1$gamma.h
#  mean(fit.m1$group==as.numeric(beta!=0))
# 
# 
#  plot(fit.m2$beta.h)
#  points(beta,col="red", pch=17)
#  abline(h=fit.m2$gamma.h,col="blue")
#  fit.m2$gamma.h
#  mean(fit.m2$group==as.numeric(beta!=0))
#  
#  
# 
#  m.fit=regmixEM(as.numeric(Y), cbind(X,Z), k = 2, addintercept = FALSE)
#  m.group=as.numeric(m.fit$posterior[,1]<0.5)
#  1-mean(m.group==as.numeric(beta!=0))

