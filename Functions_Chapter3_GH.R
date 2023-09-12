simdata       <- function(n,pn,CR,delta0){
  library(Rlab)
  #n  : Sample size
  #pn : No covariates in parametric component should be pn=(300, 1000, 5000)
  #sc : Selection of scenario, sc=1; \beta2=c(0.1) sc=2; \beta2=c(0.5)
  #rho: Correlation coefficient rho=(0.7,0.9)
  
  # Generation of correlated parametric covariates X----------------------------
  # create the variance covariance matrix
  X <-  matrix(rnorm(n*pn), ncol=pn)
  #----------------------------------------------------------------------------
  #Generation of Regression coefficients---------------------------------------
  beta   <- 0
  for (j in 1:pn){
    del     <- delta0  
    beta[j] <- 3*I(j>=1 && j<=5)+ifelse(pn>5,sqrt(delta0)*I(j>5 && j<=pn),0)
  }
  #Generation of nonparametric covariate and function---------------------------
  t <- sort(runif(n,min=-0.5,max=0.5))
  #f <- -t*sin(-t^2)+2*log(t)-t^2+7   #f yeni
  f <- 1-(t^3)-2*exp(-25*(t^2))
  #-----------------------------------------------------------------------------
  #Generation of High-dimensional PLM in equation (4)---------------------------
  e <- rnorm(n,sd=0.5)
  y <- X%*%beta+f+e
  #GENERATE CENSORED PART
  z     <- 0
  delta <- 0                 #Censorship indicator (0 if data is censored and 1 otherwise)
  c     <- matrix(999,n,1)   #Censoring variable independent-identical distributed with Y
  delta <-1-rbern(n,CR)                       #censorship indicator (0):censored, (1):observed
  for (i in 1:n){
    if (delta[i]==0){
      while (y[i]<=c[i]){
        c[i]<-rnorm(1,mean(y),sd=sd(y))
      }
    }
    else{
      c[i]<-y[i]
    }
  }                           #generating censoring variable
  for (j in 1:n){
    if (y[j]<=c[j]){
      z[j]<-y[j]
    }
    else{
      z[j]<-c[j]
    }
  }              
  #-----------------------------------------------------------------------------
  out      <- new.env()
  out$X    <- X
  out$beta <- beta
  out$t    <- t
  out$f    <- f
  out$y    <- y
  out$pn   <- pn
  out$n    <- n
  out$e    <- e 
  out$z    <- z
  out$c    <- c
  out$delta <- delta
  
  return(out)
}
Smatrix       <- function(x, df,lam){
  n = length(x);
  A = matrix(0, n, n);
  for(i in 1:n){
    y = rep(0, n); y[i]=1;
    yi = smooth.spline(x, y, df=df,spar=lam)$y;
    A[,i]= yi;
  }
  return(A)
} 
remsebeta     <- function(realbetaFM,realbetaSM,FMbeta,SMbeta){
  MSEFM      <- mean((realbetaFM-FMbeta)^2)
  MSESM      <- mean((realbetaSM-SMbeta)^2)
  remse_beta <- MSEFM/MSESM
  return(remse_beta)
}
MSEf          <- function(f,estf){
  score <- (mean((f-estf)^2)/length(f))
  return(score)
} 
lasso         <- function(x,y){
  #ymean <- mean(y)
  #y <- y-mean(y)  
  #xmean <- colMeans(x)
  #xnorm <- sqrt(n-1)*apply(x,2,sd)
  #x <- scale(x, center = xmean, scale = xnorm)
  
  la_eq <- cv.glmnet(x,y,nfolds=10,intercept=F,lambda = seq(0.01,0.3,length.out=10))
  coefs <- coef(la_eq,c(la_eq$lambda.min))
  S1a  <- coefs@i+1                   #Subset of Strong Signals (index)
  S1r  <- coefs@i
  if ((length(S1a))>n){
    S1_beta <- coefs@x[1:n]
    S1 <- S1r[1:n]
  } else {
    S1_beta <- coefs@x
    S1 <- S1r
  }
  
  
  #S1_beta <- coefs@x
  X_S1 <- x[,S1]                               #Subset of Strong signals
  X_S1c <- x[,-S1]                                #Subset of complement
  las <- new.env()
  
  las$X_S1    <- X_S1
  las$X_S1c   <- X_S1c
  las$S1_beta <- S1_beta
  las$S1      <- S1
  las$coefs   <-coefs 
  
  return(las)
}
lasso_safe    <- function(x,y){
  library(glmnet)
  #ymean <- mean(y)
  #y <- y-mean(y)  
  #xmean <- colMeans(x)
  #xnorm <- sqrt(n-1)*apply(x,2,sd)
  #x <- scale(x, center = xmean, scale = xnorm)
  
  la_eq <- cv.glmnet(x,y,nfolds=10,intercept=F,lambda = seq(0.0005,0.0075,length.out=10))
  coefs <- coef(la_eq,c(la_eq$lambda.min))
  S1a  <- coefs@i+1                             #Subset of Strong Signals (index)
  S1r  <- coefs@i
  if ((length(S1a))>n){
    S1_beta <- coefs@x[1:n]
    S1 <- S1r[1:n]
  } else {
    S1_beta <- coefs@x
    S1 <- S1r
  }
  
  #S1_beta <- coefs@x
  X_S1 <- x[,S1]                               #Subset of Strong signals
  X_S1c <- x[,-S1]                             #Subset of complement
  #-------------------------------------------------------------------------------
  #SAFE RULE----------------------------------------------------------------------
  maxlam <- 0
  ind    <- NA
  lam2   <- la_eq$lambda.min
  for (j in 1:ncol(X_S1)){
    maxlam[j] <- t(X_S1[,j])%*%y 
    if (maxlam[j]<abs(2*lam2-maxlam[j])){
      ind[j] <- j
    }
  }
  indnew    <- ind[!is.na(ind)]
  X_SAFE    <- X_S1[,indnew] 
  beta_SAFE <- S1_beta[indnew] 
  #-------------------------------------------------------------------------------  
  las <- new.env()
  las$X_S1    <- X_S1
  las$X_S1c   <- X_S1c
  las$S1_beta <- S1_beta
  las$S1      <- S1
  las$coefs   <- coefs 
  las$X_SAFE  <- X_SAFE
  las$beta_SAFE <- beta_SAFE
  las$indnew  <-indnew 
  return(las)
}
adalasso_safe <- function(x,y){
  #fit ols 
  beta.init <- solve(t(x)%*%x,tol=1e-100,0.5*diag(ncol(x)))%*%t(x)%*%y # exclude 0 intercept
  
  # calculate weights
  w  <- abs(beta.init)  
  x2 <- scale(x, center=FALSE, scale=1/w)  
  
  la_eq <- cv.glmnet(x2,y,nfolds=10,intercept=F,lambda = seq(0.005,0.045,length.out=10))
  #la_eq <- cv.glmnet(x2,y,nfolds=10,intercept=F)
  coefs <- coef(la_eq,c(la_eq$lambda.min))
  S1a  <- coefs@i+1                              #Subset of Strong Signals (index)
  S1r  <- coefs@i
  if ((length(S1a))>n){
    S1_beta <- coefs@x[1:n]
    S1 <- S1r[1:n]
  } else {
    S1_beta <- coefs@x
    S1 <- S1r
  }
  #S1_beta <- coefs@x
  X_S1 <- x[,c(S1)]                              #Subset of Strong signals
  X_S1c <- x[,-S1]                                #Subset of complement
  #-------------------------------------------------------------------------------
  #SAFE RULE----------------------------------------------------------------------
  maxlam <- 0
  ind    <- NA
  lam2   <- la_eq$lambda.min
  for (j in 1:ncol(X_S1)){
    maxlam[j] <- t(X_S1[,j])%*%y 
    if (maxlam[j]<abs(2*lam2-maxlam[j])){
      ind[j] <- j
    }
  }
  indnew <- ind[!is.na(ind)]
  X_SAFE <- X_S1[,indnew] 
  beta_SAFE <- solve(t(X_SAFE)%*%X_SAFE)%*%t(X_SAFE)%*%y 
  #-------------------------------------------------------------------------------  
  alas <- new.env()
  alas$X_S1    <- X_S1
  alas$X_S1c   <- X_S1c
  alas$S1_beta <- S1_beta
  alas$S1      <- S1
  alas$coefs   <- coefs
  alas$X_SAFE  <- X_SAFE 
  alas$beta_SAFE <- beta_SAFE
  alas$indnew   <-indnew 
  return(alas)
}
aicfunc       <- function(modelvar,lam,betahat){
  pn      <- length(betahat)
  score   <- (modelvar)+1+((2*pn+1)/(n-pn-2))
  return(score)
}
zerone        <- function(x){
  xnew <- (x-min(x))/(max(x)-min(x))
}
syndata       <- function(y,delta){
  library(pracma)
  # where y: right-censored observations, delta: censorship indicator
  n<-length(y)
  M<-0
  yg<-0
  M<-0
  delta2<-0
  #Synthetic data transformation
  y1<-cbind(y,delta)
  y2<-sortrows(y1,k=2)
  delta2<-y2[,2]
  delta2[n]<-1
  sy<-sort(y)
  for (i1 in 1:n){
    Ma<-1
    for (i2 in 1:n){
      mGt<-((n-i2)/(n-i2+1))^((delta2[i2]==0)*(sy[i2]<=sy[i1]))
      Ma<-Ma*mGt
    }
    M[i1]=Ma
    yg[i1]<-(delta[i1]*y[i1])/M[i1]
  }
  return(yg)
}
shrinkage     <- function(beta_FM,beta_SM,T,p2){
  beta_s <- beta_FM-(beta_FM-beta_SM)*(1-(p2-2)*(1/T[1]))
  return(beta_s)
}
pretest       <- function(beta_FM,beta_SM,T,p2){
  beta_PT     <- beta_FM-(beta_FM-beta_SM)*I(T[1]<=qchisq(p=0.05,p2)) 
  return(beta_PT)
} 