#FUNCTIONS
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
#-------------------------------------------------------------------------------
x <- read.table("X.txt")
y <- read.table("Y.txt")
delta0 <- read.table("cens.txt")
attach(y)
y<-(y$V1)
y <- zerone(y)
y <-ifelse(y>0.6,y=mean(y),y)
n<-length(y)
attach(delta0)
delta0<-delta0$V1
delta <- 1-delta0
xr<-data.matrix(x)
x<-zerone(matrix(c(xr[,1:444],xr[,446:549]),115,548))
x <- zerone(x)
t<-zerone(xr[,393]+rnorm(n,sd=0.05))
p<-ncol(x)
#-------------------------------------------------------------------------------                                #Number of simulation
del <- 0                   #Delta to violate sparsity
yg  <-syndata(y,delta)     #synthetic data transformation
################################################################################
l_seq  <- 15
lamseq <- seq(0.7,0.8,length.out=l_seq)
AICc   <- 0
for (j in 1:l_seq){
  Smatp       <- Smatrix(t,df=2,lamseq[j])
  xtilp       <- (diag(n)-Smatp)%*%x
  ytilp       <- (diag(n)-Smatp)%*%yg
  obj_Lassop  <- lasso_safe(xtilp,ytilp)
  xsp         <- obj_Lassop$X_SAFE 
  Hp          <- xsp%*%solve(t(xsp)%*%xsp,tol=1e-100)%*%t(xsp)
  var_Mod     <- (t(yg-Hp%*%yg)%*%(yg-Hp%*%yg))/(n-ncol(xsp)) 
  AICc[j]     <- log(var_Mod)+1+((2*(ncol(xsp)+1))/(n-ncol(xsp)-2))
}
for (j2 in 1:l_seq){
  if (min(AICc)==AICc[j2]){
    opt.lam1 <- lamseq[j2]
  }
}

#-------------------------------------------------------------------------------
Smat       <- Smatrix(t,df=2,opt.lam1) 
xtil       <- (diag(n)-Smat)%*%x
ytil       <- (diag(n)-Smat)%*%yg
obj_Lasso  <- lasso_safe(xtil,ytil)
obj_aLasso <- adalasso_safe(xtil,ytil)
#-------------------------------------------------------------------------------
#FUll Model Estimation----------------------------------------------------------
beta_FM <- obj_Lasso$beta_SAFE
xtil_FM <- obj_Lasso$X_SAFE
x_FM_UN <- x[,obj_Lasso$S1]             #UNSAFE
x_FM    <- x_FM_UN[,obj_Lasso$indnew]
S1_FM   <- obj_Lasso$S1[obj_Lasso$indnew] 
#-------------------------------------------------------------------------------
beta_SM <- obj_aLasso$beta_SAFE
xtil_SM <- obj_aLasso$X_SAFE
xtil_S2 <- obj_aLasso$X_S1c
x_SM_UN <- x[,obj_aLasso$S1]             #UNSAFE
x_SM    <- x_SM_UN[,obj_aLasso$indnew]
p1      <- length(beta_SM)
S1_SM   <- obj_aLasso$S1[obj_aLasso$indnew]
#-------------------------------------------------------------------------------
#zeros     <- rep(0,(length(S1_FM)-length(S1_SM)))
p2        <- (length(S1_FM)-length(S1_SM)) 
#beta_SM.F <- c(beta_SM,zeros)
beta_SM.F <- 0
ss <- 1
for(i in 1:length(S1_FM)){
  if(S1_FM[i]%in% S1_SM){
    beta_SM.F[i] <- beta_SM[ss]
    ss <- ss+1
  }
  else{
    beta_SM.F[i] <- 0
  }
}
#Necessary Arguments------------------------------------------------------------
yhat_SM <- x_SM%*%beta_SM+(Smat%*%(yg-x_SM%*%beta_SM))
resid_SM<- yg-yhat_SM 
var_SM  <- (1/abs(n-p1))*(t(resid_SM))%*%(resid_SM) 
U       <- diag(n)-xtil_SM%*%solve(t(xtil_SM)%*%xtil_SM)%*%t(xtil_SM)
beta2_LS<- solve(t(xtil_S2)%*%U%*%xtil_S2+0.1*diag(ncol(xtil_S2)))%*%t(xtil_S2)%*%U%*%ytil
T       <- (1/var_SM)*t(beta2_LS)%*%solve(t(xtil_S2)%*%U%*%xtil_S2+0.1*diag(ncol(xtil_S2)))%*%beta2_LS*(p/p1)
#-------------------------------------------------------------------------------
beta_s  <- shrinkage(beta_FM,beta_SM.F,T,p2) 
beta_pt <- pretest(beta_FM,beta_SM.F,T,p2)
#-------------------------------------------------------------------------------
yhat_FM <- x_FM%*%beta_FM+(Smat%*%(yg-x_FM%*%beta_FM))
yhat_S  <- x_FM%*%beta_s+(Smat%*%(yg-x_FM%*%beta_s))
yhat_PT <- x_FM%*%beta_pt+(Smat%*%(yg-x_FM%*%beta_pt))

MSEFM   <- mean((yg-yhat_FM)^2)  
MSESM   <- mean((yg-yhat_SM)^2)
MSES   <- mean((yg-yhat_S)^2)  
MSEPT   <- mean((yg-yhat_PT)^2)

varFM   <- (1/abs(n-ncol(x_FM)))*sum((yg-yhat_FM)^2)  
varSM   <- (1/abs(n-p1))*sum((yg-yhat_SM)^2)
varS    <- (1/abs(n-p1))*sum((yg-yhat_S)^2)  
varPT   <- (1/abs(n-p1))*sum((yg-yhat_PT)^2)

data.frame(varS,varPT,varSM,varFM)

ReMSE_SM <- MSEFM/MSESM 
ReMSE_S  <- MSEFM/MSES 
ReMSE_PT <- MSEFM/MSEPT 
#-------------------------------------------------------------------------------
fhat_S  <- Smat%*%(yg-x_FM%*%beta_s)
fhat_PT <- Smat%*%(yg-x_FM%*%beta_pt)  
fhat_FM <- Smat%*%(yg-x_FM%*%beta_FM)
fhat_SM <- Smat%*%(yg-x_SM%*%beta_SM)

fp_S  <- yg-x_FM%*%beta_s
fp_PT <- yg-x_FM%*%beta_pt
fp_SM <- yg-x_SM%*%beta_SM
fp_FM <- yg-x_FM%*%beta_FM

MSEf_S  <- MSEf(fp_S,fhat_S)
MSEf_PT <- MSEf(fp_PT,fhat_PT)  
MSEf_SM <- MSEf(fp_SM,fhat_SM)
MSEf_FM <- MSEf(fp_FM,fhat_FM)

data.frame(MSEf_S,MSEf_PT,MSEf_SM,MSEf_FM)

data.frame(ReMSE_S,ReMSE_PT,ReMSE_SM)

data.frame(MSES,MSEPT,MSESM,MSEFM)
BP <- (c(0.167,0.163,0.217))

par(mar=c(4,4,2,1))
par(mfrow=c(2,1))
plot(sort(t),zerone(fhat_S),type="l",xlab="393rd Gene",ylab="f(393rd)",lwd=2,main="(a) Fitted Curves for NSBC Dataset")
lines(sort(t),zerone(fhat_PT),col="chartreuse4",lty=2,lwd=2)
lines(sort(t),zerone(fhat_SM),col="coral",lty=3,lwd=2)
par(new=TRUE)
plot(t,zerone(yg),col="gray",pch=19,cex=0.55,xlab="393rd Gene",ylab="f(393rd)")
par(new=TRUE)
plot(t,zerone(y),pch=19,cex=0.55,xlab="393rd Gene",ylab="f(393rd)")
legend("topright",legend=c("f(S)","f(PT)","f(SM)","yg","z"),col=c(1,"chartreuse4","coral","gray",1),lty=c(1,2,3,NA,NA),pch=c(NA,NA,NA,19,19),cex=0.75)

library(RColorBrewer)
coul <- brewer.pal(5, "Set2") 
mybar <- barplot(BP,names=c("f(S))","f(PT)","f(SM)"),ylab="MSE(f)",xlab="Method",col=coul,ylim=c(0,0.35),main="(b) Barplots for MSE scores")
text(mybar,BP+0.02,paste("MSE :",BP,sep=""),cex=1)
grid()
legend("topleft",legend=c("f(S)","f(PT)","f(SM)"),col=coul,bty="n",pch=20,pt.cex =2,cex=0.8,inset=c(0.02,0.001))

