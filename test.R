rm(list = ls())
t1 <- Sys.time()
setwd("C:/Users/Administrator/Desktop/huangjian")
#library("MSBVAR")
library("Matrix")
library("MASS")
source("rmultnorm.R")
source("lambda.R")
source("estimation.R")
nsim=1
n=100

#生成X的参数
dx=5
mux<-rep(0,dx)
rho=0.3
sigx=rho*matrix(1,nrow=dx,ncol=dx)+(1-rho)*diag(dx)
sigx=sigx

#空矩阵
Ksim=seq(0,0,length=nsim)
MSEmu=seq(0,0,length=nsim)
MSEbeta=seq(0,0,length=nsim)
MSEmuor=seq(0,0,length=nsim)
MSEbetaor=seq(0,0,length=nsim)
betasim=matrix(0,nrow=nsim,ncol=dx)
alphasim=matrix(0,nrow=nsim,ncol=n)
alphasig=matrix(0,nrow=nsim,ncol=n)
betasig=matrix(0,nrow=nsim,ncol=dx)
pvaluesim=matrix(NA,nrow=nsim,ncol=n)
alphasigor=matrix(0,nrow=nsim,ncol=2)
alphasimor=matrix(0,nrow=nsim,ncol=2)
betasigor=matrix(0,nrow=nsim,ncol=dx)
lamall = seq(0,0,length=nsim)
set.seed(1)

beta=runif(dx,min=0.5,max=1)

#大循环
for(i in 1:nsim){
  
  e=rnorm(n,mean=0,sd=1)#误差
  x=rmultnorm(n, mux, sigx)
  u=runif(n)
  
  #各样本的mu取值
  alp=2
  #mu=(u<(1/4))*2*(-alp)+(u>=(1/4))*(u<(1/2))*(-alp)+(u>=(1/2))*(u<(3/4))*alp+(u>=(3/4))*2*alp
  #mu=(u<(1/3))*alp+(u>(2/3))*(-alp)
  mu=(u<=(1/2))*alp+(u>(1/2))*(-alp)
  
  #真实的组标签
  z1=as.numeric(u<=(1/2))
  z2=as.numeric(u>(1/2))
  z=cbind(z1,z2)
  
  #Example3
  #mu=2
  #z=matrix(1,nrow=n,ncol=1)
  
  alphatrue=c(-alp,alp)
  #alphatrue=mu
  K=length(alphatrue)

  
  y=mu+x%*%beta+e
  
  #Q=cbind(as.numeric(mu==alp),as.numeric(mu!=alp))
  #QX=cbind(Q,x)
  #theta=solve(t(QX)%*%QX)%*%t(QX)%*%y
  #muhat=theta[1:2]
  #muhat=Q%*%muhat
  ##betahat=theta[3:length(theta)]
  #MSEmu[i]=sqrt(t(muhat-mu)%*%(muhat-mu)/n)
  #MSEbeta[i]=sqrt(t(betahat-beta)%*%(betahat-beta)/ncol(x))
  
  al=4 # 1 LASSO  2 MCP  3 SCAD 4 truncated L1
  varth=1
  gam=3
  ph=1
  tau=0.5 # truncated parameter for truncated L1
  

  lamopt=lamoptimal(x,y,n,varth,gam,al,ph,alphatrue)
  lam=lamopt[[1]]
  lamall[i]=lam
  
  
  #开始估计参数
  result=estimation(x,y,n,lam,varth,gam,al,ph,alphatrue)
  muhat=result[[2]]
  betahat=result[[3]]
  Khat=result[[1]]
  Ksim[i]=Khat
  MSEmu[i]=sqrt(t(muhat-mu)%*%(muhat-mu)/n)
  MSEbeta[i]=sqrt(t(betahat-beta)%*%(betahat-beta)/ncol(x))
  alphahat=result[[7]]
  zhat=result[[8]]
  sighat=result[[9]]
  betasim[i,]=betahat
  alphasim[i,1:Khat]=alphahat
  des=cbind(zhat,x)
  deINV=solve(t(des)%*%des)
  sign=sighat*sqrt(diag(deINV))
  alphasig[i,1:Khat]=sign[1:Khat]
  betasig[i,]=sign[(Khat+1):length(sign)]
  an=seq(0,0,length=length(sign))
  
  if(Khat>1){
    for(j in 2:Khat){
      an[1]=1
      an[j]=-1
      sigj=sighat*sqrt(t(an)%*%deINV%*%an)
      zcrit=abs(alphahat[1]-alphahat[j])/sigj
      zcrit=zcrit[1]
      pvaluesim[i,(j-1)]=2*(1-pnorm(zcrit))
    }
  }
  MOD=lm(y~z+x-1)
  ORcoef=MOD$coef
  alphaor=ORcoef[1:ncol(z)]
  alphasimor[i,]=alphaor
  muor=z%*%alphaor
  betaor=ORcoef[(ncol(z)+1):(ncol(z)+ncol(x))]
  MSEmuor[i]=sqrt(t(muor-mu)%*%(muor-mu)/n)
  MSEbetaor[i]=sqrt(t(betaor-beta)%*%(betaor-beta)/ncol(x))
  sigor=sqrt(t(y-z%*%alphaor-x%*%betaor)%*%(y-z%*%alphaor-x%*%betaor)/(n-K-dx))
  des=cbind(z,x)
  deINV=solve(t(des)%*%des)
  sign=sigor*sqrt(diag(deINV))
  alphasigor[i,1:K]=sign[1:K]
  betasigor[i,]=sign[(K+1):length(sign)]
  
  #i为迭代次数
  print(i)
}