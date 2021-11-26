estimation<-function(x,y,n,lam,varth,gam,al,ph,alphatrue)
{
dx=ncol(x)
#muini=y
#betaini=seq(0,0,length=dx)
Inx=cbind(1,x)
thetain=solve(t(Inx)%*%Inx)%*%t(Inx)%*%y
#muini=rep(1,length=n)*thetain[1]
betaini=thetain[2:length(thetain)]
#betaini=solve(t(x)%*%x)%*%t(x)%*%y
muini=y-x%*%betaini

etaini=seq(0,0,length=(n*(n-1)/2))
#etatrue=seq(0,0,length=(n*(n-1)/2))
delta=matrix(0,nrow=n,ncol=(n*(n-1)/2))
e=diag(n)
for (i in 1:(n-1)){
ind=(i-1)*n-i*(i-1)/2
etaini[(ind+1):(ind+n-i)]=muini[i]-muini[(i+1):n]
#etatrue[(ind+1):(ind+n-i)]=mu[i]-mu[(i+1):n]
delta[,(ind+1):(ind+n-i)]=e[,i]-e[,(i+1):n]
}
delta=t(delta)
upsini=seq(0,0,length=(n*(n-1)/2))
weight=exp(-ph*etaini^2)
#weight=1/abs(etaini)^2

muold=muini
betaold=betaini
etaold=etaini
upsold=upsini
varthold=varth
roold=0
Px=solve(t(x)%*%x)%*%t(x)
Qx=x%*%Px
epsilon=10^(-3)#-------------收敛的阈值
step=0
rnew=2
snew=2
rr=NULL
ss=NULL
mud=NULL
while(rnew>epsilon)
{

step=step+1
rold=delta%*%muold-etaold
rold=sqrt(t(rold)%*%rold)
#Iinv=(1+n*varthold)*diag(n)-varthold*rep(1,length=n)%*%t(rep(1,length=n))-Qx
Iinv=diag(n)+n*varthold*t(delta)%*%delta-Qx#jia n
Iinv=solve(Iinv)
II=(diag(n)-Qx)%*%y+n*varthold*t(delta)%*%(etaold-(1/varthold)*upsold)#==jia n
munew=Iinv%*%II
betanew=Px%*%(y-munew)
del=delta%*%munew+(1/varthold)*upsold
lam1=lam/varthold
eta=abs(del)-lam1
eta=eta*(eta>0)
eta=sign(del)*eta

if(al==1){
lamL=lam1*weight
etaL=abs(del)-lamL
etaL=etaL*(etaL>0)
etanew=sign(del)*etaL
#etanew=eta
}else if(al==2){
etanew=(eta/(1-(varthold*gam)^(-1)))*(abs(del)<=(gam*lam))+del*(abs(del)>(gam*lam)) 
}else if(al==3){
eta1=eta*(abs(del)<=(lam+lam1))
eta2=del*(abs(del)>(gam*lam))
lam2=(lam/varthold)*(gam/(gam-1))
etat=abs(del)-lam2
etat=etat*(etat>0)
etat=sign(del)*etat 
eta3=(etat/(1-((gam-1)*varthold)^(-1)))*(abs(del)<=(gam*lam))*(abs(del)>(lam+lam1))
etanew=eta1+eta2+eta3  #SCAD
}else{
etaL=abs(del)-lam1
etaL=etaL*(etaL>0)
dell=etaL*(del/abs(del))
etanew=del*(abs(etaold)>=tau)+dell*(abs(etaold)<tau)
}

upsnew=upsold+varthold*(delta%*%munew-etanew)
varthnew=(1*varthold)*(rold>(0.25*roold))+varthold*(rold<=(0.25*roold))
varthnew=varthnew[1,1]
#varthnew=varthold

roold=rold
rnew=delta%*%munew-etanew
rnew=sqrt(t(rnew)%*%rnew)
rold=rnew
snew=t(delta)%*%(etanew-etaold)
snew=sqrt(t(snew)%*%snew)
mudiff=sqrt(t(muold-munew)%*%(muold-munew))

muold=munew
betaold=betanew
etaold=etanew
upsold=upsnew
varthold=varthnew
rr=c(rr,rnew)
ss=c(ss,snew)
mud=c(mud,mudiff)
}

seq=1:n
group=matrix(0,nrow=n,ncol=n)
alphaold=seq(0,0,length=n)
zhat=matrix(0,nrow=n,ncol=n)
K=1

while(length(seq)>0)
{
i=seq[1]
ind=(i-1)*n-i*(i-1)/2
id=which(etaold[(ind+1):(ind+n-i)]==0)
id=c(i,i+id)
id <- intersect(1:n,id)#=================
group[1:length(id),K]=id
alphaold[K]=mean(muold[id])
zhat[id,K]=1
seq=seq [! seq %in% id]
K=K+1
}

K=K-1
alphaold=alphaold[1:K]
ind=seq(1:length(alphaold))
zhat=zhat[,1:K]
zhat=as.matrix(zhat)
iz=seq(0,0,length=length(alphatrue))
if(length(alphaold)>1){
for(ii in 1:length(alphatrue)){
	diff=abs(alphaold-alphatrue[ii])
	com=cbind(diff,ind)
	com=com[sort.list(com[,1]), ]	
	iz[ii]=com[1,2]
}	
iz <- unique(iz)#=======================
alphaold1=alphaold[iz]
alphaold2=alphaold[-iz]
alphaold=c(alphaold1,alphaold2)
zhat1=zhat[,iz]
zhat2=zhat[,-iz]
zhat=cbind(zhat1,zhat2)
   }
sig=sqrt(t(y-zhat%*%alphaold-x%*%betaold)%*%(y-zhat%*%alphaold-x%*%betaold)/(n-K-dx))


return(list(K,muold,betaold,etaold,group,delta,alphaold,zhat, sig))
  }