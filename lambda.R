lamoptimal<-function(x,y,n,varth,gam,al,ph,alphatrue)
{
	lamU=0.1
	lamL=1
	lams=seq(lamL,lamU,length=20)
	ss=length(lams)
	BIC1=seq(-1,-1,length=ss)
	BIC2=seq(-1,-1,length=ss)
	GCV=seq(-1,-1,length=ss)
	deg=seq(-1,-1,length=ss)
      	Qnn=seq(-1,-1,length=ss)
	for(i in 1:ss){
		lam=lams[i]
		result=estimation(x,y,n,lam,varth,gam,al,ph,alphatrue)
		muhat=result[[2]]
            K=result[[1]]
            betahat=result[[3]]
            etahat=result[[4]]
            group=result[[5]]
            delta=result[[6]]
            df=K+ncol(x)
           # dff=degree(x,y,n,lam,varth,gam,al,ph,K,muhat,betahat,etahat,group,delta)
           # df=dff[[1]]
            deg[i]=df
            # Phat=matrix(0,nrow=n,ncol=K)
            # group=group[,1:K]
            #group=as.matrix(group)
            # for(i in 1:K){
            # gi=group[,i]
            # gi=gi[which(gi!=0)]
            #  Phat[gi,i]=1
           #   }
           # Pr=cbind(Phat,x)
           #  thetahat=solve(t(Pr)%*%Pr)%*%t(Pr)%*%y
           #  muhat=thetahat[1:ncol(Phat)]
           #  betahat=thetahat[(1+ncol(Phat)):ncol(Pr)]
           #  muhat=Phat%*%muhat
            Qn=t(y-muhat-x%*%betahat)%*%(y-muhat-x%*%betahat)/n
            Qn=Qn[1,1]
            Qnn[i]=Qn
            yhat=muhat+x%*%betahat
		BIC1[i]=log(Qn)+10*log(log(n+ncol(x)))*log(n)*df/n
            BIC2[i]=log(Qn)+log(n)*df/n
            GCV[i]=Qn/(1-df/n)^2

	}
	
	lam=lams[min(which(BIC1==min(BIC1)))]

       #lam=lams[min(which(GCV==min(GCV)))]
	
	    return(list(lam,BIC1))
  }