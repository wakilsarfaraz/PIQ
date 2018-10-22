
gknotLG = gauss.quad(9000,kind="hermite")
poiscompLG <- function(mu,sig,dat){
  t=mu+sqrt(2*sig)*gknotLG$nodes
  dp=dpois(dat,exp(t))*gknotLG$weights/sqrt(pi)
  evs=dp%*%cbind(1,t,t^2)
  un=evs[2]/evs[1]
  vn=evs[3]/evs[1]-evs[2]^2/evs[1]^2
  c(un,vn,evs[1])
}

getLapLoc<-function(mu,sig,d,cut=20){
  l=rep(0,length(d))
  l[d<cut]=d[d<cut]*sig[d<cut]+mu[d<cut]-lambert_W0(exp((d*sig+mu)[d<cut])*sig[d<cut])
  l[d>cut]=(log(d[d>cut])*d[d>cut]+mu[d>cut]/sig[d>cut])/(d[d>cut]+1/sig[d>cut])
  sh=1/(exp(l)+1/sig)
  pr=sqrt(2*pi*sh)*dpois(d,exp(l))*dnorm(l,mu,sqrt(sig))
  list(l,sh,pr)
}

makeTable <- function(data,muoffs,stabval,maxcut=100){
  maxval=max(data)
  llmu=sapply(0:maxval,function(i){
    if(i<maxcut){
      v=poiscompLG(muoffs,stabval,i)
      (v[1]/v[2]-muoffs/stabval)/(1/v[2]-1/stabval)
    }else{
      (log(i)*i+muoffs/stabval)/(i+1/stabval)
    }
  })
  llprec=sapply(0:maxval,function(i){
    if(i<maxcut){
      1/poiscompLG(muoffs,stabval,i)[2]-1/stabval
    }else{
      me=(log(i)*i+muoffs/stabval)/(i+1/stabval)
      exp(me)-1/stabval
    }
  })
  list(llmu,llprec)
}

fastFit <- function(data,mupri,muadjusts,sigadjusts,scid){
  datind=which(data>0)
  if(length(datind)>0){
  sigadjust=sigadjusts[data[datind]+1]-sigadjusts[1]
  #
  dss=scid[datind,datind,drop=F]*outer(sqrt(sigadjust),sqrt(sigadjust),'*')
  sccm=t(scid[datind,,drop=F])
  cminternal=(solve(dss+diag(length(datind))))*outer(sqrt(sigadjust),sqrt(sigadjust),'*')
  a1=sccm%*%cminternal%*%t(sccm)
  allinv=(scid-a1)
  vh=as.vector(muadjusts[data+1]*sigadjusts[data+1]+mupri)
  uinv=scid%*%vh-sccm%*%(cminternal%*%(t(sccm)%*%vh))
  list(uinv,allinv)
  }else{
    list(scid%*%as.vector(mupri+muadjusts[1]*sigadjusts[1]),scid)
  }
}

fastMu <- function(data,mupri,muadjusts,sigadjusts,scid){
  datind=which(data>0)
  if(length(datind)>0){
  sigadjust=sigadjusts[data[datind]+1]-sigadjusts[1]
  #
  dss=scid[datind,datind,drop=F]*outer(sqrt(sigadjust),sqrt(sigadjust),'*')
  sccm=t(scid[datind,,drop=F])
  cminternal=(solve(dss+diag(length(datind))))*outer(sqrt(sigadjust),sqrt(sigadjust),'*')
  #a1=sccm%*%cminternal%*%t(sccm)
  #allinv=(scid-a1)
  sapply(mupri,function(mu){
    vh=as.vector(muadjusts[data+1]*sigadjusts[data+1]+mu)
    uinv=scid%*%vh-sccm%*%(cminternal%*%(t(sccm)%*%vh))
  })
  }else{
    sapply(mupri,function(mu){
      scid%*%as.vector(mu+muadjusts[1]*sigadjusts[1])
    })
  }
}

#cin = solve(covpost+zerovar)
calcPR <- function(data,mupris,muadjusts,sigadjusts,cin){
  datind=which(data>0)
  dv= sapply(mupris,function(mu){
      as.vector((muadjusts[data+1]-mu)%*%cin%*%(muadjusts[data+1]-mu))
  })
  if(length(datind)>0){
  sigadjust=1/(sigadjusts[data[datind]+1]-sigadjusts[1])
  #
  dss=cin[datind,datind,drop=F]*outer(sqrt(sigadjust),sqrt(sigadjust),'*')
  sccm=t(cin[datind,,drop=F])
  cminternal=(solve(dss+diag(length(datind))))*outer(sqrt(sigadjust),sqrt(sigadjust),'*')
  da=sapply(mupris,function(mu){
    sum((((muadjusts[data+1]-mu)%*%sccm)^2)%*%cminternal)
    #%*%(t(sccm)%*%(muadjusts[data+1]-mu))
    #as.vector((muadjusts[data+1]-mu)%*%cin%*%(muadjusts[data+1]-mu))
  })
  #a1=sccm%*%cminternal%*%t(sccm)
  #allinv=(cin-a1)
  dv-da
  }else{
    dv
  }
}

makeblocks <- function(xs,bln){
  bind=do.call(c,lapply(1:bln,function(i){
    do.call(c,lapply(1:bln,function(j){
      ofx=length(xs)*(i-1)
      ofy=length(xs)*(j-1)
      c(list(cbind(xs+ofx,xs+ofy)),
        lapply(1:(length(xs)-1),function(i){
          cbind(ofx+xs[-(length(xs):(length(xs)-i+1))],ofy+xs[-(1:i)])
        }))
    }))
  }))
}



toepvals <- function(x,bind){
  sapply(bind,function(i){mean(x[i])})
}

toeptomat <- function(x,bind,sz){
  cmt=matrix(0,sz,sz)
  for(i in 1:length(bind)){
    cmt[bind[[i]]]=x[i]
  }
  cmt[cmt==0]=t(cmt)[cmt==0]
  ind1=1:(sz/2)
  ind2=ind1+(sz/2)
  same.strand.avg=(cmt[ind1,ind1]+cmt[ind2,ind2])/2
  cmt[ind1,ind1]=same.strand.avg
  cmt[ind2,ind2]=same.strand.avg
  cmt
}

topproj <- function(x,bind){
  spp=toepvals(x,bind)
  toeptomat(spp,bind,ncol(x))
}

topiter <- function(matinp,bind){
  matin=topproj(matinp,bind)
  ecp=eigen(matin)
  ecv=ecp$values
  while(min(ecv)< -1e-2){
    print(min(ecv))
    matin=ecp$vectors%*%diag(ecv*(ecv>0))%*%t(ecp$vectors)
    matin=topproj(matin,bind)
    ecp=eigen(matin)
    ecv=ecp$values
  }
  matin+diag(ncol(matin))*0.011
}
