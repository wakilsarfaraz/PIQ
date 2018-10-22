#assert existence of
#commonfile
source(commonfile)
#pwmfile
#load(file.path(pwmdir,paste0(pwmid,'.pwmout.RData')))
#tmpdir
#outdir

#####
# make bg
datadir = paste0(tmpdir,'/',pwmid,'/')

nsites = Reduce('+',lapply(list.files(datadir,paste0('positive.tf',pwmid,'-'),full.names=T),function(i){load(i);ncol(pos.mat)}))
if(nsites < 100){fast.mode=F}

load(file.path(datadir,paste0('background.tf',pwmid,'-',seqnames(genome)[1],'.RData')))

if( two.pass & file.exists(paste0(outdir,'/params/twopass.RData')) ){
    load(paste0(outdir,'/params/twopass.RData'))
    tmat.pos = cov.premade[1:(2*wsize),]
    tmat.neg = cov.premade[(2*wsize+1):(4*wsize),]
}

#use precomputed tables if wsize is 1000
if( (!fast.mode) & (wsize <= 1000)){
    load('precalc/cov.premade.RData')
    subsel = 1000+(-(wsize-1):wsize)
    tmat.pos = cov.premade[subsel,]
    tmat.neg = cov.premade[subsel + (2000),]
}

#require at least 30k sites to do a background estimate
if((ncol(pos.mat) > 30000) & (!fast.mode) & (wsize != 1000)){
source('covfun.r')
require(Matrix)
lf=list.files(datadir,paste0('background.tf',pwmid,'-'),full.names=T)
bgwsize=2*wsize
pmc=Reduce('+',lapply(lf,function(i){
    print(i)
    load(i)
    m=rBind(pos.mat[1:bgwsize,],neg.mat[1:bgwsize,])
    m%*%t(m)
}))
alldat=do.call(cBind,lapply(lf,function(i){
    load(i)
    rBind(pos.mat[1:bgwsize,],neg.mat[1:bgwsize,])
}))
datin=t(as.matrix(rBind(pos.mat[1:bgwsize,],neg.mat[1:bgwsize,])))
datin[datin>3]=3
stabval=4
mustab=-3
sz=ncol(datin)
bind=makeblocks(1:bgwsize,2)
tabs = makeTable(datin,mustab,stabval)
matin=matrix(tabs[[1]][as.double(datin)+1],nrow(datin))
covbase = (pmc-min(pmc))/(pmc[1,1]-min(pmc))
covin = covbase
sci=solve(covin)
muin= mustab
mupri=rep(muin,sz)%*%sci
scid = solve(sci+diag(rep(tabs[[2]][1],sz)))
for(m in 1:3){
    print(m)
    musum=rep(0,ncol(covin))
    devsum=0
    muevsum=0
    covsumT=matrix(0,ncol(covin),ncol(covin))
    sqmuT=matrix(0,ncol(covin),ncol(covin))
    muT=rep(0,ncol(covin))
    tct=0
    numex=min(nrow(datin),5000)
    for(i in sample(1:nrow(datin),numex)){
        print(i)
        ff=fastFit(datin[i,],mupri,tabs[[1]],tabs[[2]],scid)
        muall= t(ff[[1]])
        covsum = ff[[2]]
        devs = 1/(1/diag(ff[[2]])+tabs[[2]][datin[i,]+1])
        muevs= (ff[[1]]/diag(ff[[2]])-tabs[[1]][datin[i,]+1]*tabs[[2]][datin[i,]+1])*devs
        #
        sqmuT=sqmuT+t(muall)%*%muall
        muT=muT+colSums(muall)
        covsumT=covsumT+covsum
        devsum=devsum+sum(devs)
        muevsum=muevsum+sum(muevs)
        tct=tct+nrow(muall)
    }
    cch=sqmuT/tct-(t(t(muT))%*%t(muT))/(tct^2)
    cpost=(covsumT/tct+cch)
    muhat=muT/tct
    print(c(stabval,muin))
    tabs = makeTable(datin,muin,stabval)
    print(tabs)
    covin=cpost
    sci=solve(covin)
    mupri=muhat%*%sci
    scid = solve(sci+diag(rep(tabs[[2]][1],sz)))
}
svpos=svd(covin[1:(2*wsize),])
svneg=svd(covin[2*wsize+(1:(2*wsize)),])

ecsel = 1:ncol(svpos$u)
ecpv = svpos$d/max(svpos$d)*100
ecnv = svneg$d/max(svneg$d)*100
ecpv[ecpv < 1]=1
ptrans = t(t(svpos$u[,ecsel])*(sqrt(ecpv[ecsel])))
ecnv[ecnv < 1]=1
ntrans = t(t(svneg$u[,ecsel])*(sqrt(ecnv[ecsel])))
tmat.pos = cbind(ptrans, matrix(0,2*wsize,length(ecsel)))
tmat.neg = cbind(matrix(0,2*wsize,length(ecsel)), ntrans)
}else{
	rot.pos = diag(wsize*2)
	rot.neg = diag(wsize*2)
}

#####
# make svlite

converter.src <- '
Rcpp::NumericMatrix triplePos(rtriplePos);
Rcpp::NumericMatrix tripleNeg(rtripleNeg);
Rcpp::NumericVector sums(rsum);
int label = as<int>(rlabel);
int negoff = as<int>(rnegoff);
Rcpp::CharacterVector cv(sums.size());
int tripleindPos = 0;
int tripleindNeg = 0;
int len = 65536 * 10;
for(int i=0; i < sums.size(); i++){
  char str[len];
  double sumct = sums[i];
  int cp = snprintf(str, 65536, "%d %d:%f", label, 1, sumct);
  while(tripleindPos < triplePos.nrow() && triplePos(tripleindPos,1)-1 == i){
    int crd = triplePos(tripleindPos,0)+1;
    double wt = triplePos(tripleindPos,2);
    cp += snprintf(str + cp, len-cp, " %d:%f", crd, wt);
    tripleindPos++;
  }
  while(tripleindNeg < tripleNeg.nrow() && tripleNeg(tripleindNeg,1)-1 == i){
    int crd = tripleNeg(tripleindNeg,0)+1+negoff;
    double wt = tripleNeg(tripleindNeg,2);
    cp += snprintf(str + cp, len-cp, " %d:%f", crd, wt);
    tripleindNeg++;
  }
  cv[i] = str;
}
return cv;
'
converter <- cxxfunction(signature(rtriplePos="numeric",rtripleNeg="numeric",rsum="numeric",rnegoff="integer",rlabel="integer"),converter.src,plugin="Rcpp",includes="#include <stdio.h>")

sumtr <-function(x){
    x
}

makesvlite <- function(filename,label,rot.pos,rot.neg,minread=5){
    print(filename)
    load(filename)
    rct=sumtr(colSums(pos.mat)+colSums(neg.mat))+1
    pms=as.matrix(summary(pos.mat))
    pns=as.matrix(summary(neg.mat))
    cv=converter(pms,pns,rct,nrow(pos.mat),label)
    cv[rct>=(sumtr(minread)+1)]
}

if(fast.mode){

validpos = list.files(datadir,paste0('positive.tf',pwmid,'-'),full.names=T)
validneg = list.files(datadir,paste0('background.tf',pwmid,'-'),full.names=T)

allpos=do.call(c,lapply(validpos,makesvlite,label=1,minread=max(10*wsize/1000,1),rot.pos=rot.pos,rot.neg=rot.neg))
allneg=do.call(c,lapply(validneg,makesvlite,label= -1,minread= 0, rot.pos=rot.pos,rot.neg=rot.neg))
writeLines(c(allpos,allneg),file.path(tmpdir,paste0(pwmid,'-svlite.txt')))

itmax=min(5e+07,500*length(allpos))

wsz = 2*nrow(pos.mat)+2
sv.fit=sofia(file.path(tmpdir,paste0(pwmid,'-svlite.txt')),verbose=T,dimensionality=wsz,random_seed=1,lambda=1*10^6/(length(allpos)*wsz),iterations=itmax,learner_type='logreg-pegasos',eta_type='basic',loop_type='balanced-stochastic')

vpos=suppressMessages(as.vector((sv.fit$weights[2+(1:(nrow(pos.mat)))])))
vneg=suppressMessages(as.vector((sv.fit$weights[2+(1:(nrow(pos.mat)))+(nrow(pos.mat))])))
sv.rotate = c(sv.fit$weights[1:2],vpos,vneg)

save(sv.fit,sv.rotate,file=file.path(tmpdir,paste0(pwmid,'.svout.RData')))


}else{

validpos = list.files(datadir,paste0('positive.tf',pwmid,'-'),full.names=T)
validneg = list.files(datadir,paste0('background.tf',pwmid,'-'),full.names=T)

allnom = c('label','rct',paste0('x',1:ncol(tmat.pos)))
alldat=rbind(
    do.call(rbind,lapply(validpos,function(i){
        print(i)
        load(i)
        rct = sumtr(colSums(pos.mat) + colSums(neg.mat))+1
        y=data.frame(label=rep(1,sum(rct>(sumtr(10)+1))),rct[rct>(sumtr(10)+1)],x=t(as.matrix(t(tmat.pos)%*%pos.mat[,rct>(sumtr(10)+1)] + t(tmat.neg)%*%neg.mat[,rct>(sumtr(10)+1)])))
        colnames(y)=allnom
        y
    })),
    do.call(rbind,lapply(validneg,function(i){
        print(i)
        load(i)
        rct = sumtr(colSums(pos.mat) + colSums(neg.mat))+1
        y=data.frame(label=-1,rct,x=t(as.matrix(t(tmat.pos)%*%pos.mat + t(tmat.neg)%*%neg.mat)))
        colnames(y)=allnom
        y
})))

write.svmlight(alldat[,1],as.matrix(alldat[,-1]),file.path(tmpdir,paste0(pwmid,'-svlite.txt')))

itmax=min(5e+07,1000*sum(alldat[,1]==1))

sv.fit=sofia(file.path(tmpdir,paste0(pwmid,'-svlite.txt')),verbose=T,dimensionality=ncol(alldat),random_seed=1,lambda=1000/sum(alldat[,1]==1),iterations=itmax,learner_type='logreg-pegasos',eta_type='basic',loop_type='balanced-stochastic')

sv.rotate=c(sv.fit$weights[1:2],tmat.pos%*%(sv.fit$weight[-(1:2)]),tmat.neg%*%(sv.fit$weight[-(1:2)]))

save(sv.fit,sv.rotate,file=file.path(tmpdir,paste0(pwmid,'.svout.RData')))

}


#
#####
