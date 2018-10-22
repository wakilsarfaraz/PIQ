#!/usr/bin/Rscript
#Rscript script common-path jaspar-dir id output
#example:
#Rscript pwmmatch.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/pwms/jaspar.txt 141 /cluster/thashim/basepiq/tmp/pwmout.RData

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

commondir = args[1]
jaspardir = args[2]
id = as.double(args[3])
outdir = args[4]

if(file.exists(paste0(outdir,id,'.pwmout.RData'))){
  stop("pwm file already exists")
}

source(commondir)

####
# load PWMs
####

#pwmin = 'pwms/'


importJaspar <- function(file=myloc) {
  vec <- readLines(file)
  vec <- gsub("\\[|\\]", "", vec)
  start <- grep(">", vec); end <- grep(">", vec) - 1
  pos <- data.frame(start=start, end=c(end[-1], length(vec)))
  pwm <- sapply(seq(along=pos[,1]), function(x) vec[pos[x,1]:pos[x,2]])
  pwm <- sapply(seq(along=pwm), function(x) strsplit(pwm[[x]], " {1,}"))
  pwm <- lapply(seq(along=start), function(x) matrix(as.numeric(t(as.data.frame(pwm[(pos[x,1]+1):pos[x,2]]))[,-1]), nrow=4, dimnames=list(c("A", "C", "G", "T"), NULL)))
  names(pwm) <- gsub(">", "", vec[start])
  return(pwm)
}
pwmtable <- importJaspar(jaspardir)

pwmnum = id
pwmin = pwmtable[[pwmnum]]
pwmname = names(pwmtable)[pwmnum]

####
# end input script
# assert: existence of pwmin and pwmname
####


####
# motif match

pwmnorm=t(t(pwmin)/colSums(pwmin))
#informbase=colSums((log(pwmnorm+0.01)-log(1/4))*pwmnorm) #
#pwmnorm = pwmnorm[,(informbase > basecut)]
ipr=log(pwmnorm)-log(1/4)
pwmgen=(t(t(pwmin)/colSums(pwmin)))^(0.7)

#generate nkmer motif match candidates
str=apply(pwmgen,2,function(j){
    sample(1:4,nkmer,prob=j,replace=T)
})

scores=rowSums(sapply(1:ncol(str),function(i){
    ipr[str[,i],i]
}))

#filter by motif match cutoff
passcut= which(scores > motifcut)
print(length(passcut))

basenames=rownames(pwmin)
strs=sapply(passcut,function(i){
    paste0(basenames[str[i,]],collapse='')
})
ustrs=unique(strs)
uscores=as.double(scores[passcut][match(ustrs,strs)])
if(match.rc){
  ustrs=as.character(reverseComplement(DNAStringSet(ustrs)))
}

if(length(ustrs)>0){

#calculate chr offsets
chrstr = seqnames(genome)

#find motif matches
pd=PDict(ustrs)
coords.list=lapply(chrstr,function(i){
    print(i)
    mpd=matchPDict(pd,genome[[i]])
})

coords.pwm=sapply(coords.list,function(i){
    ci=countIndex(i)
    cid=which(ci>0)
    do.call(c,lapply(cid,function(j){
        rep(uscores[j],ci[j])
    }))
})

allpwm=do.call(c,coords.pwm)
pwmcut2=sort(allpwm,decreasing=T)[min(length(allpwm),maxcand)]
rm(allpwm)
print(pwmcut2)

coords=lapply(1:length(coords.list),function(i){
    unlist(coords.list[[i]])[coords.pwm[[i]] > pwmcut2]
})

coords.pwm=lapply(coords.pwm,function(i){i[i>pwmcut2]})

#coords=lapply(coords.list,unlist)

clengths=sapply(coords,length)
coords.short=coords[clengths>0]
names(coords.short)=chrstr[clengths>0]
ncoords=chrstr[clengths>0]#names(coords)
coords2=sapply(coords.short,flank,width=wsize,both=T)


save(coords,coords.pwm,ipr,pwmin,pwmname,ustrs,chrstr,clengths,coords.short,ncoords,coords2,file=paste0(outdir,id,'.pwmout.RData'))

}else{
clengths=0
save(clengths,file=paste0(outdir,id,'.pwmout.RData'))
}

#
#####
