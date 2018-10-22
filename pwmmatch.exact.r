#!/usr/bin/Rscript
#Rscript script common-path jaspar-dir id output
#example:
#Rscript pwmmatch.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/pwms/jaspar.txt 141 /cluster/thashim/basepiq/tmp/pwmout.RData

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

commonfile = args[1]
jaspardir = args[2]
pwmid = as.double(args[3])
outdir = args[4]

outdir=paste0(outdir,'/')
source(commonfile)
if(!overwrite & file.exists(paste0(outdir,pwmid,'.pwmout.RData'))){
  stop("pwm file already exists")
}



####
# load PWMs
####

#pwmin = 'pwms/'


importJaspar <- function(file=myloc) {
  vec <- readLines(file)
  vec <- gsub("\t"," ",vec)
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

pwmnum = pwmid
pwmin = pwmtable[[pwmnum]] + 1e-20
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

#chr names
chrstr = seqnames(genome)

if(exists('blacklist') & !is.null(blacklist)){
    blacktable=read.table(blacklist)
}

if(exists('whitelist') & !is.null(whitelist)){
    whitetable=read.table(whitelist)
}

#####
# fw motif match

pwuse = ipr

coords.list = lapply(chrstr,function(i){
    print(i)
    gi=genome[[i]]
    if(remove.repeatmask & !is.null(masks(gi))){
        active(masks(gi)) <- rep(T,length(masks(gi)))
    }
    if(exists('blacklist') & !is.null(blacklist)){
        blacksel= blacktable[,1]==i
        if(sum(blacksel)>0){
            flsize = wsize*flank.blacklist
            ir=intersect(IRanges(1,length(gi)),reduce(IRanges(blacktable[blacksel,2]-flsize,blacktable[blacksel,3]+flsize)))
            mask=Mask(length(gi),start(ir),end(ir))
            if(is.null(masks(gi)))
                masks(gi) = mask
            else
                masks(gi) = append(masks(gi),mask)
        }
    }
    if(exists('whitetable')){
        whitesel=whitetable[,1]==i
        if(sum(whitesel)>0){
            wchr=whitetable[whitesel,,drop=F]
            ir=IRanges(wchr[,2],wchr[,3])
            air=IRanges(1,length(gi))
            nir=setdiff(air,ir)
            rir=reduce(IRanges(start(nir)-wsize,end(nir)+wsize))
            maskr=intersect(rir,air)
            mask = Mask(length(gi),start(maskr),end(maskr))
            if(is.null(masks(gi)))
                masks(gi) = mask
            else
                masks(gi) = append(masks(gi),mask)
        }else{
            mask = Mask(length(gi),1,length(gi))
            if(is.null(masks(gi)))
                masks(gi) = mask
            else
                masks(gi) = append(masks(gi),mask)
        }
    }
    mpwm=matchPWM(pwuse,gi,min.score=motifcut)
    pscore=PWMscoreStartingAt(pwuse,as(gi,"DNAString"),start(mpwm))
    list(mpwm,pscore)
})

if(sum(sapply(coords.list,function(i){length(i[[2]])}))>0){

allpwm=do.call(c,lapply(coords.list,function(i){i[[2]]}))
pwmcut2=sort(allpwm,decreasing=T)[min(length(allpwm),maxcand)]
rm(allpwm)
print(pwmcut2)

coords=lapply(1:length(coords.list),function(i){
    as(coords.list[[i]][[1]],'IRanges')[coords.list[[i]][[2]] >= pwmcut2]
})

coords.pwm=lapply(coords.list,function(i){i[[2]][i[[2]] >= pwmcut2]})

#coords=lapply(coords.list,unlist)

clengths=sapply(coords,length)
print(sum(clengths))
coords.short=coords[clengths>0]
names(coords.short)=chrstr[clengths>0]
ncoords=chrstr[clengths>0]#names(coords)
coords2=sapply(coords.short,flank,width=wsize,both=T)

save(coords,coords.pwm,ipr,pwmin,pwmname,chrstr,clengths,coords.short,ncoords,coords2,file=paste0(outdir,pwmid,'.pwmout.RData'))

}else{
clengths=0
save(clengths,file=paste0(outdir,pwmid,'.pwmout.RData'))
}

#
#####

#####
# RC motif match

pwuse = reverseComplement(ipr)

coords.list = lapply(chrstr,function(i){
    print(i)
    gi=genome[[i]]
    if(remove.repeatmask & !is.null(masks(gi))){
        active(masks(gi)) <- rep(T,length(masks(gi)))
    }
    if(exists('blacklist') & !is.null(blacklist)){
        blacksel= blacktable[,1]==i
        if(sum(blacksel)>0){
            flsize = wsize*flank.blacklist
            ir=intersect(IRanges(1,length(gi)),reduce(IRanges(blacktable[blacksel,2]-flsize,blacktable[blacksel,3]+flsize)))
            mask=Mask(length(gi),start(ir),end(ir))
            if(is.null(masks(gi)))
                masks(gi) = mask
            else
                masks(gi) = append(masks(gi),mask)
        }
    }
    if(exists('whitetable')){
        whitesel=whitetable[,1]==i
        if(sum(whitesel)>0){
            wchr=whitetable[whitesel,,drop=F]
            ir=IRanges(wchr[,2],wchr[,3])
            air=IRanges(1,length(gi))
            nir=setdiff(air,ir)
            rir=reduce(IRanges(start(nir)-wsize,end(nir)+wsize))
            maskr=intersect(rir,air)
            mask = Mask(length(gi),start(maskr),end(maskr))
            if(is.null(masks(gi)))
                masks(gi) = mask
            else
                masks(gi) = append(masks(gi),mask)
        }else{
            mask = Mask(length(gi),1,length(gi))
            if(is.null(masks(gi)))
                masks(gi) = mask
            else
                masks(gi) = append(masks(gi),mask)
        }
    }
    mpwm=matchPWM(pwuse,gi,min.score=motifcut)
    pscore=PWMscoreStartingAt(pwuse,as(gi,"DNAString"),start(mpwm))
    list(mpwm,pscore)
})

if(sum(sapply(coords.list,function(i){length(i[[2]])}))>0){

allpwm=do.call(c,lapply(coords.list,function(i){i[[2]]}))
pwmcut2=sort(allpwm,decreasing=T)[min(length(allpwm),maxcand)]
rm(allpwm)
print(pwmcut2)

coords=lapply(1:length(coords.list),function(i){
    as(coords.list[[i]][[1]],'IRanges')[coords.list[[i]][[2]] >= pwmcut2]
})

coords.pwm=lapply(coords.list,function(i){i[[2]][i[[2]] >= pwmcut2]})

#coords=lapply(coords.list,unlist)

clengths=sapply(coords,length)
print(sum(clengths))
coords.short=coords[clengths>0]
names(coords.short)=chrstr[clengths>0]
ncoords=chrstr[clengths>0]#names(coords)
coords2=sapply(coords.short,flank,width=wsize,both=T)

save(coords,coords.pwm,ipr,pwmin,pwmname,chrstr,clengths,coords.short,ncoords,coords2,file=paste0(outdir,pwmid,'.pwmout.rc.RData'))

}else{
clengths=0
save(clengths,file=paste0(outdir,pwmid,'.pwmout.rc.RData'))
}
