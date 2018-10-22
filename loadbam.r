#assert existence of
#commonfile
source(commonfile)
#bamfile
load(bamfile)
#pwmfile
#load(paste0(pwmdir,pwmid,'.pwmout.RData'))
#tmpdir

bamnames = names(allreads[[1]])

coords2=sapply(coords.short,flank,width=wsize,both=T)

obschrnames=names(allreads)
preads=allreads[[obschrnames[1]]][[1]]$plus
cutat=10

#stablize the variance (helps when there are few high coverage sites).
tfun <- function(x){
    y = x
    x[x>cutat]=cutat
    y[x>0] = sqrt(x[x>0])
    y
}

unlink(paste0(tmpdir,'/',pwmid,'/','*tf',pwmid,'*'))

use.w=!is.null(whitelist)
if(use.w){
    wtable=read.table(whitelist)
    white.list=lapply(levels(wtable[,1]),function(i){
        wtchr=wtable[wtable[,1]==i,]
        ir=IRanges(wtchr[,2],wtchr[,3])
    })
    names(white.list)=levels(wtable[,1])
}

dir.create(paste0(tmpdir,'/',pwmid),recursive=T)
makeTFmatrix <- function(coords,prefix='',offset=0){
    cwidth = width(coords[[1]][1])
    obschrnames=names(allreads)
    validchr = obschrnames[which(obschrnames%in%ncoords)]
    readcov=lapply(1:length(bamnames),function(j){
        sapply(validchr,function(i){length(allreads[[i]][[j]]$plus)+length(allreads[[i]][[j]]$minus)})/seqlengths(genome)[validchr]
    })
    readfact = lapply(readcov,function(i){i/i[1]})
    slen = seqlengths(genome)
    scrd =sapply(coords,length)
    minbgs=floor(max(10000,sum(scrd))*(scrd/sum(scrd)));
    for(chr in validchr){
        chrlen = slen[chr]
        print(chr)
        if(prefix=='background.'){
            nsites = max(length(coords[[chr]]),minbgs[chr])
            coind = sample(start(coords[[chr]]),nsites,replace=T)+offset
            if(use.w){
                wchr=white.list[[chr]]
                wlarge=wchr[width(wchr)>(2*wsize+1)]
                csamp = sample(1:length(wlarge),nsites,prob=(width(wlarge)-(2*wsize)),replace=T)
                starts = start(wlarge)[csamp]
                ends = end(wlarge)[csamp]
                coind=floor((ends-starts - 2*wsize) * runif(length(csamp)))+(starts+wsize)
            }
            chrcoord=sort(IRanges(start=coind-wsize,width=2*wsize))
        }else{
            chrcoord=coords[[chr]]
        }
        pos.mat = do.call(rBind,lapply(1:length(bamnames),function(i){
            pluscoord=allreads[[chr]][[i]]$plus
            if(length(pluscoord)>0){
                rre = rle(sort(pluscoord))
                irp=IRanges(start=rre$values,width=1)
                fos=findOverlaps(chrcoord,irp)
                uquery=queryHits(fos)
                querycoord=rre$values[subjectHits(fos)]
                uoffset = querycoord-start(chrcoord)[uquery]+1
                rval= rre$lengths[subjectHits(fos)] / readfact[[i]][chr]
                pos.triple = cbind(round(uquery),round(uoffset),tfun(rval))
                pos.mat=sparseMatrix(i=round(uoffset),j=round(uquery),x=tfun(rval),dims=c(2*wsize,length(chrcoord)),giveCsparse=T)
            }else{
                pos.triple=cbind(1,1,0)
                pos.mat=Matrix(0,nrow=2*wsize,ncol=length(chrcoord))
            }
            pos.mat
        }))
    #
        neg.mat = do.call(rBind,lapply(1:length(bamnames),function(i){
            minuscoord=allreads[[chr]][[i]]$minus
            if(length(minuscoord)>0){
                rre = rle(sort(minuscoord))
                irp=IRanges(start=rre$values,width=1)
                fos=findOverlaps(chrcoord,irp)
                uquery=queryHits(fos)
                querycoord=rre$values[subjectHits(fos)]
                uoffset = querycoord-start(chrcoord)[uquery]+1
                rval= rre$lengths[subjectHits(fos)] / readfact[[i]][chr]
                neg.triple = cbind(round(uquery),round(uoffset),tfun(rval))
                neg.mat=sparseMatrix(i=round(uoffset),j=round(uquery),x=tfun(rval),dims=c(2*wsize,length(chrcoord)),giveCsparse=T)
            }else{
                neg.triple=cbind(1,1,0)
                neg.mat=Matrix(0,nrow=2*wsize,ncol=length(chrcoord))
            }
            neg.mat
        }))
#
        save(pos.mat,neg.mat,file=paste0(tmpdir,'/',pwmid,'/',prefix,'tf',pwmid,'-',chr,'.RData'))
	gc()
    }
}

makeTFmatrix(coords2,'positive.')
makeTFmatrix(coords2,'background.',10000)

#
#####
