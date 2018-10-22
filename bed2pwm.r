#!/usr/bin/Rscript
#Rscript script common-path jaspar-dir id output
#example:
#Rscript pwmmatch.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/pwms/jaspar.txt 141 /cluster/thashim/basepiq/tmp/pwmout.RData

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

commonfile = args[1]
bedin = args[2]
bedname = args[3]
bedid = as.double(args[4])
outdir = args[5]

pwmname = bedname
pwmid = bedid

outdir=paste0(outdir,'/')
source(commonfile)
if(!overwrite & file.exists(paste0(outdir,pwmid,'.pwmout.RData'))){
  stop("pwm file already exists")
}



####
# load PWMs
####

require('GenomicRanges')
rt = read.table(bedin)

gr = GRanges(rt[,1],ranges=IRanges(start=rt[,2],end=rt[,3]),strand=rt[,6],score=rt[,5])
chrstr = unique(format(rt[,1],justify='none'))

gstr = gr[strand(gr)=='+']

coords=lapply(chrstr,function(i){
    gsub=gstr[seqnames(gstr)==i]
    ranges(gsub)
})
coords.pwm=lapply(chrstr,function(i){
    gsub=gstr[seqnames(gstr)==i]
    score(gsub)
})

clengths=sapply(coords,length)
print(sum(clengths))
coords.short=coords[clengths>0]
names(coords.short)=chrstr[clengths>0]
ncoords=chrstr[clengths>0]#names(coords)
coords2=sapply(coords.short,flank,width=wsize,both=T)

pwmin=matrix(rep(1,16),4,4)
ipr=matrix(rep(0,16),4,4)

save(coords,coords.pwm,ipr,pwmin,pwmname,chrstr,clengths,coords.short,ncoords,coords2,file=paste0(outdir,pwmid,'.pwmout.RData'))


gstr = gr[strand(gr)=='-']

coords=lapply(chrstr,function(i){
    gsub=gstr[seqnames(gstr)==i]
    ranges(gsub)
})
coords.pwm=lapply(chrstr,function(i){
    gsub=gstr[seqnames(gstr)==i]
    score(gsub)
})

clengths=sapply(coords,length)
print(sum(clengths))
coords.short=coords[clengths>0]
names(coords.short)=chrstr[clengths>0]
ncoords=chrstr[clengths>0]#names(coords)
coords2=sapply(coords.short,flank,width=wsize,both=T)

pwmin=matrix(rep(1,16),4,4)
ipr=matrix(rep(0,16),4,4)

save(coords,coords.pwm,ipr,pwmin,pwmname,chrstr,clengths,coords.short,ncoords,coords2,file=paste0(outdir,pwmid,'.pwmout.rc.RData'))
