#!/usr/bin/Rscript
#Rscript script commonfile bamout
#example:
#Rscript bam2rdata.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/tmp/bams.RData /cluster/cwo/dnase_seq/bams/D0_50-100_130801.bwa.mapq20.mm10.bam

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

commonfile = args[1]
bamout = args[2]

source(commonfile)

if(exists('readGAlignmentPairs')){
    readBam = readGAlignmentPairs
}else{
    readBam = readBamGappedAlignmentPairs
}

#nucleosome peaks at 210 bp, 400bp
fragrange = c(0,100,150,250,350,1000)
fragnames = paste(fragrange[-length(fragrange)],fragrange[-1],sep='-')

bamname = args[3]
bamnames=bamname

flags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isDuplicate=F,isNotPassingQualityControls=F),what=c('mapq'))
reads = readBam(bamname,param=flags)
subs = strand(first(reads))=='+'
lreads= first(reads)[subs]
rreads= last(reads)[subs]
fraglen = end(rreads)-start(lreads)
lstart = start(lreads)
rstart = end(rreads)

obschrnames=levels(c(seqnames(reads)))
allreads=lapply(obschrnames,function(chr){
	print(chr)
        lx = lapply(1:(length(fragrange)-1),function(fr){
            qsel = (mcols(lreads)$mapq > mapq)
            if(any(is.na(qsel))) qsel = T
            fsel = fraglen > fragrange[fr] & fraglen <= fragrange[fr+1]
            select = (seqnames(lreads)==chr) & qsel & fsel
            pluscoord=lstart[which(select)]
            qsel = (mcols(rreads)$mapq > mapq)
            if(any(is.na(qsel))) qsel = T
            select = (seqnames(rreads)==chr) & qsel & fsel
            minuscoord=rstart[which(select)]
            list(plus=pluscoord,minus=minuscoord)
        })
        names(lx)=fragnames
        lx
})

names(allreads)=obschrnames
save(allreads,file=bamout)
