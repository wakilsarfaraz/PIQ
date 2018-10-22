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

if(exists('readGAlignments')){
    readBam = readGAlignments
}else{
    readBam = readBamGappedAlignments
}

#if only one bam..
if(length(args)==3){
bamname = args[3]
bamnames=bamname

plusflags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isMinusStrand=F,isDuplicate=F,isNotPassingQualityControls=F),what=c('mapq'))
plusstrand = readBam(bamname,param=plusflags)

minusflags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isMinusStrand=T,isDuplicate=F,isNotPassingQualityControls=F),what=c('mapq'))
minusstrand = readBam(bamname,param=minusflags)

obschrnames=levels(c(seqnames(plusstrand),seqnames(minusstrand)))
allreads=lapply(obschrnames,function(chr){
	print(chr)
        qsel = (mcols(plusstrand)$mapq > mapq)
        if(any(is.na(qsel))) qsel = T
        select = (seqnames(plusstrand)==chr) & qsel
	pluscoord=start(plusstrand[select])
        qsel = (mcols(minusstrand)$mapq > mapq)
        if(any(is.na(qsel))) qsel = T
        select = (seqnames(minusstrand)==chr) & qsel
        minuscoord=start(minusstrand[select])
        lx=list(list(plus=pluscoord,minus=minuscoord))
        names(lx)=bamnames
        lx
})
names(allreads)=obschrnames

save(allreads,file=bamout)

}else{
#else multiple replicates to merge..
    bamnames = args[3:length(args)]

    bamlist=lapply(bamnames,function(bamname){
        print(bamname)
        plusflags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isMinusStrand=F,isDuplicate=F,isNotPassingQualityControls=F),what=c('mapq'))
        plusstrand = readBam(bamname,param=plusflags)
        #
        minusflags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isMinusStrand=T,isDuplicate=F,isNotPassingQualityControls=F),what=c('mapq'))
        minusstrand = readBam(bamname,param=minusflags)
        list(plusstrand,minusstrand)
    })

    obschrnames=unique(do.call(c,lapply(bamlist,function(i){
        levels(c(seqnames(i[[1]]),seqnames(i[[2]])))
    })))

    allreads=lapply(obschrnames,function(chr){
	print(chr)
        lx=lapply(bamlist,function(bam){
            qsel = (mcols(bam[[1]])$mapq > mapq)
            if(any(is.na(qsel))) qsel = T
            select = (seqnames(bam[[1]])==chr) & qsel
            pluscoord=start(bam[[1]][select])
            qsel = (mcols(bam[[2]])$mapq > mapq)
            if(any(is.na(qsel))) qsel = T
            select = (seqnames(bam[[2]])==chr) & qsel
            minuscoord=start(bam[[2]][select])
            list(plus=pluscoord,minus=minuscoord)
        })
        names(lx)=bamnames
        lx
    })
    names(allreads)=obschrnames

    save(allreads,file=bamout)

}
