#!/bin/bash

cd ..
outdir=testout
mkdir $outdir
wget http://piq.csail.mit.edu/data/bams/k562.bam -O $outdir/k562.bam

mkdir $outdir/pwms
Rscript pwmmatch.exact.r common.r pwms/jasparfix.txt 139 $outdir/pwms

Rscript bam2rdata.r common.r $outdir/k562.RData $outdir/k562.bam

mkdir $outdir/tmp
mkdir $outdir/out
Rscript pertf.r common.r $outdir/pwms/ $outdir/tmp/ $outdir/out/ $outdir/k562.RData 139
