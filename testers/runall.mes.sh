#!/bin/bash

####
# Run specific stuff

jobname="mm10mes"
bampath="http://piq.csail.mit.edu/data/bams/mES.bam"

#####
# Below is an example of how to do DNase-seq calls using a sge based cluster
#####

tmpdir="/scratch/tmp/"
basedir="$(readlink -f $(pwd -P)/../)"
baseoutdir="/cluster/thashim/PIQ/"

pushd $basedir
jobid="$(date +"%y%m%d")-$(git rev-parse --short HEAD)"
idname="$jobid-$jobname"
popd

#Location of the common.r file that does package load and defines parameters.
commonfile="$basedir/common.mm.r"
outdir="$baseoutdir/$idname.calls/"
bamin="$baseoutdir/bams/$jobname.bam"
bamfile="$baseoutdir/rdata/$jobname.RData"

#Jaspar file directory path
jaspardir="$basedir/pwms/jasparfix.txt"
#Directory in which the PWM files should be dumped
pwmdir="$baseoutdir/$idname.pwms/"

wget $bampath -O $bamin 1> NUL 2> NUL

#bam proc
pushd $basedir
./bam2rdata.r $commonfile $bamfile "$bamin"
popd

mkdir $pwmdir
#Number of motifs to process
for pwmid in {1..1316}
do
	  echo "$basedir/pwmmatch.exact.r "$commonfile $jaspardir $pwmid $pwmdir | /usr/bin/qsub -wd $basedir/ -e $pwmdir/errfile.txt -o /dev/null -N piqpwm-$pwmid-$idname
done
cp $jaspardir $pwmdir
cp $commonfile $pwmdir
git log > $pwmdir/gitlog.txt
git diff > $pwmdir/gitdiff.txt

mkdir $outdir
##
# now call
for pwmid in {1..1316}
do
	echo "$basedir/pertf.r "$commonfile $pwmdir $tmpdir $outdir $bamfile $pwmid | /usr/bin/qsub -wd $basedir/ -e $outdir/err.txt -o /dev/null -N piqcal-$pwmid-$idname -hold_jid piqpwm-$pwmid-$idname
done
cp $jaspardir $outdir
cp $commonfile $outdir
git log > $outdir/gitlog.txt
git diff > $outdir/gitdiff.txt

####
# dump to web
####

until [ $(qstat | wc -l) == 0 ]
do
  echo -n "."
  sleep 1
done

dumpdir="/cluster/www/piq/data/$idname.calls"
mkdir $dumpdir
tar -zcf $dumpdir/$idname.calls.tar.gz $outdir
