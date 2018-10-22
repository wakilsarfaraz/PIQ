#####
# Below is an example of how to do DNase-seq calls using a sge based cluster
#####

#Location of the common.r file that does package load and defines parameters.
commonfile="/cluster/thashim/basepiq/common.r"

#Directory in which the PWM files were outputted using pwmmatch
pwmdir="/cluster/thashim/PIQ/140320.hg19.exact.pwm/"

#Temporary read/write dir. Should be fast and local to avoid bottlenecking
tmpdir="/scratch/tmp/"

#Output directory to which all call files are written
outdir="/cluster/thashim/PIQ/140320.hg19.exact.call/"

#Processed bamfile created by bam2rdata script.
bamfile="/cluster/thashim/PIQ/k562_comb.RData"

#Number of motifs to process
for pwmid in {1..1316}
do
	echo "/cluster/thashim/basepiq/pertf.r "$commonfile $pwmdir $tmpdir $outdir $bamfile $pwmid | /usr/bin/qsub -wd /cluster/thashim/basepiq/ -e $outdir/err.txt -o /dev/null
done