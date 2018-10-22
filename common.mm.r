source('common.global.r')

#####
# Check before running

#minimum map quality
mapq = -1

# load genome file
suppressMessages(bis("BSgenome.Mmusculus.UCSC.mm10"))
genome = BSgenome.Mmusculus.UCSC.mm10

#####
# Run mode options

#if a file indicating completion is found; exit.
overwrite=T
#purity of calls
purity.cut = 0.7

#####
# Motif options

#motif score cutoff (log)
motifcut = 5
#max candidate sites (default 500k)
maxcand = 100000

#####
# Approx motif match options

#number of kmer samples to draw
nkmer = 5000000


######
# Optional: location of a whitelist directory. PWM matches will only be used if their window completely fits in the whitelist.
# whitelist should be in .bed format.

# WARNING: Do not use the whitelist to subset regions of interest in the genome such as DNase hypersensitive regions.
# This breaks the PIQ background model and will result in reduced performance
# This is only meant to be used to remove regions that have mapping, repeat or experiment problems.
# In all other cases, run PIQ with no whitelist, and subselect binding sites after calling.

#whitelist = '/cluster/thashim/PIQ/capture.mm10.bed'
whitelist = NULL

#example - encode hg19 blacklists
#system('wget -O /cluster/thashim/PIQ/enc.hg19.blacklist.bed.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz')
#blacklist = '/cluster/thashim/PIQ/enc.hg19.blacklist.bed.gz'
#disable blacklist by setting to NULL
blacklist = NULL

#should blacklist just prevent motif matches at the blacklist, or also drop [-wsize,wsize] regions around them.
flank.blacklist = T

###
# blacklist repeatmask, does not work for yeast.
remove.repeatmask = F

####
# plot options

# change the window size; useful to see a close-up for the footprint
plot.wind = wsize
plot.bothstrand = T
