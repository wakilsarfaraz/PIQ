IMPORTANT - READ ME
======================================

Frequently asked questions
-----

If you have a question [read me first](https://bitbucket.org/thashim/piq-single/wiki/FAQ)

common.r
-----
MODIFY genome to fit your data, match bis("BSgenome.*") to your genome

As of bioconductor 2.14 some of the package names have changed. Description below:

Organism | package name | description
----------|---------------|--------------
Human | BSgenome.Hsapiens.UCSC.hg19.masked | Human genome + repeatmask
Human | BSgenome.Hsapiens.UCSC.hg19| Human genome (set mapq=0 or blacklist)
Mouse | BSgenome.Mmusculus.UCSC.mm10| Mouse genome (set mapq=0 or blacklist)
Yeast | BSgenome.Scerevisiae.UCSC.sacCer2 | Yeast (set mapq=0)

In general we suggest repeatmasking the genome when using unique maps (mapq>=1). hg19.masked does this by default for mm10 and sacCer2, the masked genomes do not have repeatmask, and we'd suggest blacklisting them via a blacklist or use non-unique maps. 




input files
----
use jasparfix.txt in the pwms folder (all JASPAR including PBM hits) for a comprehensive list. If for whatever reason you only want JASPAR CORE hits, use jaspar.txt.


How to run
=======================================

Running a single motif
-----------------------

0. Make sure you can execute common.r without errors (this make sure your R is set up correctly).

1. Generate the PWM hits across genome (does *not* depend on choice of BAM)

    `Rscript pwmmatch.exact.r /cluster/thashim/basepiq/common.r /cluster/thashim/basepiq/pwms/jasparfix.txt 139 /cluster/thashim/PIQ/motif.matches/`

    This uses the genome and PWM cutoffs in common.r with the 139th motif in jaspar.txt (CTCF) and writes the matches as a binary R file called 139.RData in tmppiq.

2. Convert BAM to internal binary format (does *not* depend on choice of motif).

    `Rscript bam2rdata.r /cluster/thashim/basepiq/common.r /cluster/thashim/PIQ/d0.RData /cluster/cwo/dnase_seq/bams/D0_50-100_130801.bwa.mapq20.mm10.bam`

    This takes the D0_50-100_130801.bwa.mapq20.mm10.bam file and retains only the read end locations, storing it into d0.RData.

3. Combine BAM + PWM to make calls, depends on choice of BAM and PWM. If you have multiple simultaneous runs, make sure each run gets a unique tmp folder.

    `Rscript pertf.r /cluster/thashim/basepiq/common.r /cluster/thashim/PIQ/motif.matches/ /scratch/tmp/ /cluster/thashim/130130.mm10.d0/ /cluster/thashim/PIQ/d0.RData 139`

    This takes settings in common.r to call the motif match 139.RData (CTCF, from above) using data d0.RData (from bam2rdata) and writing the output into 130130.mm10.d0. The tmp folder stores some large temporary matrices and can be wiped after the run.

Running multiple motifs
----------------------

Bash scripts for running multiple motifs are provided as part of the utils/sgepwm.sh and utils/sgecalls.sh

Special use cases
----------------------

#### Multiple replicates

If you have replicates, please merge them using samtools merge.

#### Multiple experiments

If you have multiple experiments in the same experimental condition but different assays, bam2rdata.r will take multiple bamfiles as part of its argument and use them simultaneously.

#### Control experiments

If you have control experiments such as genomic DNA or naked DNase-I digestion and want to know if your results are significant with respect to control. Do the following:

    `Rscript pertf.bg.r /cluster/thashim/basepiq/common.r /cluster/thashim/PIQ/motif.matches/ /scratch/tmp/ /cluster/thashim/130130.mm10.d0/ /cluster/thashim/PIQ/d0.RData /cluter/thashim/PIQ/control.RData 139`

    This will use the profile / footprint learned from the training data (d0.RData) and output scores over the control data.


#### Train and test separately

The control experiments case can be used to learn footprints from a separate dataset than used for calls. This may be useful if

* You have multiple low-coverage DNase experiments in different conditions. Combine data to train but make calls separately.

* You have differential DNase-seq and want to make sure the same footprints are used across experiments to reduce false positives.


Parameters
===================================

Check before every run
----
* **Genome type**

Quality control / blacklisting
----
* **mapq**: reject reads in the bam with quality scores below mapq. Default setting enforces unique maps for bwa.
* **remove.repeatmask**: by default any BSgenome mask is enabled, for human this removes any repeatmasked regions. For mouse this disables all assembly gaps but not repeatmask
* **blacklist / whitelist**: if you have reason to believe there are regions that are not representative, use blacklist and whitelist in BED file format. **do not use whitelist to select hypersensitive regions**

Changing motif match candidates
----

Note: Increasing motif match candidate count may adversely affect performance due to inability to distinguish non-target TF binding with similar DNase profiles.

* **motifcut**: log-odds score motif cutoff vs uniform unigram background. default is generally good
* **maxcand**: never allow a PWM to match more than this many sites. Helps control degenarate PWM matchs.





About the output
======================================

For each PWM, PIQ outputs 3 files: two .csv files containing call scores and a .pdf containing summary and goodness of fit info.

*-calls.all.csv
--------------------------
This file contains ALL instances of motif matches from pwmmatch.r and its associated scores. This file is useful if you want to re-do the cutoffs yourself or if you want to look at instances where the motif is matched but the factor is unbound.



* **column 1**: row id, useful if going from calls.csv and mapping to calls.all.csv
* **column 2**: chromosome id, taken from the define BSgenome package
* **column 3**: coordinate uses same coordinates as BAM file
* **column 4**: pwm score, interpreted as log-probability vs uniform ATCG background.
* **column 5**: shape, interpreted as the logistic odds between PWM match and background, or the log-odds of this example being drawn from a PWM match rather than a background item. Do *not* interpret this as a probability by passing through a logit transform.
* **column 6**: score, weighted score of both PWM and shape (distance from linear separator in page 3, panel 4 of diag.pdf) same units as column 4 and 5 of log-odds
* **column 7**: purity, or an esimate of PPV: If you cutoff at purity value of X, then X% of those will be true binding sites, assuming a 50-50 prior on a candidate site being a binding site.


*-calls.csv
------
This file contains instances of motif match whose scores had purity above the purity.cuotff value in common.r. Purity is a proxy for PPV. A purity of 0.7 means 70% of the entries in calls.csv will be bound (estimated using background binding sites).

*-diag.pdf
-----
This file contains diagnostic pdf output of the motif and its binding sites.

### page 1 (input motif used to match)

### page 2 (global statistics)
* **Panel 1** (top left) is the PPV estimate vs rank of binding site as sorted by 'score' column in output, this helps estimate the number of sites bound.
* **Panel 2** (top right) is the observed best tradeoff between PWM score and DNase profile, higher values mean more weight is put on the DNase-shape compared to sequence.
* **Panel 3** (bottom left) is the correlation corrected footprint and linear decision bound for whether something is bound.
* **Panel 4** (bottom right) is the observed counts in each strand versus a set of random background sites (green)

### page 3 (site to site statistics)
* **Panel 1** shows distribution of scores vs random background taken by sampling; difference in two distributions are indicative of detectable binding. Vertical line indicates purity cutoff

* **Panel 2** is dnase reads (x axis) vs dnase score from piq (yaxis), showing whether the binding site calls are based on shape or overall counts. Lower coverage experiments will rely less on shape.
* **Panel 3** is the histogram representation of panel 2.
* **Panel 4** plots pwm score (xaxis) vs dnase score (yaxis). Black points are 'real' pwm match sites and their score. Red points are random background sites with dnase score taken from their reads, and assigned a random PWM score from the real pwm match set via permutation. Bold black points are those points called as 'bound'.