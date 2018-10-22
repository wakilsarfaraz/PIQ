#####
# utils

old.repos <- getOption("repos")
on.exit(options(repos = old.repos))
new.repos <- old.repos
new.repos["CRAN"] <- "http://cran.stat.ucla.edu"
options(repos = new.repos)

source("http://www.bioconductor.org/biocLite.R")
#if(!suppressMessages(require("BiocInstaller",quietly=T))){
#    install.packages("BiocInstaller", repos="http://www.bioconductor.org/packages/2.13/bioc")
#}

ris <- function(x){if(!require(x,character.only=T,quietly=T,warn.conflicts=F)){
    install.packages(x)
    require(x,character.only=T,quietly=T,warn.conflicts=F)
}}
bis <- function(x){if(!require(x,character.only=T,quietly=T,warn.conflicts=F)){
    biocLite(x)
    require(x,character.only=T,quietly=T,warn.conflicts=F)
}}

wipetemp <- function(){
    x <- readline("this will wipe existing temp files (y/n)")
    if(x=="y"){
        print("wiping temp")
        unlink("tmp/*")
    }else{
        print("keeping temp")
    }
}

set.seed(1)

####
# dependencies

suppressMessages(require('BiocInstaller',quietly=T))
suppressMessages(bis("Rsamtools"))
suppressMessages(bis("GenomicAlignments"))
suppressMessages(bis('Biostrings'))
suppressMessages(bis('seqLogo'))

suppressMessages(ris('RSofia'))
suppressMessages(ris('statmod'))
suppressMessages(ris('Rcpp'))
suppressMessages(ris('inline'))
suppressMessages(require(Matrix,quietly=T))

####
# params


#window size
wsize = 300
#fast mode
fast.mode = T

#seqbias mode
avoid.seqbias = F
