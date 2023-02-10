# USDA-fermenter_data
## R script to run dada2 algorithm on trimmed paired-end reads
## NORMAL error learning - for non-binned quality scores
## ---- SET VARIABLES ----
## modify these variables:
# path to filtered and cleaned reads
CLEANEDPATH = "/gpfs/scratch/rkl5473/filtered"
list.files(CLEANEDPATH)
# path to output tables
#OUTPATH = "/gpfs/group/evk5387/default/Novogene/metformin"
OUTPATH = "/gpfs/scratch/rkl5473/"
# database to silva training set
# downloaded from: https://benjjneb.github.io/dada2/training.html
DB = "/storage/home/rkl5473/silva_nr_v132_train_set/silva_nr_v132_train_set.fa"
# paired end patterns 
FILTEREDF = "F_filt.fastq.gz"
FILTEREDR = "R_filt.fastq.gz"

# ---- install and read data ----
### packages must be previously installed
require(dada2)
require(tidyr)
require(phyloseq)
require(magrittr)
#list.files(CLEANEDPATH)
## test that pathway works
#if(!list.files(CLEANEDPATH)) {
#  cat("Can't read file pathway or files are not present")
#}
## ---- core dada algorithm ----
# get forward and reverse reads
forward <- sort(list.files(CLEANEDPATH, pattern = FILTEREDF, full.names = TRUE))
reverse <- sort(list.files(CLEANEDPATH, pattern = FILTEREDR, full.names = TRUE))
# check to make sure that the lengths of both files are the same and that they match
fwdNames <- sapply(strsplit(basename(forward), FILTEREDF), `[`, 1)
revNames <- sapply(strsplit(basename(reverse), FILTEREDR), `[`, 1)
# error catch
if(length(fwdNames) != length(revNames)) {
  stop("The number of forward and reverse files do not match.")
} else {
  
  if(any(!fwdNames%in% revNames)) {
    
    stop("Forward and reverse reads are out of order.")
  }
}
# learn errors for forward and reverse
errF <- learnErrors(
  forward,
  multithread = TRUE,
  verbose = TRUE
)
errR <- learnErrors(
  reverse,
  multithread = TRUE,
  verbose = TRUE
)
# save progress
save.image(file = paste0(OUTPATH, "/error-learning.RData"))
# plot forward errors
pdf(paste0(OUTPATH, "/forward-errorplot.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()
# plot reverse errors
pdf(paste0(OUTPATH, "/reverse-errorplot.pdf"))
plotErrors(errR, nominalQ = TRUE)
dev.off()
