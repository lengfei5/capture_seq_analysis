##################################################
## Project:capture-seq analysis
## Script purpose: make gene annotation for the customized geneome mm9_G11D
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jan  9 11:45:33 2018
##################################################

## inputs are bed files and a table for gene and transcripts mapping (Refseq mm9 here)
DIR.cwd="../../../annotations/mm9/"
bed = read.table(paste0(DIR.cwd, "Refseq_mm9_downloaed_UCSC.bed"), sep = "\t", header=FALSE, as.is = c(2, 3, 5, 7))
mapping = read.table(paste0(DIR.cwd, "Refseq_geneName_transcripts_mm9_downloaed_UCSC.txt"), sep = "\t", header = FALSE)
