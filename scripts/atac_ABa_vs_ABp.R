##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep 26 14:51:15 2018
##########################################################################
##########################################################################
peak.files = list.files(path = "../R6548_atac/Peaks/macs2",
                       pattern = "*.xls", full.names = TRUE)
# import tbx peaks
peak.files = c(peak.files, "data/tbx_90min_peaks_merged_macs2_p_5_filtered_N2_gene_assignment_TSS_WBcel235_analysis_v2.bed")

DIR.cwd = getwd();
outDir = paste0(DIR.cwd, "/results/DB_ABa_ABp")
if(!dir.exists(outDir)) dir.create(outDir);

###############################
# libraries and functions
###############################
library("ChIPseeker");
library("rtracklayer")
#library(UpSetR);library("ChIPpeakAnno")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2");
library("GenomicFeatures");
library("GenomicRanges")
#library("DiffBind");library("IRanges");
source('functions_chipSeq.R')

###############################
# make design matrix and import peak files
###############################
source('functions_chipSeq.R')
design.matrix = make.design.matrix.from.peaks.files(peak.files = peak.files)

cat("-- import peak files as GRanges objects \n")
peaks.list = c()
for(n in 1:nrow(design.matrix)){
  #cat(n, '\n')
  xx = readPeakFile(design.matrix$file.path[n], as = "GRanges");
  if(seqlevelsStyle(xx) != "UCSC") seqlevelsStyle(xx) = "UCSC";
  
  peaks.list = c(peaks.list, xx)
  #eval(parse(text = paste0("pp.", n, " = xx")))
}

###############################
# peak overlapping checking
###############################
cat("-- compare peak overlapping and make plots \n")
#kk = c(1:6)

pdf(paste0(outDir, "/Comparison_peaks_for_ABa_ABp_overlapping_atac_tbx.pdf"),
    width = 12, height = 8)
source("functions_chipSeq.R") #
#par(mfrow=c(1,2))
Comparison.overlapping.peaks(design.matrix, peaks.list, toCompare="condition", PLOT.p10 = TRUE)

dev.off()

###############################
# save the ABa, ABp unique peaks regions
###############################
library(rtracklayer)

ABa = peaks.list[[which(design.matrix$factor.condition == "Aba_90min")]]
ABp = peaks.list[[which(design.matrix$factor.condition == "Abp_90min")]]
p = ABa 
p <- p[mcols(p)[,"X.log10.pvalue."] > pval]
p <- p[!overlapsAny(p, ABp)]

export(p, con = paste0("results/motif_analysis/peaks_bed/ABa_unique_peaks.bed"))

p = ABp
p <- p[mcols(p)[,"X.log10.pvalue."] > pval]
p <- p[!overlapsAny(p, ABa)]
export(p, con = paste0("results/motif_analysis/peaks_bed/ABp_unique_peaks.bed"))
###############################
# run motif analysis for those two unique peak sets 
###############################




