##########################################################################
##########################################################################
# Project: Daniel's captured-seq data
# Script purpose: to identify PRC-unrelated geneomic regions based on which the normalization will be done 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep  3 14:24:13 2018
##########################################################################
##########################################################################
version.DATA = '20180903'
version.analysis = paste0(version.DATA, "_20180903")

resDir = paste0("results/normalization")
if(!dir.exists(resDir)){dir.create(resDir)}

peakDir = "../../Jorge/Analysis_ChIP_seq/Results/peaks_GCc_merged_with_geneAssignment"

###############################
# important peak coordinates from Jorge's ChIP-seq data and identify PRC-related regions by merging them 
###############################
library("GenomicFeatures");
library("GenomicRanges")
library("ChIPseeker");
xlist <-list.files(path=peakDir, pattern = "AN312*", full.names = TRUE)

peaks = c()
for(n in 1:length(xlist)){
  cat(n, '\n')
  xx = readPeakFile(xlist[n], as = "GRanges");
  if(n==1){
    peaks = xx
  }else{
    peaks= GenomicRanges::union(peaks, xx, ignore.strand=TRUE)
  }
  #peaks.list = c(peaks.list, xx)
  #eval(parse(text = paste0("pp.", n, " = xx")))
}

merged = reduce(peaks, drop.empty.ranges=TRUE, min.gapwidth=2000L, with.revmap=TRUE, with.inframe.attrib=FALSE)
pp = merged
#merged = as.data.frame(merged)
df <- data.frame(seqnames=seqnames(pp), starts=start(pp)-1, ends=end(pp), names=c(rep(".", length(pp))),
                 scores=c(rep(".", length(pp))), strands=strand(pp))
#kk = which(!is.na(match(df[,1], chroms)));
#df = df[kk, ]
bedname = paste0(resDir, '/mergedPRCsubunits_peaks_fromJorgeChipseq_mm10.bed')

write.table(df, file=bedname, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

###############################
# 1) next step: find the complement with bedtools
# 2) convert mm10 coordinate to mm9
###############################

###############################
# 3) convert mm9 to mm9_G11dual genome 
# because the G11 dual reporter ocuppied the chr15:78,845,690- 78,845,691 and took 9923 bp (length)
# we just need to change the cooridates if they are behind chr15:78,845,690- 78,845,691 by add 9923bp
###############################
mm9bed = paste0(resDir, "/Complement_mergedPRCsubunits_peaks_fromJorgeChipseq_liftovered_mm9.bed")
aa = read.table(mm9bed, sep = "\t")
jj = which(aa$V1=="chr15" & aa$V2>=78845691)
#xx = aa[jj, ]
aa$V2[jj] = aa$V2[jj] + 9923 
aa$V3[jj] = aa$V3[jj] + 9923

#bb = rbind(aa[-jj, ], xx)
#bb = data.frame(bb)

#bb = bb[order(bb$V1, bb$V2), ]
mm9G11dual = paste0(resDir, "/Complement_mergedPRCsubunits_peaks_fromJorgeChipseq_mm9G11_dual.bed")
write.table(aa, file=mm9G11dual, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

########################################################
########################################################
# Section: calculate scaling factors using PRC-unrelated regions for chip-seq and capture-seq data  
########################################################
########################################################





