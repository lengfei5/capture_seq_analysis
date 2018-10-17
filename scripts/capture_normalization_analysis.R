##################################################
## Project: Capture-seq analysis
## Script purpose: use csaw to do differential binding analysis for capture-seq 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jan 15 14:56:54 2018
##################################################
library(csaw)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)

DIR.bams = "../../R6329_R6532_R6533_chipseq_captured/alignments/BAMs_unique_rmdup"
resDir = '../results/DB_capture/'
NormDir = "../results/normalization/";
version.analysis = "_20181017"

#resDir = NormDir;

if(!dir.exists(resDir)) system(paste0('mkdir -p ', resDir))

bam.files = list.files(path = DIR.bams, pattern = "*.bam$", full.names = TRUE)
#bam.files <- c("es_1.bam", "es_2.bam", "tn_1.bam", "tn_2.bam")

## design matrxi (to define later)
#design <- model.matrix(~factor(c('es', 'es', 'tn', 'tn')))
#colnames(design) <- c("intercept", "cell.type")

########################################################
########################################################
# Section I : Normalization
# calculate the scaling factors for normalization
# 1) we quantify read counts within windows using csaw 
# 2) select windows not releveant to PRC binding (and within the baits for capture-seq data)
#3 ) identify the median as scaling factors
########################################################
########################################################
regions.sel = read.table(paste0(NormDir, "Complement_mergedPRCsubunits_peaks_fromJorgeChipseq_mm9G11_dual.bed"),
                         header = FALSE, sep = "\t")

sizes = regions.sel[, 3] - regions.sel[,2]
hist(log10(sizes))

baits = read.delim(file="../../../Oliver/capture_seq/baits/baits_mm9G11D_withG11D.bed", sep = "\t", header = FALSE)
baits = data.frame(baits, stringsAsFactors = FALSE)
colnames(baits) = c("chr", "start", "end")
baits$score=1;
size.baits = baits$end - baits$start
hist(log10(size.baits))

#frag.len <- 110
win.width <- 2000
chr.selected = c(paste0("chr", c(1:19)))
param <- readParam(pe="none", dedup = TRUE, minq=30, restrict = chr.selected)

binned <- windowCounts(bam.files, bin=TRUE, width=win.width, param=param)

## save the big file in case it is lost 
save(binned, baits, file = paste0(NormDir, "bam_unique_rmdup_readCounts_windows_csaw_forNormalization.Rdata"))

###############################
# select the PRC-unrelated regions and select just baits regions for captured-seq data 
###############################
load(file = paste0(NormDir, "bam_unique_rmdup_readCounts_windows_csaw_forNormalization.Rdata"))

Select.PRC.unrelated.Regions = TRUE
#Select.PRC.unrelated.Regions.inBaits.forCapturedData = TRUE

design = as.data.frame(colData(binned))
find.sampleID = function(x){
  #x = design$bam.files[1]
  x = unlist(strsplit(as.character(basename(x)), "_"))
  x = x[-c(length(x), (length(x)-1))]
  return(x[length(x)])
}

design$sampleID = sapply(design$bam.files, find.sampleID)
design$type = "captured"
design$type[match(c(71079:71090), design$sampleID)] = "chipseq"
design$condition = sapply(design$bam.files, function(x) unlist(strsplit(as.character(basename(x)), "_"))[1]) 
design$IP = sapply(design$bam.files, function(x) unlist(strsplit(as.character(basename(x)), "_"))[2]) 
design = data.frame(design, stringsAsFactors = FALSE)

if(Select.PRC.unrelated.Regions){
  sels = data.frame(regions.sel[, c(1:3)], stringsAsFactors = FALSE)
  colnames(sels) = c("chr", "start", "end")
  sels$score=1;
  df = makeGRangesFromDataFrame(sels)
  
  # filter windows out of baits
  suppressMessages(keep <- overlapsAny(rowRanges(binned), df, type = "within"))
  sum(keep)
  
  filtered.binned = binned[keep, ]
}
#binned.chipseq = filtered.binned[, which(design$type=="chipseq")]

save(filtered.binned, regions.sel, baits, design, file = paste0(NormDir, "normalization_factors_for_chipseq_captured_seq.Rdata"))


########################################################
########################################################
# Section: II 
# Differential Binding analysis and test the size facotrs of normalization for ChIP-seq data 
########################################################
########################################################
load(file = paste0(NormDir, "normalization_factors_for_chipseq_captured_seq.Rdata"))
DIR.peaks = "../../R6329_R6532_R6533_chipseq_captured/Peaks/macs2_broad"
peak.list = list.files(path = DIR.peaks, pattern = "*.xls", full.names = TRUE)
peak.list = peak.list[grep("710", peak.list)]

for(prot in c("Cbx7", "Ring1B")){
  
  # prot = "Cbx7"
  source("functions_chipSeq.R")
  
  peaks = merge.peaks.macs2(peak.list[grep(prot, peak.list)]);
  bams = design$bam.files[which(design$type=="chipseq" & design$IP==prot)]
  design.matrix = design[which(design$type=="chipseq" & design$IP==prot), ]
  
  source("functions_chipSeq.R")
  counts = quantify.signals.within.peaks(peaks, bam.list = bams)
  #colnames(counts) = basename(bams)
  
  Calculate.Scaling.factors.for.Chipseq = FALSE
  if(Calculate.Scaling.factors.for.Chipseq){
    
    binned.chipseq = filtered.binned[, which(design$type=="chipseq")]
    source("functions_analysis_captured.R")
    
    norms.chipseq = calcNormFactors.for.caputred.using.csaw(dd = binned.chipseq, method = "DESeq2", cutoff.average.counts = 300);
    source("functions_analysis_captured.R")
    norms.chipseq = calcNormFactors.for.caputred.using.csaw(dd = binned.chipseq, method = "DESeq2", cutoff.average.counts = 300);
    
    norms = norms.chipseq$size.factors[match(design.matrix$bam.files, norms.chipseq$bam.files)] 
    #design.matrix = data.frame(design.matrix)
  }
  
  kk = which(colnames(design.matrix) == "condition")
  
  pdfname = paste0(resDir, "Data_Qulity_Assessment_DB_analysis_ChIPseq_", prot, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  
  source("functions_chipSeq.R")
  res = DB.analysis(counts, design.matrix[, kk], size.factors = NULL, Threshold.read.counts = 50)
  
  dev.off()
  
}

###############################
#  make bigwig file with the size factors identified by previous step for ChIP-seq data
###############################
load(file = paste0(NormDir, "normalization_factors_for_chipseq_captured_seq.Rdata"))
library(GenomicAlignments)
library(rtracklayer)

bw.Dir = paste0(resDir, "normalizedBigWigs/")
if(!dir.exists(bw.Dir)) system(paste0('mkdir -p ', bw.Dir))

#norms = rbind(norms.chipseq, norms.captured)
Normalize.chipseq.make.BigWig = FALSE
if(Normalize.chipseq.data){
  norms = norms.chipseq
  norms = data.frame(norms)
  
  plot(norms$totals, norms$size.factors, log = 'xy')
  
  param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE, isSecondaryAlignment=FALSE),
                        mapqFilter = 30)
  
  for(n in 1:nrow(norms)){
    # n = 2
    cat(norms$bam.files[n], "---")
    ga = readGAlignments(norms$bam.files[n], param=param)
    xx = coverage(ga)
    xx = xx*norms$totals[n]/mean(norms$totals)/norms$size.factors[n]
    #xx = 
    bwname = basename(norms$bam.files[n])
    bwname = gsub(".bam", "_normalized.bw",bwname)
    cat(bwname, "\n")
    export.bw(xx, con = paste0(bw.Dir, bwname))
  }
  
}

########################################################
########################################################
# Section III: 
# identify signals in the baits for capture-seq data
########################################################
########################################################
load(file = paste0(NormDir, "normalization_factors_for_chipseq_captured_seq.Rdata"))

calculate.size.factors.by.select.PRC.unrelated.Regions.inBaits = TRUE

###############################
# quantify counts for baits with csaw
###############################
jj = which(design$type=="captured" & design$totals>10000)
bams = design$bam.files[jj]
design.matrix = design[jj, ]

#frag.len <- 110
win.width <- 200
chr.selected = unique(baits$chr)
chr.selected = chr.selected[order(chr.selected)]
param <- readParam(pe="none", dedup = TRUE, minq=30, restrict = chr.selected)

data <- windowCounts(bams, ext=110, width=win.width, param=param, bin = TRUE)

## for normalization
#binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
save(data, file = paste0(resDir, "readCounts_windows_200bp_binned_csaw_for_captured.Rdata"))

# check the windows and counts
head(assay(data))
head(rowRanges(data))

load(file = paste0(resDir, "readCounts_windows_200bp_binned_csaw_for_captured.Rdata"))

## filtering uninteresting windows
Filter.With.Baits = TRUE
if(Filter.With.Baits){
  # import bait coordinates
  library(GenomicRanges)
  library(rtracklayer)
  baits = read.delim(file="../../../Oliver/capture_seq/baits/baits_mm9G11D_withG11D.bed", sep = "\t", header = FALSE)
  baits = data.frame(baits, stringsAsFactors = FALSE)
  colnames(baits) = c("chr", "start", "end")
  baits$score=1;
  df = makeGRangesFromDataFrame(baits)
  
  # filter windows out of baits
  suppressMessages(keep <- overlapsAny(rowRanges(data), df))
  sum(keep)
  
  filtered.data = data[keep, ]
  
}else{ # otherwise filter with the cutoff of average counts  
  abundances <- aveLogCPM(asDGEList(data), prior.count = 1)
  summary(abundances)
  
  cutoff.average.counts = 5;
  keep <- abundances > aveLogCPM(cutoff.average.counts, lib.size=mean(data$totals))
  summary(keep)
  #keep.simple <- abundances > -1
  filtered.data <- data[keep,]
  #summary(keep.simpl)
}

if(calculate.size.factors.by.select.PRC.unrelated.Regions.inBaits){
  
  binned.captured = filtered.binned[, which(design$type=="captured")]
  
  df = makeGRangesFromDataFrame(baits)
  
  # filter windows out for captured-seq data
  suppressMessages(keep <- overlapsAny(rowRanges(binned.captured), df, type = "within"))
  sum(keep)
  
  binned.captured = binned.captured[keep, ]
  
  #kk = which(binned.captured$totals>100000)
  #binned.captured = binned.captured[, kk]
  
  source("functions_analysis_captured.R")
  norms.captured = calcNormFactors.for.caputred.using.csaw(dd = binned.captured, method = "DESeq2", cutoff.average.counts = 1000);
  norms.captured$library.sizes = colSums(assay(binned.captured))
  
  norms.captured = norms.captured[match(design.matrix$bam.files, norms.captured$bam.files), ]
  
  par(mfrow=c(1,3))
  plot(norms.captured$totals, norms.captured$size.factors, log='xy')
  plot(norms.captured$totals, norms.captured$library.sizes, log='xy')
  plot(norms.captured$library.sizes, norms.captured$size.factors, log='xy')
  
}

###############################
# Section: start to use DESeq2 to test differential binding
###############################
coordiantes = as.data.frame(rowRanges(filtered.data)) 
ggs = apply(coordiantes[, c(1:3)], 1, function(x) gsub(" ", "",paste(as.character(unlist(x)), sep = "",  collapse = "_"), fixed = TRUE))
#ggs = paste0(coordiantes[, c(1:3)], collapse = "_") 
aa = data.frame(ggs, assay(filtered.data), stringsAsFactors = FALSE)
colnames(aa) = c('gene', basename(filtered.data$bam.files))

### test for differential binding
#genotype <- factor(sapply(filtered.data$bam.files, function(x) unlist(strsplit(as.character(basename(x)), "_"))[1]))
#IP <- factor(sapply(filtered.data$bam.files, function(x) unlist(strsplit(as.character(basename(x)), "_"))[2]))
#design = data.frame(samples=basename(filtered.data$bam.files),genotype,IP)
rownames(design.matrix) = basename(design.matrix$bam.files)
design.matrix$size.factors = norms.captured$size.factors
design.matrix$library.sizes = norms.captured$library.sizes/median(norms.captured$library.sizes)

## quality controls for all samples except inputs, because inputs seems to be not very good 
#sels = which(design$IP != "Input")

pdfname = paste0(resDir, "Data_Qulity_Assessment_all_captured_signals_inBaits_librarySizes_nonPRC_regions_each_IP.pdf")
pdf(pdfname, width = 12, height = 10)

source("RNAseq_Quality_Controls.R")

for(prot in unique(design.matrix$IP)){
  #sels = c(1:nrow(design.matrix))
  cat("IP --- ", prot, "\n")
  sels = which(design.matrix$IP == prot)
  read.count = aa[, c(sels+1)];
  
  Check.RNAseq.Quality(read.count=read.count, 
                       design.matrix = data.frame(sample=colnames(read.count), design.matrix[sels, c(7)]), 
                       norms = design.matrix$library.sizes[sels],
                       keep.All = TRUE)
}

dev.off()


require('DESeq2')

design.matrix = design
countData = as.matrix(aa[, -1])
rownames(countData) = aa$gene;

conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
#dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
sizeFactors(dds) = norms;

fpm = fpm(dds, robust = TRUE)

## save windonw informations and associate them wiht baits
coordinates = data.frame(coordiantes, stringsAsFactors = FALSE);
coordinates$baits = NA;
for(n in 1:nrow(baits))
{
  #n = 6;
  kk1 = which(as.character(coordiantes[,1])==as.character(baits$chr[n]))
  kk2 = which(as.numeric(coordiantes[,2])>= as.numeric(baits$start[n])) 
  kk3 = which(as.numeric(coordiantes[,3])<=as.numeric(baits$end[n]))
  kk = intersect(kk1, kk2)
  kk = intersect(kk, kk3)
  coordinates$baits[kk] = paste0(c(as.character(baits[n, 1]), baits[n, c(2,3)]), collapse = "_")   
}


## Summary plots of comparisions across 3 conditions (G11D, G11D.TetR.Cbx7 and G11.TetR_Rybop) for the 5 IPs
all.baits = unique(coordinates$baits)
all.baits = all.baits[which(!is.na(all.baits))]

cols = c("darkgray", "blue", "darkorange")
ccs = c("G11D", "G11D.TetR.Cbx7", "G11D.TetR.Rybp");
lwd = 2.0

for(sel.bait in all.baits)
{
  pdfname = paste0(resDir, "Comparisons_across_three_conditions_for_each_IP_for_bait_",sel.bait, ".pdf")
  pdf(pdfname, width = 20, height = 14)
  for (ip in unique(design$IP))
  {
    #sel.bait = "chr15_78819673_78905009";
    #ip = "H3K27me3";
    
    cat('bait -- ', sel.bait, '-- IP --', ip, "\n")
     
    kk = which(coordinates$baits==sel.bait)
    jj = which(design$IP==ip)
     
     par(mfrow=c(2,2))
     for(bandwith4kernel in c(120, 500, 2000, 5000))
     {
       plot(0.1, 1, xlim= range(coordinates$start[kk], na.rm = TRUE), ylim=range(fpm[kk,jj], na.rm = TRUE), 
            main=paste0(ip, "--", sel.bait, "--bandwith--", bandwith4kernel, "bp"), xlab="position", ylab='normalized read counts', log='', type='n')
       
       for(j in jj)
       {
         if(design$genotype[j]==ccs[1]) col= cols[1]
         if(design$genotype[j]==ccs[2]) col= cols[2]
         if(design$genotype[j]==ccs[3]) col= cols[3]
         
         x = coordinates$start[kk]; y = fpm[kk, j];
         
         if(bandwith4kernel==120){
           points(x, y, col=col, type='l', lwd=lwd)
         }else{
           lines(ksmooth(x, y, "normal", bandwidth = bandwith4kernel), col = col, lwd=lwd,  type='l') 
         }
       }
       legend("topright", ccs,  col = cols, lwd=2, bty = 'n', cex = 1.2)
     }
     
  }
  dev.off()
}
