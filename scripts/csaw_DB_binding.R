##################################################
## Project: Capture-seq analysis
## Script purpose: use csaw to do differential binding analysis for capture-seq 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jan 15 14:56:54 2018
##################################################
library(csaw)
library(edgeR)

DIR.bams = "../../R6329_R6532_R6533_chipseq_captured/alignments/BAMs_All"
DIR.OUT = '../results/normalization/'
resDir = DIR.OUT;

if(!dir.exists(DIR.OUT)) system(paste0('mkdir -p ', DIR.OUT))

bam.files = list.files(path = DIR.bams, pattern = "*.bam$", full.names = TRUE)
#bam.files <- c("es_1.bam", "es_2.bam", "tn_1.bam", "tn_2.bam")

## design matrxi (to define later)
design <- model.matrix(~factor(c('es', 'es', 'tn', 'tn')))
colnames(design) <- c("intercept", "cell.type")

########################################################
########################################################
# Section: Normalization
# calculate the scaling factors for normalization
# 1) we quantify read counts within windows using csaw 
# 2) select windows not releveant to PRC binding (and within the baits for capture-seq data)
#3 ) identify the median as scaling factors
########################################################
########################################################


#library(csaw)
frag.len <- 110
win.width <- 120
chr.selected = c(paste0("chr", c(1:18)))
param <- readParam(pe="none", dedup = TRUE, minq=30, restrict = chr.selected)

data <- windowCounts(bam.files, ext=frag.len, width=win.width, param=param, bin = TRUE)

## for normalization
binned <- windowCounts(bam.files, bin=TRUE, width=100000, param=param)
#binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)

save(data, binned, file = paste0(DIR.OUT, "readCounts_windows_csaw.Rdata"))

# check the windows and counts
head(assay(data))
head(rowRanges(data))

# counting read for G11 reporter 
#my.regions <- GRanges(c("chr15"),
#                      IRanges(c(78819673), c(78905009)))

#reg.counts <- regionCounts(bam.files, my.regions, ext=frag.len, param=param)
#head(assay(reg.counts))
load(file = paste0(DIR.OUT, "readCounts_windows_csaw.Rdata"))

## filtering uninteresting windows
Filter.With.Baits = TRUE
if(Filter.With.Baits)
{
  # import bait coordinates
  library(GenomicRanges)
  library(rtracklayer)
  baits = read.delim(file="../../Oliver/capture_seq/baits/baits_mm9G11D_withG11D.bed", sep = "\t", header = FALSE)
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

Filtering.binned = FALSE
if(Filtering.binned)
{
  abundances <- aveLogCPM(asDGEList(binned), prior.count = 1)
  summary(abundances)
  
  cutoff.average.counts = 5;
  binned.keep <- abundances > aveLogCPM(cutoff.average.counts, lib.size=mean(binned$totals))
  summary(binned.keep)
  #keep.simple <- abundances > -1
  binned.filtered <- binned[binned.keep,]
  
}

lib.sizes.full = binned$totals;
norms <- normOffsets(binned, type="scaling", lib.sizes=lib.sizes.full)

filtered.data$norm.factors = norms

par(mfrow=c(1, 1))
plot(norms, filtered.data$totals/10^6, log='')
abline(0, 1, lwd=2.0, col='red')

# to visulize the normalization
nb.for.test = 5
par(mfrow=c(1, nb.for.test), mar=c(5, 4, 2, 1.5))
adj.counts <- cpm(asDGEList(binned), log=TRUE)
normfacs <- filtered.data$norm.factors

for (i in 1:nb.for.test) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                xlab="A", ylab="M", main=paste("1 vs", i+1))
  all.dist <- diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")

}

##################################################
## Section: start to use DESeq2 to test differential binding
##################################################
coordiantes = as.data.frame(rowRanges(filtered.data)) 
ggs = apply(coordiantes[, c(1:3)], 1, function(x) gsub(" ", "",paste(as.character(unlist(x)), sep = "",  collapse = "_"), fixed = TRUE))
#ggs = paste0(coordiantes[, c(1:3)], collapse = "_") 
aa = data.frame(ggs, assay(filtered.data), stringsAsFactors = FALSE)
colnames(aa) = c('gene', basename(filtered.data$bam.files))

### test for differential binding
genotype <- factor(sapply(filtered.data$bam.files, function(x) unlist(strsplit(as.character(basename(x)), "_"))[1]))
IP <- factor(sapply(filtered.data$bam.files, function(x) unlist(strsplit(as.character(basename(x)), "_"))[2]))
design = data.frame(samples=basename(filtered.data$bam.files),genotype,IP)
rownames(design) = design$samples
norms = filtered.data$norm.factors;

## quality controls for all samples except inputs, because inputs seems to be not very good 
sels = which(design$IP != "Input")
read.count = aa[, c(sels+1)];
source("RNAseq_Quality_Controls.R")

pdfname = paste0(resDir, "Data_Qulity_Assessment_all.pdf")
pdf(pdfname, width = 12, height = 10)

Check.RNAseq.Quality(read.count=read.count, design.matrix = data.frame(sample=colnames(read.count), design[sels, c(3, 2)]),
                     keep.All = TRUE, norms = norms[sels])

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

#design <- model.matrix(~genotype+IP)
#rownames(design) <- colnames(y)
#design
#y = asDGEList(filtered.data)

#y = estimateDisp(y, design)
#summary(y$trended.dispersion)

#fit <- glmQLFit(y, design, robust=TRUE)
#summary(fit$var.post)
#par(mfrow=c(1,2))
#o <- order(y$AveLogCPM)
#plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
#     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
#     ylab=("Biological coefficient of variation"))
#plotQLDisp(fit)

## extrac the pairwise comparisons using contrast in edgeR
#results <- glmQLFTest(fit, contrast=c(0, 1))
#head(results$table)
#rowData(filtered.data) <- cbind(rowData(filtered.data), results$table)
# par(mfrow=c(2,2), mar=c(5,4,2,2))
# adj.counts <- cpm(y, log=TRUE)
# for (top in c(1000)) {
#   out <- plotMDS(adj.counts, main=top, col=c("blue", "blue", "red", "red"),
#                  labels=basename(filtered.data$bam.files), top=top)
# }




