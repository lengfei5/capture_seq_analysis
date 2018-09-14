#setwd("~/clustertmp/Jorge_Arturo.Zepeda_Martinez/")
library(ChIPseeker)
library(rtracklayer)
#library(UpSetR)
library("ChIPpeakAnno")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2")
library("GenomicFeatures")

venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  
  Venn(SetNames=SetNames, Weight=Weight)
}

venn_cnt2barplot <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts")

  cnt <- venn_cnt[,(n+1):length(colnames(venn_cnt))]
  rownames(cnt) <- cbind(apply(venn_cnt[,1:(n-1)], 1, paste, collapse="|"))
  cnt <- cnt[rownames(cnt) != "0|0|0", ]
  cnt <- melt(cnt)
  colnames(cnt) <- c("intersections", "condition", "peaks")
  legend.title <- paste(colnames(venn_cnt[,1:(n-1)]), collapse="|")

  p <- ggplot(cnt, aes(x=condition, y=peaks, fill=intersections)) +  geom_bar(stat="identity", colour="grey50", alpha = 0.9) + scale_fill_brewer(palette="Set3") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank()) + labs(fill=legend.title)

  plot(p)  
}


filterPeaks <- function(files, pval = 10, qval=0, remove.blacklist=FALSE, barplot = FALSE, fileType="macs2", input1B=TRUE, addChr=FALSE)
{
  # files = macs2.files;remove.blacklist=FALSE; barplot = FALSE; fileType="macs2"; input1B=TRUE; addChr=FALSE; pval = 100; qval=0;
  bar <- NULL
  peaks <- NULL
  for (f in files)
  {
    p <- readPeakFile(f, as = "GRanges")
    if (addChr)
    {
      seqlevels(p) <- paste0("chr", seqlevels(p))
    }
    if (input1B)
    {
      start(p) <- start(p) - 1
    }
    all <- length(p)
    if (fileType == "macs14")
    {
      p <- p[mcols(p)[,"X.10.log10.pvalue."] > pval]
    } else if (fileType == "macs2") {
      if (qval == 0)
      {
        p <- p[mcols(p)[,"X.log10.pvalue."] > pval]
      } else {
        p <- p[mcols(p)[,"X.log10.qvalue."] > qval]
      }
    }
    p100 <- length(p)
    if(remove.blacklist)
    {
      p <- p[!overlapsAny(p, blacklist)]
      p.bl <- p100 - length(p)
    }
    all <- all - p100
    p100 <- length(p)
    
    if(remove.blacklist){
      bar <- cbind(bar, c(all, p100, p.bl))
    }else{
      bar <- cbind(bar, c(all, p100))
    }
    
    peaks <- c(peaks, p)
  }
  
  names(peaks) <- names(files)
  
  if (barplot)
  {
    if (qval !=0)
    {
      rownames(bar) <- c("all", paste0("q", pval), "blacklisted")
    } else {
      rownames(bar) <- c("all", paste0("p", pval), "blacklisted")
    }
    colnames(bar) <- names(peaks)
    p <- ggplot(melt(bar), aes(x=X2, y=value, fill=X1)) + geom_bar(stat="identity", colour="black") + scale_fill_brewer() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.title=element_blank())
    plot(p)
  }
  
  return(peaks)  
  
}


plotV <- function(peaks, main = "")
{
    if (length(peaks) > 5)
    {
        peaks <- peaks[1:5] #up to five!!!
    }
    ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=names(peaks), maxgap=0, minoverlap =1, main=main, connectedPeaks="keepAll")
 
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    #plot(v, doWeights=FALSE)
    #plot(v, doEuler = TRUE)
    venn_cnt2barplot(ol.peaks$vennCounts)
    return(ol.peaks)
}


doAllPeak <- function(set)
{
  if (set == "macs14")
  {
    allPeaks <- filterPeaks(macs14.files, barplot=TRUE)
    allPeaks.50 <- filterPeaks(macs14.files, pval=50, blacklist=blacklist)
  } else if (set == "macs2")
  {
    allPeaks <- filterPeaks(macs2.files, pval = 10, blacklist=blacklist, barplot=TRUE, fileType="macs2")
    allPeaks.50 <- filterPeaks(macs2.files, pval = 5, blacklist=blacklist, barplot=TRUE, fileType="macs2")
  } else if (set == "macs2qval")
  {
    allPeaks <- filterPeaks(macs2.files, qval = 10, blacklist=blacklist, barplot=TRUE, fileType="macs2")
    allPeaks.50 <- filterPeaks(macs2.files, qval = 5, blacklist=blacklist, barplot=TRUE, fileType="macs2")
  } else if (set == "rep")
  {
    allPeaks <- filterPeaks(rep.files, blacklist=blacklist, barplot=TRUE)
    allPeaks.50 <- filterPeaks(rep.files, pval=50, blacklist=blacklist)
  }
  
  promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000) #ChIPseeker method
  tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  
  tagMatrixList <- tagMatrixList[!sapply(lapply(tagMatrixList, nrow), is.null)]
  tagMatrixList[sapply(tagMatrixList, nrow) == 0] <- NULL
  
  #promoter <- promoters(genes(txdb), upstream=5000, downstream=5000) #GRanges
  #tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  #plotAvgProf(tagMatrixList, xlim=c(-5000, 5000))
  #tagHeatmap(tagMatrixList, xlim=c(-5000, 5000), color=NULL)
  
  peakAnnoList <- lapply(allPeaks, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 200), verbose=FALSE)
  
  return(list(allPeaks, allPeaks.50, peakAnnoList, tagMatrixList))
  
}


doAllPlot <- function(set, allPeaksList)
{
    allPeaks <- allPeaksList[[1]]
    allPeaks.50 <- allPeaksList[[2]]
    peakAnnoList <- allPeaksList[[3]]
    tagMatrixList <- allPeaksList[[4]]
    
    if (set == "rep")
      {
        plotV(allPeaks, main = "Ints (< 1e-100)") #different order due to error otherwise ...
        plotV(allPeaks.50, main = "Ints (< 1e-50)")
      } else {
    
        plotV(allPeaks, main = "Ints (< 1e-100)") #different order due to error otherwise ...
        plotV(allPeaks.50, main = "Ints (< 1e-50)")
      }

    #covplot(allPeaks[[9]], weightCol="V5", chrs = c(paste0("chr", c(1:19, "X", "Y"))))


    #peakAnnoList <- list("Normal"=annotatePeak(allPeaks[[1]],  TxDb=txdb,tssRegion=c(-2000, 200), verbose=FALSE), 
    #"TSS"=annotatePeak(allPeaks[[1]],  TxDb=txdb,tssRegion=c(0, 0), verbose=FALSE), "5kb"=annotatePeak(allPeaks[[1]],  
    #TxDb=txdb,tssRegion=c(-5000,5000), verbose=FALSE), "3kb"=annotatePeak(allPeaks[[1]],  TxDb=txdb, verbose=FALSE))
    
    print(plotAnnoBar(peakAnnoList))
    print(plotDistToTSS(peakAnnoList))
    
    for (n in names(peakAnnoList))
      {
        par(mfrow=c(1,1))
        vennpie(peakAnnoList[[n]])
        text(x=0, y=-1, n)
        par(mfrow=c(1,1))
        plotAnnoPie(peakAnnoList[[n]])
#        par(mfrow=c(1,1))
#        upsetplot(peakAnnoList[[n]]) #only in TB-3.2.1-dev #  vennpie=TRUE,
      }
}

#for (set in c("macs14", "macs2", "macs2qval"))
#  {
#    pdf(paste0("peaksAnalysis.", set, ".pdf"), width = 10, height=10)
#    doAllPeak(set=set)
#    dev.off()
#  }


#ChipQC
#library(ChIPQC)
#samples <- read.delim("Jorge_Arturo.Zepeda_Martinez.chipqc.txt")
#experiment = ChIPQC(samples)
#experiment
#ChIPQCreport(experiment)
PLOT.Quality.Controls.Summary = function(stat, index)
{
  par(mfrow=c(2,2))
  ss = stat[index, ]
  samples = sapply(ss$filename, find.samples.conditions, ID='samples')
  conditions =  sapply(ss$filename, find.samples.conditions, ID='conditions')
  ss$conditions = conditions
  ss = data.frame(ss, stringsAsFactors = FALSE)
  ss = ss[with(ss, order(samples, conditions)), ]
  
  cols.conditions = data.frame(unique(conditions), c(1:length(unique(conditions))))
  cols = cols.conditions[match(ss$conditions, cols.conditions[,1]), 2]
  
  barplot(as.numeric(ss$nb.UniqueRead)/10^6, col=cols, main="nb of uniquely-mapped Reads", horiz=TRUE, names.arg = ss$filename, 
          las=1, xlim = c(0, 50))
  abline(v=c(10, 20), col='red', lty=1, lwd=2.0);abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
  
  pbc = ss$nb.UniqueRmdupRead/ss$nb.UniqueRead
  barplot(pbc, col=cols, main="NRF", horiz=TRUE, names.arg = ss$filename, 
          las=1, xlim = c(0, 1))
  abline(v=c(0.9), col='red', lty=1, lwd=2.0);
  
  
  plot(ss$NSC, ss$RSC, type='n', xlab='NSC', ylab='RSC', ylim=range(ss$RSC), xlim=range(ss$NSC))
  colfunc <- colorRampPalette(c("red", "green"))
  cols = colfunc(5)
  #jj = which(stat$samples !='Input')
  #stat = stat[jj,]
  for(n in 1:nrow(ss))
  {
    points(ss$NSC[n], ss$RSC[n], type='p', col=cols[(ss$Quality[n]+3)], pch=16, cex=2.0)
    #else  points(stat$NSC.1.05.[n], stat$RSC.0.8.[n], type='p', col=cols[(stat$QualityTag..2.verylow..1.low.0.medium.1.high.2.veryhigh.[n]+3)], pch=17, cex=2.0)
    text(ss$NSC[n], ss$RSC[n], n, cex=1.6, offset=0.5, pos = 1)
  }
  abline(h=c(0.8), col='red', lty=3, lwd=2.0);abline(v=c(1.05), col='red', lty=3, lwd=2.0)
  legend('topleft', legend=rev(c('very low', 'low', 'medium', 'high', 'very high')), col=rev(cols), pch=16, bty='n', cex=1.5)
  
  plot(ss$NSC, ss$RSC, type='n',ylim=c(0, 2.0), xlim=c(1, 1.5), xlab=NA, ylab=NA, axes=FALSE)
  
  if(nrow(ss)<=40){
    legend('topleft', legend=paste(c(1:nrow(ss)), ss$filename, sep='-'), col=cols[ss$Quality+3], pch=16, bty='n', cex=1.2)
  }else{
    legend('topleft', legend=paste(c(1:40), ss$filename[1:40], sep='-'), col=cols[ss$Quality+3][1:40], pch=16, bty='n', cex=1.)
    legend('topright', legend=paste(c(41:nrow(ss)), ss$filename[41:nrow(ss)], sep='-'), col=cols[ss$Quality+3][41:nrow(ss)], pch=16, bty='n', cex=1.2)
  }
  
}

Comparison.overlapping.peaks = function(design.matrix, peaks.list, toCompare="factor.condition", pval=10, PLOT.p10 = FALSE, qval=0.1) 
{
  if(toCompare == "factor"){
    cat("DB analysis for each factor across condition \n")
    cat("Should merge peaks for all replicates \n")
    
  }
  if(toCompare=="factor.condition"){
    cat("Quality control by compare peaks from replicates \n")
    sels2compare = unique(design.matrix$factor.condition)
  }

  for(nn in sels2compare)
  {
    #nn = "H3K27me3_AN312"
    #nn = "Eed_ko_AK119ub"
    cat(nn, '\n')
    kk = which(design.matrix$factor.condition==nn)
    #if(DB.Analysis){
    #  kk = which(ff$sample==nn)
    #}else{
     
    #}
    
    if(length(kk)>1)
    {
      peaks = c()
      peaks10 = c()
      peaknames = design.matrix$file.name[kk]
      for(k in kk) 
      {
        p = peaks.list[[k]]
        #p = readPeakFile(ff$file.path[k], as = "GRanges");
        #eval(parse(text = paste0("p = pp.", k)));
        p10 <- p[mcols(p)[,"X.log10.pvalue."] > pval];
        p = reduce(p); 
        p10 = reduce(p10);
        peaks= c(peaks, p);
        peaks10= c(peaks10, p10);
        #peaknames = c(peaknames, ff$file.name[k])
        #if(DB.Analysis){
        #  peaknames = c(peaknames, ff$sample.condition[k])
        #}else{
        #}
        #test = basename(macs2.files[k])
        #test = find.sample.names(test, n2=2)
        #test = unlist(strsplit(as.character(test), '_'))
        #test =test[length(test)]
        #nb.reads = bamstatistics[grep(test, bamstatistics[,2]), 5]
        #peaknames = c(peaknames, paste(ff$file.name[k], test,  sep='_', collapse = '_'))
      }
      
      ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=peaknames, connectedPeaks="keepAll", main=nn)
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      if(PLOT.p10){
        ol.peaks <- makeVennDiagram(peaks10, NameOfPeaks=paste(peaknames, '_p10', sep=''), connectedPeaks="keepAll", main=nn)
        v <- venn_cnt2venn(ol.peaks$vennCounts)
        try(plot(v))
      }
    }
  }
  
}

find.sample.names = function(x, n1=1, n2=3)
{
  temp = unlist(strsplit(as.character(x), '_'));
  kk = which(temp=='macs2')
  return(paste(temp[c(n1:n2)], sep='', collapse = '_'))
}


find.samples.conditions = function(x, ID='samples')
{
  if(ID=='samples') 
  { 
    j = 2
    temp = unlist(strsplit(as.character(x), '_'));
    return(temp[j])
  }
  if(ID=='conditions') 
  {
    j = 1;
    temp = unlist(strsplit(as.character(x), '_'));
    temp = temp[j]
    temp = unlist(strsplit(as.character(temp), '[.]'));
    return(temp[length(temp)])
    #kk = which(temp=='macs2')
    #return(paste(temp[c(n1:(kk-n2))], sep='', collapse = '_'))
  }
}

merge.peaks.macs2 = function(peak.list){
  
}

Assess.ChIPseq.Quality.4DB = function(dds, batch = FALSE) ### Input is DEseq2 object from read counts within peak union across condition
{
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr"); 
  library("ggplot2")
  require(lattice);
  #require(ggplot2)
  cc.uniq = unique(conds);
  cols = match(conds, cc.uniq)
  
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  xx = fpm;
  
  ### boxplot (distribution) of gene expression for each sample
  par(mfrow=c(1,1))
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(15,3,2,0.2), tcl = -0.3)
  for(n in 1:ncol(xx)){
    kk = which(xx[,n]>0);
    if(n==1) boxplot(log2(xx[kk, n]), horizontal = FALSE, las=2, at=(n), ylim=c(-6, 23), xlim=c(0, (ncol(xx)+1)), names=as.character(colnames(xx)[n]),
                     las=1, width=0.6, ylab='log2(fpm)', col=cols[n], main="Distribution of normalized signals (cpm)")
    else boxplot(log2(xx[kk, n]), horizontal = FALSE, las=1, add=TRUE, at=(n), names=colnames(xx)[n], width=0.6, col=cols[n])
    mtext(colnames(xx)[n], side = 1, at=n, las=2)
  }
  
  ### pairwise correlations for fpm
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  library(corrplot)
  col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                             "cyan", "#007FFF", "blue","#00007F"))
  xx = as.matrix(fpm)
  xx[which(xx==0)] = NA
  M <- cor(xx, use = "na.or.complete")
  #corrplot(M, method="circle", type = 'upper', order="hclust")
  corrplot(M, method="ellipse", order="hclust", tl.cex=1.2, cl.cex=0.7, tl.col="black", addrect=ceiling(ncol(xx)/2), col=col1(100), rect.col=c('green'), rect.lwd=2.0)
  
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr")
  library("ggplot2")
  #colnames(fpm) = paste(colnames(fpm), '.fpm', sep='')
  #if(nrow(dds)<1000) {vsd =  rlog(dds, blind = FALSE)
  #}else{vsd <- vst(dds, blind = FALSE)}
  
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    #as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  colnames(df)[1:2] <- c("x", "y")  
  vsd.transform=ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  #print(vsd.transform)
  sampleDists <- dist(t(assay(vsd)))
  #sampleDists
  #rld <- rlog(dds, blind=FALSE)
  sampleDistMatrix <- as.matrix( sampleDists)
  #rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  #plot(distClustering)
  
  if(batch){
    pca=plotPCA(vsd, intgroup = c('batch', 'cell'), returnData=FALSE)
  }else{
    pca=plotPCA(vsd, intgroup = c('conds'), returnData=FALSE)
  }
  plot(pca);
  
  Show.sample.names.PCA.Clusters = TRUE
  if(Show.sample.names.PCA.Clusters){
    pca2save = as.data.frame(plotPCA(vsd, intgroup =c('conds'), returnData = TRUE))
    #p = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape=time)) + geom_point(size=3)
    #p + geom_text(hjust = 0.5, nudge_y = 0.1, size=2.5) 
    ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=conds)) + geom_point(size=3) +
      geom_text(hjust = 0.7, nudge_y = 1, size=2.5)  
    plot(ggp);
  }
  
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
  {
    usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(0, 1, 0, 1)) 
    r <- abs(cor(x, y)) 
    txt <- format(c(r, 0.123456789), digits=digits)[1] 
    txt <- paste(prefix, txt, sep="") 
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
    
    test <- cor.test(x,y) 
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " ")) 
    
    text(0.5, 0.5, txt, cex = cex * r) 
    text(.8, .8, Signif, cex=cex, col=2) 
  }
  
  panel.fitting = function (x, y, bg = NA, pch = par("pch"), cex = 0.01, col='black') 
  {
    #x = yy[,1];y=yy[,2];
    #kk = which(x>0 & y>0); x=x[kk];y=y[kk]
    lims = range(c(x, y), na.rm = TRUE)
    points(x, y, pch = 1, col = col, cex = cex, xlim=lims, ylim=lims)
    abline(0, 1, lwd=1.0, col='red')
    R = cor(x, y, use="na.or.complete", method='pearson')
    text(lims[2]*0.2, lims[2]*0.9, paste0('R = ', signif(R, d=2)), cex=0.5, col='red')
    jj = which(!is.na(x) & !is.na(y))
    fit = lm(y[jj] ~ x[jj])
    #slope=summary(fit)$coefficients[1]
    slope = fit$coefficients[2]
    intercept = fit$coefficients[1]
    pval=summary(fit)$coefficients[4]
    abline(intercept, slope, lwd=1.2, col='darkblue', lty=3)
    #text(lims[2]*0.1, lims[2]*0.7, paste0('slop = ', signif(slope, d=2)), cex=1., col='blue')
    #text(lims[2]*0.1, lims[2]*0.6, paste0('pval = ', signif(pval, d=2)), cex=1., col='blue')
    #ok <- is.finite(x) & is.finite(y)
    #if (any(ok)) 
    #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
    #        col = col.smooth, ...)
  }
  
  ###  pairwise correlation and fitting (to chek if there is batch effect)
  yy = as.matrix(fpm)
  #yy = as.matrix(assay(vsd))
  #yy = as.matrix(xx[, c(1:3)])
  yy[which(yy==0)] = NA;
  yy = log2(yy)
  
  cc.partA = c("AN312", "515D10", "515D10H3")
  cc.partB = c("AN312", "515D10", "D10A8", "D10D8")
  cc.partC = c("924E12", "515D10", "E12F01")
  
  kk = match(dds$conds, cc.partA); ii = which(!is.na(kk));
  if(length(unique(dds$conds[ii]))>2) pairs(yy[, ii], lower.panel=NULL, upper.panel=panel.fitting)
  
  kk = match(dds$conds, cc.partB); ii = which(!is.na(kk));
  if(length(unique(dds$conds[ii]))>2) pairs(yy[, ii], lower.panel=NULL, upper.panel=panel.fitting)
  
  kk = match(dds$conds, cc.partC); ii = which(!is.na(kk));
  if(length(unique(dds$conds[ii]))>2) pairs(yy[, ii], lower.panel=NULL, upper.panel=panel.fitting)
  
}

