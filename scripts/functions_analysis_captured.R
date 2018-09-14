calcNormFactors.for.caputred.using.csaw = function(dd, cutoff.average.counts = 20, method = "DESeq2")
{
  # dd = binned.chipseq;
  if(method == "edgeR"){
    abundances <- aveLogCPM(asDGEList(dd), prior.count = 1)
    summary(abundances)
    
    #cutoff.average.counts = 100;
    keep <- abundances > aveLogCPM(cutoff.average.counts, lib.size=mean(dd$totals))
    #keep = filterWindows(dd, type = "proportion")$filter > 0.3
    cat(summary(keep))
    #keep.simple <- abundances > -1
    dd.filtered <- dd[keep,]
    
    #lib.sizes.full = dd.filtered$totals;
    dd.filtered <- normOffsets(dd.filtered, se.out = dd.filtered)
    
    #filtered.data$norm.factors = norms
    norms = as.data.frame(colData(dd.filtered));
    totals = apply(assay(dd.filtered), 2, sum)
    norms$total.filtered = totals
    
    test.by.plot = FALSE
    if(test.by.plot){
      par(mfrow=c(1, 1))
      plot(norms, dd.filtered$totals/10^6, log='xy', xlim=c(0.1, 2))
      abline(0, 1, lwd=2.0, col='red')
      
      # to visulize the normalization
      nb.for.test = 2
      par(mfrow=c(1, nb.for.test), mar=c(5, 4, 2, 1.5))
      adj.counts <- cpm(asDGEList(dd.filtered), log=TRUE)
      normfacs <- dd.filtered$norm.factors
      
      for (i in 1:nb.for.test) {
        cur.x <- adj.counts[,1]
        cur.y <- adj.counts[,1+5]
        smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                      xlab="A", ylab="M", main=paste("1 vs", i+1))
        all.dist <- diff(log2(normfacs[c(i+1, 1)]))
        abline(h=all.dist, col="red")
      }
      
      par(mfrow=c(1, 1))
      plot(adj.counts[c(1:10000), c(1, 5)], cex=0.2, log='')
      abline(0, 1, col='red', lwd=2.0)
      
    }
    
  }
  
  if(method == "DESeq2"){
    require('DESeq2');
    countData = assay(dd);
    condition <- factor(rep("A", ncol(countData)))
    
    dds <- DESeqDataSetFromMatrix(assay(dd), data.frame(condition),  design = ~ 1 )
    dds <- dds[ rowSums(counts(dds)) > cutoff.average.counts, ]
    
    cat("after filtering--", nrow(dds), "\n")
    
    dds <- estimateSizeFactors(dds)
    
    norms = data.frame(as.data.frame(colData(dd)), size.factors = sizeFactors(dds))
    
    xx = norms$totals
    yy = norms$size.factors
    xx = xx/xx[1]; yy = yy/yy[1]
    xx = 1/xx; yy = 1/yy;
    plot(xx, yy, log = 'xy', xlab='library size', ylab="size.factors")
    abline(0, 1, col='red', lwd=2.0)
    text(xx, yy, basename(norms$bam.files), cex=0.6, pos = 1, offset = 0.5)
    #fpm = fpm(dds, robust = TRUE)
    #vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  }
  
  return(norms)
  
}