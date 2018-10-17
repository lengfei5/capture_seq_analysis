############# 
##### Collection of functions for RNA-seq data
############# 
Check.RNAseq.Quality = function(read.count, design.matrix, keep.All = FALSE, norms=NULL, Threshold.read.counts=10)
{
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr"); 
  library("ggplot2")
  
  ##########################################
  # order the samples according to the second column except the sample name 
  ##########################################
  o1 = order(design.matrix[,2])
  design.matrix = design.matrix[o1, ]
  read.count = read.count[, o1];
  if(!is.null(norms)) norms = norms[o1];
  
  #load(file=paste0('Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  # kk = grep('Ab', colnames(bb))
  if(ncol(design.matrix)>2){cc = apply(design.matrix[, -1], 1, paste0, collapse="_")
  }else{cc = design.matrix[, -1]}
  #o1 = order(cc)
  #read.count = read.count[o1,]
  #cc = cc[o1]
  raw = as.matrix(read.count)
  #xx = raw
  dim(raw)
  raw[which(is.na(raw))] = 0
  xx = raw;
  
  par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,16,2,0.8)+0.1, tcl = -0.3)
  par(mfrow=c(1,1))
  
  total = apply(raw, 2, sum)
  cc.uniq = unique(cc);
  cols = match(cc, cc.uniq)
  #cols = (c(1:length(unique(cc)))-1)%/%3+1
  barplot(total/10^6, horiz = TRUE, names.arg = colnames(raw), las=1, col = cols, main='Total nb of reads quantified for features', xlab='number of reads (Million)')
  abline(v=c(1, 2, seq(5, 20, by=5)), col='red', lty=1, lwd=2.0);#abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
  xx = raw
  xx[which(xx==0)] = NA
  
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
  boxplot(log10(xx), las=3, col=cols, ylab='log10(nb of reads)', main='Read distribution for features')
  
  ### make DESeq object using read counts and design matrix
  countData = ceiling(raw)
  conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
  eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
  #dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
  
  if(!keep.All){
    dds <- dds[ rowSums(counts(dds)) > Threshold.read.counts, ]
  }
  
  if(is.null(norms)){
    dds <- estimateSizeFactors(dds)
  }else{
    sizeFactors(dds) = norms;
  }
  
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  ### boxplot (distributions of fpm) for all samples
  xx = as.matrix(fpm)
  par(mfrow=c(1,1))
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
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
  
  ### count transformation using vsd
  library("dplyr")
  library("ggplot2")
  dds <- estimateSizeFactors(dds)
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    #as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  ## illustration of count transformation using two samples
  colnames(df)[1:2] <- c("x", "y")  
  vsd.transform=ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  print(vsd.transform)
  
  ### heatmap clustering samples using sample distance  
  sampleDists <- dist(t(assay(vsd)))
  #sampleDists
  #rld <- rlog(dds, blind=FALSE)
  sampleDistMatrix <- as.matrix( sampleDists )
  #rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  
  ## heatmap clustering samples using poisson distance 
  library("PoiClaClu")
  poisd <- PoissonDistance(t(counts(dds)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- colnames(dds)
  colnames(samplePoisDistMatrix) <- NULL
  pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors)
  
  ### clustering samples using PCA analysis
  pca=plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = FALSE)
  print(pca)
  
  Show.sample.names.PCA.Clusters = FALSE 
  if(Show.sample.names.PCA.Clusters){
    pca2save = as.data.frame(plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = TRUE))
    #p = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape=time)) + geom_point(size=3)
    #p + geom_text(hjust = 0.5, nudge_y = 0.1, size=2.5) 
    ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape=time)) + geom_point(size=3) +
      geom_text(hjust = 0.7, nudge_y = 1, size=2.5)  
    plot(ggp);
  }
  
  ###  pairwise correlation and fitting (to chek if there is batch effect)
  if(ncol(fpm)<20)
  {
    yy = as.matrix(fpm)
    yy[which(yy==0)] = NA;
    yy = log2(yy)
    pairs(yy, lower.panel=NULL, upper.panel=panel.fitting)
  }
  
  #ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  #  geom_point(size=3) +
  #  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #  coord_fixed()
  #write.table(processed, file=paste0('Table_quanseq_WT_tbx_negative_positive', version.analysis, '.txt'), sep='\t', quote=FALSE, col.names = TRUE, row.names = FALSE)
  #kk = grep('fpm', colnames(process
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

panel.fitting = function (x, y, bg = NA, pch = par("pch"), cex = 0.5, col='black') 
{
  #x = yy[,1];y=yy[,2];
  #kk = which(x>0 & y>0); x=x[kk];y=y[kk]
  lims = range(c(x, y), na.rm = TRUE)
  points(x, y, pch = 1, col = col, cex = cex, xlim=lims, ylim=lims)
  abline(0, 1, lwd=1.5, col='red')
  R = cor(x, y, use="na.or.complete", method='pearson')
  text(lims[2]*0.2, lims[2]*0.9, paste0('R = ', signif(R, d=2)), cex=1., col='red')
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

