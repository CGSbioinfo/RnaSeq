suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
#suppressMessages(library(GGally))
suppressMessages(library(ggplot2))
#suppressMessages(library(vsn))


multipleComparison=function(data,comparisons,pairedDesign, gtf.file){
  # Initialize vector to record number of DE genes
  significant01=c()
  significant05=c()
  
  # Get annotation
  GTF <- import.gff(gtf.file, format="gtf", feature.type="gene")
  df=GTF$gene_name
  names(df)=GTF$gene_id
  
  for (i in 1:nrow(comparisons)){
    print(i)
    # Subset data and setup environment
    temp_comparison=c(as.character(comparisons[i,1]), as.character(comparisons[i,2]))
    temp_sample_info=sample_info[which(sample_info$Group %in% temp_comparison),]
    temp_sample_info=droplevels(temp_sample_info)
    temp_data=data[,which(colnames(data)%in%temp_sample_info$SampleID)]
    group=c(as.character(temp_sample_info$Group))
    newd=paste0(temp_comparison, collapse = "__VS__")
    dir.create(paste0(outdir,'/',newd), showWarnings=FALSE)
    
    if (pairedDesign==TRUE){
      print('paired')
      temp_colData<-data.frame(Group=temp_sample_info$Group, Sibship=temp_sample_info$Sibship)
      rownames(temp_colData)=temp_sample_info$SampleID
      dds <- DESeqDataSetFromMatrix(countData= temp_data, colData=temp_colData, design=~Group)
      design(dds)<-formula(~Sibship+Group)
      dds$Group=relevel(dds$Group, ref=temp_comparison[2])
    } else if (pairedDesign==FALSE) {
      print('non-paired')
      temp_colData<-data.frame(Group=temp_sample_info$Group)
      rownames(temp_colData)=temp_sample_info$SampleID
      dds <- DESeqDataSetFromMatrix(countData= temp_data, colData=temp_colData, design= ~Group)
      dds$Group=relevel(dds$Group, ref=temp_comparison[2])
    }
    print(colData(dds))
  
    
    #  Prefiltering #
    #---------------#
    dds <- dds[rowSums(counts(dds)) > 1, ]
    print(dds)
    dds <- DESeq(dds)

    # Normalizing dds #
    #-----------------#
    #dds <- estimateSizeFactors(dds)
    #pdf(paste0(outdir,'/',newd,'/size_factors_', newd, ".pdf"))
    #plot(sizeFactors(dds), colSums(counts(dds)))
    #abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
    #dev.off()
    
    # Exploring data  #
    #-----------------#
    rld <- rlog(dds, blind=FALSE)

    # Heatmap genes 
    pdf(paste0(outdir,'/',newd,'/Heatmap_',newd, ".pdf"), width=9,height=9)
    heatmap.2(cor(assay(rld)), margins=c(14,14), scale=c('none'), density.info='density', trace='none')
    dev.off()

    #y=as.matrix(dist(t(assay(rld))))
    #pdf(paste0(outdir,'/',newd,'/Heatmap_dist_', newd,".pdf"), width=9,height=9)
    #heatmap.2(y, margins=c(14,14), scale=c('none'), density.info='density', trace='none')
    #dev.off()

    # PCA
    pca_data=plotPCA(rld, intgroup=c('Group'), returnData=TRUE)
    percentVar=round(100*attr(pca_data,'percentVar'))
    print(percentVar)
    print(pca_data)
    #pdf(paste0(outdir,'/',newd,'/PCA_', newd,".pdf"), width=9,height=9)
    p=ggplot(pca_data, aes(PC1,PC2, color=Group, label=rownames(pca_data))) + geom_point() + 
     geom_text(show_guide=F) + xlab(paste0('PC1: ', percentVar[1], '% variance'))  + ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
     theme(panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.3,linetype=2), panel.grid.minor=element_line(colour='grey',size=.3,linetype=2)) + xlim(range(pca_data$PC1)[1]-3, range(pca_data$PC1)[2]+3)
    ggsave(p,filename=paste0(outdir,'/',newd,'/PCA_', newd,".pdf"), width=9,height=9)
    #dev.off()
    
    #scatter plot
    png(paste0(outdir,'/',newd,'/ScatterPlot_', newd, ".png"), width = 1360, height = 1360)
    pairs.panels(assay(rld), smooth=FALSE)
    dev.off()
    
    # Design
    #if (pairedDesign==TRUE){
    #  sibship=factor(temp_sample_info$Sibship)
    #  treatment=factor(temp_sample_info$Group)
    #  design(dds) <- ~sibship + treatment # this is not tested !!
    #} else {
    #  design(dds) <- ~Group
    #}
    
    # DE
    #dds <- DESeq(dds)
    res <- results(dds)
    res <- res[order(res$padj),]
    pdf(paste0(outdir, '/',newd,'/MAplot_', newd, ".pdf"))
    plotMA(res)
    dev.off()

    #resMLE=results(dds,addMLE=TRUE)

    #head(res)
    res <- as.data.frame(res)

    detags <- rownames(res)
    res = cbind(res,counts(dds,normalized=TRUE)[detags,]) # chech individual cpm values for 
    res=cbind(GeneName=df[rownames(res)], res)
    write.csv(res,paste0(outdir, '/',newd,'/Results_', newd, ".csv"))
    
    significant01=c(significant01, sum(res$padj <0.01, na.rm=TRUE))
    names(significant01)[i]=newd
    significant05=c(significant05, sum(res$padj <0.05, na.rm=TRUE))
    names(significant05)[i]=newd

    # Dispersion
    pdf(paste0(outdir, '/',newd,'/Dispersion_', newd, ".pdf"))
    plotDispEsts(dds)
    dev.off()
  }
  summ=cbind(significant01,significant05)
  colnames(summ)=c('p<0.01',"p<0.05")
  write.csv(summ,paste0(outdir,"/significant_summ.csv"),quote=F)  
}

########### scatterplots
panel.cor.scale <- function(x, y, digits=2, prefix="")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex * abs(r))
}
panel.cor <- function(x, y, digits=4, prefix="")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 2)
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE)
{if (smooth ){
  if (scale) {
    pairs(x,diag.panel=panel.hist,lower.panel=panel.smooth)
  }
  else {pairs(x,diag.panel=panel.hist,lower.panel=panel.smooth)
  } #else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
}
else #smooth is not true
{ if (scale) {pairs(x,diag.panel=panel.hist)
} else {pairs(x,diag.panel=panel.hist, cex.labels = 1.3) }
} #end of else (smooth)
} #end of function

############## PCA modified
plotPCA_4pcs=function (object, ...) 
{
    .local <- function (object, intgroup = "condition", ntop = 500, 
        returnData = TRUE) 
    {
        group=colData(object)$Group
        if (class(object)=="DESeqTransform"){
	object<-assay(object)	
	}
        rv <- rowVars(object)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
            length(rv)))]
        pca <- prcomp(t(object[select, ]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], Group=group, name = colnames(object))
        if (returnData) {
            attr(d, "percentVar") <- percentVar[1:4]
            return(d)
        }
    }
    .local(object, ...)
}

pca_plot_single_label=function(x, pcs, col, shape, percentVar){
  plus_margin=diff(range(x[,pcs[1]]))/2
  plus_marginy=diff(range(x[,pcs[2]]))/10
  plot=ggplot(as.data.frame(x), aes_string(x=pcs[1], y=pcs[2], col=col, shape=shape)) +
    geom_point(size=.8) + geom_text(label=as.character(rownames(x)),show.legend = FALSE, vjust=1.3, size=4) + 
    xlab(paste0(pcs[1],' ', percentVar[pcs[1]], "%")) + ylab(paste0(pcs[2], ' ', percentVar[pcs[2]], "%")) +
    ggtitle(paste0(pcs[1],' vs ',pcs[2])) +
    theme(plot.margin= unit(c(1, 1, 1, 1), "cm"),
          panel.border = element_rect(fill=NA),
          panel.grid.major =element_line(colour = 'grey', linetype='dashed'),
          legend.text = element_text(size=15),
          legend.title= element_text(size=15),
          plot.title = element_text(size=15, face='bold'),
          axis.text = element_text(size=15),
          axis.title = element_text(size=15))  +
    guides(colour = guide_legend(override.aes = list(size=2.2)))  +
    xlim(min(x[,pcs[1]])-plus_margin,max(x[,pcs[1]])+plus_margin) +
    ylim(min(x[,pcs[2]])-plus_marginy, max(x[,pcs[2]])+plus_marginy) +
    scale_shape_manual(values=1:nrow(x))
  return(plot)
}

pca_plot_single=function(x, pcs, col, percentVar){
  plus_margin=diff(range(x[,pcs[1]]))/10
  plus_marginy=diff(range(x[,pcs[2]]))/15
  plot=ggplot(as.data.frame(x), aes_string(x=pcs[1], y=pcs[2], col=col)) +
    geom_point(size=2) +
    xlab(paste0(pcs[1],' ', percentVar[pcs[1]], "%")) + ylab(paste0(pcs[2], ' ', percentVar[pcs[2]], "%")) +
    ggtitle(paste0(pcs[1],' vs ',pcs[2])) +
    theme(plot.margin= unit(c(1, 1, 1, 1), "cm"),
          panel.border = element_rect(fill=NA),
          panel.grid.major =element_line(colour = 'grey', linetype='dashed'),
          legend.text = element_text(size=15),
          legend.title= element_text(size=15),
          plot.title = element_text(size=15, face='bold'),
          axis.text = element_text(size=15),
          axis.title = element_text(size=15))  +
    guides(colour = guide_legend(override.aes = list(size=2.2)))  +
    xlim(min(x[,pcs[1]])-plus_margin,max(x[,pcs[1]])+plus_margin) +
    ylim(min(x[,pcs[2]])-plus_marginy, max(x[,pcs[2]])+plus_marginy)
  return(plot)
}

pca_plot_pairs=function(x, plot_name, column_col){
  pca=prcomp(t(exprs(x)))
  percentVar<-pca$sdev^2/sum(pca$sdev^2)
  percentVar=round(percentVar*100)
  names(percentVar)=paste0('PC',1:length(percentVar))
  pca_points<-pca$x
  pca_points<-data.frame(SampleName=row.names(pca_points),pData(x),pca_points)
  
  pcs=c('PC1','PC2')
  plot1=pca_plot_single(pca_points,pcs,column_col, percentVar=percentVar)
  
  pcs=c('PC1','PC3')
  plot2=pca_plot_single(pca_points,pcs,column_col, percentVar=percentVar)
  
  pcs=c('PC1','PC4')
  plot3=pca_plot_single(pca_points,pcs,column_col, percentVar=percentVar)
  
  pcs=c('PC2','PC3')
  plot4=pca_plot_single(pca_points,pcs,column_col, percentVar=percentVar)
  
  pcs=c('PC2','PC4')
  plot5=pca_plot_single(pca_points,pcs,column_col, percentVar=percentVar)
  
  pcs=c('PC3','PC4')
  plot6=pca_plot_single(pca_points,pcs,column_col, percentVar=percentVar)
  
  png(plot_name, height=880, width=1680)
  grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,ncol=3)
  dev.off()
}

################# clean y
cleanY = function(y, mod1, svs) {
  X = cbind(mod1, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod1)
  #dat_transformed=dat0 - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

