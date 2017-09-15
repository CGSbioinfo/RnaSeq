suppressMessages(library(edgeR))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(gridExtra))
suppressMessages(library(cqn))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

multipleComparison=function(data,comparisons,pairedDesign, min.count, min.nsamples, gtf.file, gc.length.correction=FALSE, gc.length.table=NULL){
  # Initialize vector to record number of DE genes
  significant01=c()
  significant05=c()
  significant01_gc=c()
  significant05_gc=c()
  overlap_up_0.05=c()
  overlap_down_0.05=c()
  overlap_0.05=c()
  numbers_of_genes=c()
  
  
  # Get annotation
  GTF <- import.gff(gtf.file, format="gtf")
  df=GTF$gene_name
  names(df)=GTF$gene_id
  
  for (i in 1:nrow(comparisons)){
    
    # Subset data and setup environment
    temp_comparison=c(as.character(comparisons[i,1]), as.character(comparisons[i,2]))
    
    cat(paste0("Performing comparison: ", temp_comparison, collapse = "__VS__"))
    cat("\n")
    temp_sample_info=sample_info[which(sample_info$Group %in% temp_comparison),]
    temp_sample_info=droplevels(temp_sample_info)
    temp_sample_info$Group=factor(temp_sample_info$Group, levels=c(as.character(comparisons[i,1]), as.character(comparisons[i,2])))
    
    temp_data=data[,which(colnames(data)%in%temp_sample_info$SampleID)]
    temp_data=temp_data[!rowSums(temp_data)==0,]
    
    group=group <- temp_sample_info$Group
    newd=paste0(temp_comparison, collapse = "__VS__")
    dir.create(paste0(outdir,'/',newd), showWarnings=FALSE)
    
    dge <- DGEList(counts=temp_data, group=group)
    
    # Filtering dge #
    #---------------#
    
    cat("Filtering reads\n")
    keep <- rowSums(cpm(dge)>min.cpm) >= min.nsamples    
    cat(paste0("Number kept ",length(keep),"\n"))
    dge <- dge[keep,]
    dge$samples$lib.size <- colSums(dge$counts)

    # Normalizing dge #
    #-----------------#
    cat("Normalising (for non gc/length corrected comparisons)\n")
    dge_norm <- calcNormFactors(dge)
    
    ##### loading gc content table
    
    # if(gc.length.correction){
    #   cat("Calculating gc and length correction\n")
    #   gc.length.data<-read.table(gc.length.table, header=T, row.names = 1)
    #   gc.length.data<-gc.length.data[row.names(dge$counts),]
    # 
    #   # NAs in the length are set to 1 (length needs to be at least 1)
    #   gc.length.data$length<-log2(gc.length.data$length)
    #   if(length(gc.length.data[is.na(gc.length.data$length),]$length)>0){
    #     gc.length.data[is.na(gc.length.data$length),]$length<-1
    #   }
    #   #cat(length(gc.length.data[is.na(gc.length.data$length),]$gccontent))
    #   if(length(gc.length.data[is.na(gc.length.data$gccontent),]$gccontent)>0){
    #     gc.length.data[is.na(gc.length.data$gccontent),]$gccontent<-median(gc.length.data$gccontent, na.rm = T)
    #   }
    #   cat("Correcting gc and length bias (skipping edgeR norm)\n")
    #   cqn.subset<-cqn(dge$counts, lengths=gc.length.data$length, x=gc.length.data$gccontent, sizeFactors = dge$samples$lib.size, verbose=T)
    # }
    
    # Library size  #
    #---------------#
    cat("Plotting library size\n")
    png(paste0(outdir,'/',newd,'/Libsize_',newd, ".png"), height=750, width=750, res=120)
    par(mar=c(14.8,4.1,4.1,2.1))
    barplot(dge$samples$lib.size*1e-6, names=colnames(temp_data), las=2, ylab="Library size (millions)") # library sizes
    dev.off()

    # Exploring data  #
    #-----------------#
    
    cat("Plotting MDS\n")
    mycoldf=cbind(unique(as.character(dge$samples$group)),1:length(unique(as.character(dge$samples$group))))
    mycol=c()
    for (j in 1:length(as.character(dge$samples$group))){
      mycol=c(mycol,mycoldf[which(mycoldf[,1]==as.character(dge$samples$group)[j]),2])
    }
    png(paste0(outdir,'/',newd,'/MDS_', newd, ".png"), height=750, width=750, res=120)
    plotMDS(dge, col=mycol, method="bcv")
    dev.off()
    
    # Heatmap
    cat("Ploting heatmap\n")
    y = cpm(dge,prior.count = 1, log=TRUE)
    png(paste0(outdir,'/',newd,'/Heatmap_', newd, ".png") )
    heatmap.2(cor(y),scale=c('none'), density.info='density', trace='none', margins=c(15,15))
    dev.off()
    
    # scatter plot
    cat("Plotting scatter plot\n")
    y = cpm(dge,prior.count = 1, log=TRUE)
    png(paste0(outdir,'/',newd,'/ScatterPlot_', newd, ".png"), width = 1360, height = 1360, res=120)
    pairs.panels(y, smooth=FALSE)
    dev.off()
    
    # Dispersion
    cat("Calculating dispersion for non corrected data\n")
    if (pairedDesign==TRUE){
      sibship=factor(temp_sample_info$Sibship)
      treatment=factor(temp_sample_info$Group, levels=c(temp_comparison)) # Note that the first group in the levels vector is the baseline for the comparison
      design<-model.matrix(~sibship+treatment)
      dge_norm<-estimateGLMCommonDisp(dge_norm, design)
      dge_norm<-estimateGLMTrendedDisp(dge_norm, design)
      dge_norm<-estimateGLMTagwiseDisp(dge_norm, design)
    } else {
      dge_norm <- estimateCommonDisp(dge_norm) # Common dispersio estimates the overall BCV of the dataset, averaged over all genes # BCV is the sqroot of the common dispersion
      dge_norm <- estimateTagwiseDisp(dge_norm) # Estimaes gene specific dispersions
    }
    
    ## Plot BCV for non corrected data
    png(paste0(outdir,'/',newd,'/Dispersion_non_corrected_', newd, ".png"), height=800, width=850, res=120)
    plotBCV(dge_norm, cex=0.4, pch=20)
    dev.off()
    
    # Differential Expression for non corrected data
    cat("Comparing non corrected data\n")
    if (pairedDesign==TRUE){
      ### model for non corrected data
      fit<-glmFit(dge_norm,design)
      lrt<-glmLRT(fit, coef=ncol(design))
      results=topTags(lrt, n=dim(lrt)[1])[[1]]
    } else {
      et <- exactTest(dge_norm, pair=c(temp_comparison[1], temp_comparison[2])) # Note that the first group listed in the pair is the baseline for the comparison-
      results=topTags(et, n=dim(et)[1])[[1]]
    }
    
    detags <- rownames(results) # chech individual cpm values for top genes
    results = cbind(results,cpm(dge,log=FALSE)[detags,]) # chech individual cpm values for top genes
    results=cbind(GeneName=df[rownames(results)], results)
    write.csv(results,paste0(outdir, '/',newd,'/Results_', newd, ".csv"))
    
    cat("plotting MA plots\n")
    if (pairedDesign==TRUE){
      de <- decideTestsDGE(lrt, p=0.05, adjust="BH") # total number of DE genes at 5% FDR
      detags <- rownames(dge_norm)[as.logical(de)]
      png(paste0(outdir, '/',newd,'/MAplot_', newd, ".png"), height=800, width=850, res=120)
      plotSmear(lrt, de.tags=detags)
      abline(h = c(-1, 1), col = "blue")
      dev.off()
    } else {
      de <- decideTestsDGE(et, p=0.05, adjust="BH") # total number of DE genes at 5% FDR
      detags <- rownames(dge_norm)[as.logical(de)]
      png(paste0(outdir, '/',newd,'/MAplot_', newd, ".png"), height=800, width=850, res=120)
      plotSmear(et, de.tags=detags, cex=0.4, pch=20)
      abline(h = c(-1, 1), col = "blue")
      dev.off()
    }
    
    if(gc.length.correction==FALSE){
      results$significant<-results$FDR<0.05 
      results$up<-results$logFC>0
      results$significant_up<-results$significant & results$up
      results$down<-results$logFC<0
      results$significant_down<-results$significant & results$down
      
      results$gc<-gc.length.data[row.names(results),]$gc
      results$length<-gc.length.data[row.names(results),]$length
      results_gc_sig<-results[,c("significant", "gc")]
      results_lenght_sig<-results[,c("significant", "length")]
      
      results_gc_sig_up<-results[,c("significant_up", "gc")]
      results_lenght_sig_up<-results[,c("significant_up", "length")]
      
      results_gc_sig_down<-results[,c("significant_down", "gc")]
      results_lenght_sig_down<-results[,c("significant_down", "length")]
      
      melt_gc<-melt(results_gc_sig,measure.vars = "gc")
      melt_lenght<-melt(results_lenght_sig, measure.vars="length")
      p_gc<-ggplot(melt_gc, aes(x=value, fill=significant))+geom_density(alpha=0.25)+ggtitle("GC content all genes")
      p_length<-ggplot(melt_lenght, aes(x=value, fill=significant))+geom_density(alpha=0.25)+ggtitle("Length all genes")
      
      melt_gc_up<-melt(results_gc_sig_up,measure.vars = "gc")
      melt_lenght_up<-melt(results_lenght_sig_up, measure.vars="length")
      p_gc_up<-ggplot(melt_gc_up, aes(x=value, fill=significant_up))+geom_density(alpha=0.25)+ggtitle("GC content up regulated genes")
      p_length_up<-ggplot(melt_lenght_up, aes(x=value, fill=significant_up))+geom_density(alpha=0.25)+ggtitle("Length up regulated gene")
      
      melt_gc_down<-melt(results_gc_sig_down,measure.vars = "gc")
      melt_lenght_down<-melt(results_lenght_sig_down, measure.vars="length")
      p_gc_down<-ggplot(melt_gc_down, aes(x=value, fill=significant_down))+geom_density(alpha=0.25)+ggtitle("GC content down regulated genes")
      p_length_down<-ggplot(melt_lenght_down, aes(x=value, fill=significant_down))+geom_density(alpha=0.25)+ggtitle("Length down regulated gene")
      
      grid= grid.arrange(p_gc,p_length, p_gc_down, p_length_down, p_gc_up, p_length_up, ncol=2, nrow=3)
      ggsave(file=paste0(outdir,"/",newd,"/gc_length_plots_",newd,".png"), plot=grid, width=24,height=12)
      
      significant01=c(significant01, table(results$FDR<0.01)["TRUE"])
      names(significant01)[i]=newd
      significant05=c(significant05, table(results$FDR<0.05)["TRUE"])
      names(significant05)[i]=newd
      numbers_of_genes=c(numbers_of_genes, dim(results)[1])
    }
    
    #####------------------------------------------#####
    ##### Processing data for gc/length correction #####
    #####------------------------------------------#####
    
    if(gc.length.correction==TRUE){
      cat("Calculating gc and length correction\n")
      gc.length.data<-read.table(gc.length.table, header=T, row.names = 1)
      gc.length.data<-gc.length.data[row.names(dge$counts),]
      
      # NAs in the length are set to 1 (length needs to be at least 1)
      gc.length.data$length<-log2(gc.length.data$length)
      if(length(gc.length.data[is.na(gc.length.data$length),]$length)>0){
        gc.length.data[is.na(gc.length.data$length),]$length<-1
      }
      #cat(length(gc.length.data[is.na(gc.length.data$length),]$gccontent))
      if(length(gc.length.data[is.na(gc.length.data$gccontent),]$gccontent)>0){
        gc.length.data[is.na(gc.length.data$gccontent),]$gccontent<-median(gc.length.data$gccontent, na.rm = T)
      }
      cat("Correcting gc and length bias (skipping edgeR norm)\n")
      cqn.subset<-cqn(dge$counts, lengths=gc.length.data$length, x=gc.length.data$gccontent, sizeFactors = dge$samples$lib.size, verbose=T)
      cat("Comparing corrected data\n")
      if(pairedDesign == TRUE){
        dge$offset<-cqn.subset$glm.offset
        dge$offset[is.na(dge$offset)]<-0
        sibship=factor(temp_sample_info$Sibship)
        treatment=factor(temp_sample_info$Group, levels=c(temp_comparison)) # Note that the first group in the levels vector is the baseline for the comparison
        design<-model.matrix(~sibship+treatment)
        dge<-estimateGLMCommonDisp(dge, design)
        dge<-estimateGLMTrendedDisp(dge, design)
        dge<-estimateGLMTagwiseDisp(dge, design)
        png(paste0(outdir,'/',newd,'/Dispersion_corrected_', newd, ".png"), height=800, width=850, res=120)
        plotBCV(dge, cex=0.4, pch=20)
        dev.off()
        fit_gc<-glmFit(dge,design)
        lrt_gc<-glmLRT(fit_gc, coef=ncol(design))
        results_gc=topTags(lrt_gc, n=dim(lrt_gc)[1])[[1]]
      }else{
        dge$offset<-cqn.subset$glm.offset
        dge$offset[is.na(dge$offset)]<-0
        treatment=factor(temp_sample_info$Group, levels=c(temp_comparison)) # Note that the first group in the levels vector is the baseline for the comparison
        design<-model.matrix(~treatment)
        dge<-estimateGLMCommonDisp(dge, design)
        dge<-estimateGLMTrendedDisp(dge, design)
        dge<-estimateGLMTagwiseDisp(dge, design)
        png(paste0(outdir,'/',newd,'/Dispersion_corrected_', newd, ".png"), height=800, width=850, res=120)
        plotBCV(dge, cex=0.4, pch=20)
        dev.off()
        fit_gc<-glmFit(dge,design)
        cat(paste0("Coef being used is ", ncol(design)))
        lrt_gc<-glmLRT(fit_gc, coef=ncol(design))
        results_gc=topTags(lrt_gc, n=dim(lrt_gc)[1])[[1]]
      }
      detags_gc <- rownames(results_gc) # chech individual cpm values for top genes
      results_gc = cbind(results_gc,cpm(dge,log=FALSE)[detags_gc,]) # chech individual cpm values for top genes
      results_gc=cbind(GeneName=df[rownames(results_gc)], results_gc)
      write.csv(results_gc,paste0(outdir, '/',newd,'/Results_gc_length_corrected_', newd, ".csv"))
    
    
      cat("processing results\n")
    
    
      
      #-------------------# 
      #Process gc/length corrected data
      
      cat("plotting MA plot gc corrected data\n")
      
      de_gc <- decideTestsDGE(lrt_gc, p=0.05, adjust="BH") # total number of DE genes at 5% FDR
      detags_gc <- rownames(dge)[as.logical(de_gc)]
      png(paste0(outdir, '/',newd,'/MAplot_gc_lenght_corrected_', newd, ".png"), height=800, width=850, res=120)
      plotSmear(lrt_gc, de.tags=detags_gc)
      abline(h = c(-1, 1), col = "blue")
      dev.off()
      
      
      cat("Plotting gc and length graph for the significant results\n")
      results$significant<-results$FDR<0.05 
      results$up<-results$logFC>0
      results$significant_up<-results$significant & results$up
      results$down<-results$logFC<0
      results$significant_down<-results$significant & results$down
      
      results$gc<-gc.length.data[row.names(results),]$gc
      results$length<-gc.length.data[row.names(results),]$length
      results_gc_sig<-results[,c("significant", "gc")]
      results_lenght_sig<-results[,c("significant", "length")]
      
      results_gc_sig_up<-results[,c("significant_up", "gc")]
      results_lenght_sig_up<-results[,c("significant_up", "length")]
      
      results_gc_sig_down<-results[,c("significant_down", "gc")]
      results_lenght_sig_down<-results[,c("significant_down", "length")]
      
      melt_gc<-melt(results_gc_sig,measure.vars = "gc")
      melt_lenght<-melt(results_lenght_sig, measure.vars="length")
      p_gc<-ggplot(melt_gc, aes(x=value, fill=significant))+geom_density(alpha=0.25)+ggtitle("No GC bias correction")
      p_length<-ggplot(melt_lenght, aes(x=value, fill=significant))+geom_density(alpha=0.25)+ggtitle("No length bias correction")
      
      melt_gc_up<-melt(results_gc_sig_up,measure.vars = "gc")
      melt_lenght_up<-melt(results_lenght_sig_up, measure.vars="length")
      p_gc_up<-ggplot(melt_gc_up, aes(x=value, fill=significant_up))+geom_density(alpha=0.25)+ggtitle("No GC bias correction up regulated genes")
      p_length_up<-ggplot(melt_lenght_up, aes(x=value, fill=significant_up))+geom_density(alpha=0.25)+ggtitle("No length bias correction up regulated gene")
      
      melt_gc_down<-melt(results_gc_sig_down,measure.vars = "gc")
      melt_lenght_down<-melt(results_lenght_sig_down, measure.vars="length")
      p_gc_down<-ggplot(melt_gc_down, aes(x=value, fill=significant_down))+geom_density(alpha=0.25)+ggtitle("No GC bias correction down regulated genes")
      p_length_down<-ggplot(melt_lenght_down, aes(x=value, fill=significant_down))+geom_density(alpha=0.25)+ggtitle("No length bias correction down regulated gene")
      
      results_gc$significant<-(results_gc$FDR<0.05)
      results_gc$gc<-gc.length.data[row.names(results_gc),]$gc
      results_gc$length<-gc.length.data[row.names(results_gc),]$length
      results_gc_sig_correct<-results_gc[,c("significant", "gc")]
      results_lenght_sig_correct<-results_gc[,c("significant", "length")]
      melt_gc_correct<-melt(results_gc_sig_correct,measure.vars = "gc")
      melt_lenght_correct<-melt(results_lenght_sig_correct, measure.vars="length")
      p_gc_correct<-ggplot(melt_gc_correct, aes(x=value, fill=significant))+geom_density(alpha=0.25)+ggtitle("GC bias corrected plot")
      p_length_correct<-ggplot(melt_lenght_correct, aes(x=value, fill=significant))+geom_density(alpha=0.25)+ggtitle("Length bias corrected plot")
      
      results_gc$up<-(results_gc$logFC>0)
      results_gc$significant_up<-results_gc$significant & results_gc$up
      results_gc$gc<-gc.length.data[row.names(results_gc),]$gc
      results_gc$length<-gc.length.data[row.names(results_gc),]$length
      results_gc_sig_correct_up<-results_gc[,c("significant_up", "gc")]
      results_lenght_sig_correct_up<-results_gc[,c("significant_up", "length")]
      melt_gc_correct_up<-melt(results_gc_sig_correct_up,measure.vars = "gc")
      melt_lenght_correct_up<-melt(results_lenght_sig_correct_up, measure.vars="length")
      p_gc_correct_up<-ggplot(melt_gc_correct_up, aes(x=value, fill=significant_up))+geom_density(alpha=0.25)+ggtitle("GC bias corrected plot up regulated genes")
      p_length_correct_up<-ggplot(melt_lenght_correct_up, aes(x=value, fill=significant_up))+geom_density(alpha=0.25)+ggtitle("Length bias corrected plot up regulated genes")
      
      results_gc$down<-(results_gc$logFC<0)
      results_gc$significant_down<-results_gc$significant & results_gc$down
      results_gc$gc<-gc.length.data[row.names(results_gc),]$gc
      results_gc$length<-gc.length.data[row.names(results_gc),]$length
      results_gc_sig_correct_down<-results_gc[,c("significant_down", "gc")]
      results_lenght_sig_correct_down<-results_gc[,c("significant_down", "length")]
      melt_gc_correct_down<-melt(results_gc_sig_correct_down,measure.vars = "gc")
      melt_lenght_correct_down<-melt(results_lenght_sig_correct_down, measure.vars="length")
      p_gc_correct_down<-ggplot(melt_gc_correct_down, aes(x=value, fill=significant_down))+geom_density(alpha=0.25)+ggtitle("GC bias corrected plot down regulated genes")
      p_length_correct_down<-ggplot(melt_lenght_correct_down, aes(x=value, fill=significant_down))+geom_density(alpha=0.25)+ggtitle("Length bias corrected plot down regulated genes")
      
      
      #grid_gc=grid.arrange(p_gc, p_gc_correct,ncol=2)
      grid_gc_length= grid.arrange(p_gc,p_gc_correct,p_length, p_length_correct, ncol=2, nrow=2)
      grid_gc_length_up= grid.arrange(p_gc_up,p_gc_correct_up,p_length_up, p_length_correct_up, ncol=2, nrow=2)
      grid_gc_length_down= grid.arrange(p_gc_down,p_gc_correct_down,p_length_down, p_length_correct_down, ncol=2, nrow=2)
      ggsave(file=paste0(outdir,"/",newd,"/gc_length_bias_",newd,".png"), plot=grid_gc_length, width=24,height=12)
      ggsave(file=paste0(outdir,"/",newd,"/gc_length_bias_",newd,"_UP.png"), plot=grid_gc_length_up, width=24,height=12)
      ggsave(file=paste0(outdir,"/",newd,"/gc_length_bias_",newd,"_DOWN.png"), plot=grid_gc_length_down, width=24,height=12)
      
      significant01=c(significant01, table(results$FDR<0.01)["TRUE"])
      names(significant01)[i]=newd
      significant05=c(significant05, table(results$FDR<0.05)["TRUE"])
      names(significant05)[i]=newd
      significant01_gc=c(significant01_gc, table(results_gc$FDR<0.01)["TRUE"])
      names(significant01_gc)[i]=newd
      significant05_gc=c(significant05_gc, table(results_gc$FDR<0.05)["TRUE"])
      names(significant05_gc)[i]=newd
      
      ### get overlap between corrected and not corrected data
      up_corrected<-results_gc[results_gc$significant_up,]
      up_not_corrected<-results[results$significant_up,]
      overlap_up_0.05=c(overlap_up_0.05,dim(inner_join(up_corrected, up_not_corrected, by="GeneName"))[1])
      down_corrected<-results_gc[results_gc$significant_down,]
      down_not_corrected<-results[results$significant_down,]
      overlap_down_0.05=c(overlap_down_0.05,dim(inner_join(down_corrected, down_not_corrected, by="GeneName"))[1])
      corrected<-results_gc[results_gc$significant,]
      not_corrected<-results[results$significant,]
      overlap_0.05=c(overlap_0.05,dim(inner_join(corrected, not_corrected, by="GeneName"))[1])
      numbers_of_genes=c(numbers_of_genes,dim(results_gc)[1])
    }


    png(paste0(outdir, '/',newd,'/pval_hist_', newd, ".png"), height=750, width=850, res=120)
    hist(results$PValue, breaks=100)
    dev.off()

  }
  if(gc.length.correction==TRUE){
    summ=cbind(significant01,significant05, significant01_gc, significant05_gc, overlap_up_0.05, overlap_down_0.05, overlap_0.05, numbers_of_genes)
    colnames(summ)=c('p<0.01 before gc/length correction',"p<0.05 before gc/length correction", 'p<0.01 after gc/length correction',"p<0.05 after gc/length correction", "Up regulated (0.05) overlap cor/non cor", "Down regulated (0.05) overlap cor/non cor", "Overlap (0.05) cor/non cor", "Total number of genes")
    write.csv(summ,paste0(outdir,"/significant_summ.csv"),quote=F)  
  }else{
    summ=cbind(significant01,significant05, numbers_of_genes)
    colnames(summ)=c('p<0.01 significant genes',"p<0.05 significant genes", "Total number of genes")
    write.csv(summ,paste0(outdir,"/significant_summ.csv"),quote=F)  
  }
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
    pairs(x,diag.panel=panel.hist,lower.panel=panel.smooth, pch=20)
  }
  else {pairs(x,diag.panel=panel.hist,lower.panel=panel.smooth, pch=20)
  } #else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
}
else #smooth is not true
{ if (scale) {pairs(x,diag.panel=panel.hist)
} else {pairs(x,diag.panel=panel.hist, cex.labels = 1.3, pch=20) }
} #end of else (smooth)
} #end of function
