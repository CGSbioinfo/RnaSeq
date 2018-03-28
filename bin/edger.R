#!/usr/local/bin/Rscript

suppressMessages(library(edgeR))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
source('bin/edger_functions.R')

####################################################
# In the multipleComparisons function call, changed
# min.counts to min.cpm to match argument file
####################################################

# Read arguments from standard file
arguments_file = commandArgs(TRUE)[1]
arguments_file=read.csv(arguments_file, header=FALSE)

read_args_line=function(x){
  x=gsub(' *$', '',x)
  x=gsub('.*=', '',x)
  x=gsub('^ ','',x)
  return(x)
}
indir=read_args_line(grep('^indir =', arguments_file$V1, value=TRUE))
outdir=read_args_line(grep('^outdir =', arguments_file$V1, value=TRUE))
sample_info=read_args_line(grep('^sample_info =', arguments_file$V1, value=TRUE))
comparisons=read_args_line(grep('^comparisons =', arguments_file$V1, value=TRUE))
min.cpm=as.numeric(read_args_line(grep('^min.cpm =', arguments_file$V1, value=TRUE)))
min.nsamples=as.numeric(read_args_line(grep('^min.nsamples =', arguments_file$V1, value=TRUE)))
design=read_args_line(grep('^design =', arguments_file$V1, value=TRUE))
gtf.file=read_args_line(grep('^gtfFile =', arguments_file$V1, value=TRUE))
gc.length.table = read_args_line(grep('^gc.length.table =', arguments_file$V1, value= TRUE))
gc.length.correction = read_args_line(grep('^gc.length.correction =', arguments_file$V1, value= TRUE))

dir.create(outdir, recursive=TRUE, showWarnings = FALSE)


# Reading counted reads files #
#-----------------------------#

files=list.files(indir, pattern = "_count.txt", full.names=TRUE)

data=read.table(files[1], row.names=1)
colnames(data)=files[1]

for (i in 2:length(files)){
  name=files[i]
  assign(name, read.table(name))
  sample=eval(as.name(name))
  index=match(rownames(data), sample[,1])
  data=cbind(data, sample[index,2])
  colnames(data)[i]=name
}
colnames(data)=gsub('_count.txt','',colnames(data))
colnames(data)=gsub('.*/','',colnames(data))

# Keep reads mapping to features
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),]



# Getting information about the groups #
#--------------------------------------#
sample_info = read.csv(sample_info)
sample_info=sample_info[order(sample_info$Group),]
group <- c(as.character(sample_info$Group))


# Matching the data matrix to the sample_info table #
#---------------------------------------------------#
data=data[,match(as.character(sample_info$SampleID),colnames(data))]

# Remove undetected genes in all samples #
#----------------------------------------#
#data_all_samples=data[!rowSums(data)==0,]

# Creating a DGE object to analyse data with edgeR #
#--------------------------------------------------#
#dge <- DGEList(counts=data_all_samples, group=group)

# Filter low counts genes #
#-------------------------#
#keep <- rowSums(cpm(dge)>min.cpm) >= min.nsamples
#dge.all.samples <- dge[keep,]
#dge.all.samples$samples$lib.size <- colSums(dge.all.samples$counts)

# Normalizing dge #
#-----------------#
#dge.all.samples <- calcNormFactors(dge.all.samples)

# Exploring data  #
#-----------------#
#mycoldf=cbind(unique(as.character(dge.all.samples$samples$group)),1:length(unique(as.character(dge.all.samples$samples$group))))
#mycol=c()
#for (i in 1:length(as.character(dge.all.samples$samples$group))){
#  mycol=c(mycol,mycoldf[which(mycoldf[,1]==as.character(dge.all.samples$samples$group)[i]),2])
#}
#pdf(paste0(outdir,'/MDSPlot_allSamples_bcv', ".pdf"))
#plotMDS(dge.all.samples, col=mycol, method="bcv")# , xlim=c(-1.4,1.4),ylim=c(-.7,.7))# 
#dev.off()
#pdf(paste0(outdir,'/MDSPlot_allSamples_logFc', ".pdf"))
#plotMDS(dge.all.samples, col=mycol, method="logFC")# , xlim=c(-1.4,1.4),ylim=c(-.7,.7))# 
#dev.off()

#y = cpm(dge.all.samples,prior.count = 1, log=TRUE)
#pdf(paste0(outdir,'/Heatmap_allSamples', ".pdf"), width=9,height=9)
#heatmap.2(cor(y),scale=c('none'), density.info='density', trace='none', margins=c(8,8), cex.lab=0.8)
#dev.off()

#pdf(paste0(outdir,'/Libsize', ".pdf"), width=9)
#par(mar=c(14.8,4.1,4.1,2.1))
#barplot(dge.all.samples$samples$lib.size*1e-6, names=colnames(data), las=2, ylab="Library size (millions)") # library sizes
#dev.off()


# Doing Differential Gene Expression  #
#-------------------------------------#
comparisons=read.csv(comparisons, header=TRUE)
if (design=='pairedSamples'){
  pairedDesign=TRUE
} else if(design=='non-pairedSamples'){
  pairedDesign=FALSE
}

if (length(sample_info$Group)==length(unique(sample_info$Group))){
  print('No replicate analysis')
  multipleComparison_no_replicates(data,comparisons,pairedDesign, min.cpm, min.nsamples, gtf.file)
} else {
  multipleComparison(data,comparisons,pairedDesign, min.cpm, min.nsamples, gtf.file, gc.length.correction, gc.length.table)
}
