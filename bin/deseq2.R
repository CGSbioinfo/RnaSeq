#!/usr/local/bin/Rscript

print('Loading libraries')
suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
#suppressMessages(library(GGally))
suppressMessages(library(ggplot2))
#suppressMessages(library(vsn))
source("bin/deseq2_functions.R")

print('Reading arguments')
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
design=read_args_line(grep('^design =', arguments_file$V1, value=TRUE))
gtf.file=read_args_line(grep('^gtfFile =', arguments_file$V1, value=TRUE))

print('Creating output directory')
dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

# Reading counted reads files #
#-----------------------------#

print('Reading counted reads files')
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
total_counts=colSums(data)
print(dim(data))
print(head(data))
# Excluding reads not falling in features, ambiguous, too low qual, multiple locations
data.notcounted=data[which(rownames(data)=='__no_feature'):dim(data)[1],]

# Subsetting including reads
print('Reads falling to features:')
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),]
total.data.counted=colSums(data)
print(dim(data))

# Getting information about the groups #
#--------------------------------------#
print('Reading group info')
sample_info = read.csv(sample_info)
sample_info=sample_info[order(sample_info$Group),]
group <- c(as.character(sample_info$Group))

# Matching the data matrix to the sample_info table #
#---------------------------------------------------#
data=data[,match(as.character(sample_info$SampleID),colnames(data))]

# Creating a DDS object to analyse data with DESEQ2 #
#--------------------------------------------------#
print('Creating a DDS object')
colData<-data.frame(Group=sample_info$Group)
rownames(colData)=sample_info$SampleID

if (length(colData$Group)==length(unique(colData$Group))){
  print('No replicate analysis')
  dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design= ~1)
} else {
  dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design= ~Group)
}
#q()
#dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design= ~Group)
print(dds)

print('Prefiltering: removing features with zero counts')
# Prefiltering
dds <- dds[rowSums(counts(dds)) >= 1, ]
dds <- DESeq(dds)
print(dds)

# Exploring data  #
#-----------------#
# Transforming values
print('Transforming values to rlog')
rld <- rlog(dds, blind=FALSE)
#rld_blind <- rlog(dds, blind=TRUE)

# Heatmap genes
print('Plotting heatmap')
#select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
pdf(paste0(outdir,'/Heatmap_allSamples', ".pdf"), width=9,height=9)
heatmap.2(cor(assay(rld)), margins=c(12,12), scale=c('none'), density.info='density', trace='none')
dev.off()

#y=as.matrix(dist(t(assay(rld))))
#pdf(paste0(outdir,'/Heatmap_allSamples_dist', ".pdf"), width=9,height=9)
#heatmap.2(y, margins=c(12,12), scale=c('none'), density.info='density', trace='none')
#dev.off()

# PCA
print('Plotting PCA')
pca_data=plotPCA(rld, intgroup=c('Group'), returnData=TRUE)
percentVar=round(100*attr(pca_data,'percentVar'))
pdf(paste0(outdir,'/PCA_allSamples',".pdf"), width=9,height=9)
ggplot(pca_data, aes(PC1,PC2, color=Group, label=rownames(pca_data))) + geom_point() + geom_text(show_guide=F) + 
  xlab(paste0('PC1: ', percentVar[1], '% variance'))  + ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
  theme(panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.3,linetype=2), panel.grid.minor=element_line(colour='grey',size=.3,linetype=2)) + xlim(min(pca_data$PC1)-2,max(pca_data$PC1)+2)
dev.off()                                                      
#pdf(paste0(outdir,'/PCA_allSamples', ".pdf"), width=9,height=9)
#plotPCA(rld, intgroup=c('Group'))
#dev.off()

#Libsize
#pdf(paste0(outdir,'/Libsize', ".pdf"), width=9)
#par(mar=c(12.8,5.1,4.1,2.1))
#barplot(colSums(counts(dds, normalize=FALSE)), names=colnames(counts(dds, normalize=FALSE)), las=2, ylab="") # library sizes
#dev.off()

# Doing Differential Gene Expression  #
#-------------------------------------#
#comparisons=read.csv(comparisons, header=TRUE)
#if (design=='pairedSamples'){
#  pairedDesign=TRUE
#} else if(design=='non-pairedSamples'){
#  pairedDesign=FALSE
#}

#multipleComparison(data,comparisons,pairedDesign, gtf.file)
