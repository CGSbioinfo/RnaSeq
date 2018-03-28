#!/usr/local/bin/Rscript

print('Loading libraries')
suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
#suppressMessages(library(GGally))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(matrixStats))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
source("bin/deseq2_functions.R")


#############################################################################################
#Commented out the section about 'clones' and added a space between the 'plot_device' and '='
#############################################################################################

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
plot_device = read_args_line(grep('^plot_device =', arguments_file$V1, value=TRUE))
print(plot_device)

print('Creating output directory')
dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

# Read comparisons file
#comparisonsFile=commandArgs(TRUE)[2]
comparisons=read.csv(comparisons, header=TRUE)

# Reading counted reads files #
#-----------------------------#
print('Reading counted reads files')
files=list.files(indir, pattern = "_count.txt", full.names=TRUE)

data=read.table(files[1], row.names=1)
colnames(data)=files[1]

for (j in 2:length(files)){
  name=files[j]
  assign(name, read.table(name))
  sample=eval(as.name(name))
  index=match(rownames(data), sample[,1])
  data=cbind(data, sample[index,2])
  colnames(data)[j]=name
}
colnames(data)=gsub('_count.txt','',colnames(data))
colnames(data)=gsub('.*/','',colnames(data))

print('Keep reads mapping to features:')
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),]


# Getting information about the groups #
#--------------------------------------#
print('Reading group info')
sample_info = read.csv(sample_info)
sample_info=sample_info[order(sample_info$Group),]
group <- c(as.character(sample_info$Group))

# Matching the data matrix to the sample_info table #
#---------------------------------------------------#
data=data[,match(as.character(sample_info$SampleID),colnames(data))]

# remove sample?
removeSample=commandArgs(TRUE)[3]
print(removeSample)
if (!is.na(removeSample)){
 sample_info=sample_info[-which(sample_info$SampleID%in%removeSample),]
 data=data[,-which(colnames(data)%in%removeSample)]
}
#sample_info=sample_info[-which(sample_info$SampleID%in%removeSample),]

# Create a DDS object
print('Creating a DDS object')
colData<-sample_info
rownames(colData)=sample_info$SampleID

dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design= ~Group)

# Prefiltering
dds <- dds[rowSums(counts(dds)) >= 1, ]
dds <- DESeq(dds)

# Transform values to rlog
print('Transforming values to rlog')
rld <- rlog(dds, blind=FALSE)

png(paste0(outdir,"/Heatmap_allSamples.png"), height=850, width=900, res=120)
heatmap.2(cor(assay(rld)), margins=c(12,12), scale=c('none'), density.info='density', trace='none')
dev.off()

print('Plotting PCA')
pca_data=plotPCA_4pcs(rld, intgroup=c('Group'), returnData=TRUE)
percentVar=round(100*attr(pca_data,'percentVar'))
names(percentVar)=colnames(pca_data)[1:4]
pca_data$Time=gsub(".*-","",pca_data$Group)
pca_data$KO=gsub("-.*","",pca_data$Group)

# KO
p1=pca_plot_single(pca_data, pcs=c('PC1','PC2'), col='KO', percentVar)
p2=pca_plot_single(pca_data, pcs=c('PC1','PC3'), col='KO', percentVar)
p3=pca_plot_single(pca_data, pcs=c('PC1','PC4'), col='KO', percentVar)
p4=pca_plot_single(pca_data, pcs=c('PC2','PC3'), col='KO', percentVar)
p5=pca_plot_single(pca_data, pcs=c('PC2','PC4'), col='KO', percentVar)
p6=pca_plot_single(pca_data, pcs=c('PC3','PC4'), col='KO', percentVar)

#png(paste0(outdir,"/PCA_allSamples_KO.png"), height=880, width=1680, res=80 )
grid<-grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
ggsave(filename=paste0(outdir,'/PCA_allSamples_KO.', plot_device), width=40, height=20,unit="cm", plot=grid)
#dev.off()


# Time 
p1=pca_plot_single(pca_data, pcs=c('PC1','PC2'), col='Time', percentVar)
p2=pca_plot_single(pca_data, pcs=c('PC1','PC3'), col='Time', percentVar)
p3=pca_plot_single(pca_data, pcs=c('PC1','PC4'), col='Time', percentVar)
p4=pca_plot_single(pca_data, pcs=c('PC2','PC3'), col='Time', percentVar)
p5=pca_plot_single(pca_data, pcs=c('PC2','PC4'), col='Time', percentVar)
p6=pca_plot_single(pca_data, pcs=c('PC3','PC4'), col='Time', percentVar)

#png(paste0(outdir,"/PCA_allSamples_time.png"), height=880, width=1680, res=80 )
grid<-grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
ggsave(filename=paste0(outdir,'/PCA_allSamples_time.', plot_device), width=40, height=20,unit="cm", plot=grid)
#dev.off()


# Number of reads
#pca_data=left_join(pca_data, data.frame(name=names(total.data.counted), total.data.counted), by='name')
#p1=pca_plot_single(pca_data, pcs=c('PC1','PC2'), col='total.data.counted', percentVar)
# p2=pca_plot_single(pca_data, pcs=c('PC1','PC3'), col='total.data.counted', percentVar)
# p3=pca_plot_single(pca_data, pcs=c('PC1','PC4'), col='total.data.counted', percentVar)
# p4=pca_plot_single(pca_data, pcs=c('PC2','PC3'), col='total.data.counted', percentVar)
# p5=pca_plot_single(pca_data, pcs=c('PC2','PC4'), col='total.data.counted', percentVar)
# p6=pca_plot_single(pca_data, pcs=c('PC3','PC4'), col='total.data.counted', percentVar)
# 
# #png(paste0(outdir,"/PCA_allSamples_nreads.png"), height=880, width=1680, res=80 )
# grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
# ggsave(filename=paste0(outdir,'/PCA_allSamples_nreads.', plot_device), width=30, height=20,unit="cm", plot=grid)



# Clone identity
#pca_data=left_join(pca_data, data.frame(name=sample_info$SampleID, clone=sample_info$Clone), by='name')
#p1=pca_plot_single(pca_data, pcs=c('PC1','PC2'), col='clone', percentVar)
#p2=pca_plot_single(pca_data, pcs=c('PC1','PC3'), col='clone', percentVar)
#p3=pca_plot_single(pca_data, pcs=c('PC1','PC4'), col='clone', percentVar)
#p4=pca_plot_single(pca_data, pcs=c('PC2','PC3'), col='clone', percentVar)
#p5=pca_plot_single(pca_data, pcs=c('PC2','PC4'), col='clone', percentVar)
#p6=pca_plot_single(pca_data, pcs=c('PC3','PC4'), col='clone', percentVar)

#png(paste0(outdir,"/PCA_allSamples_clone.png"), height=880, width=1680, res=80 )
grid<-grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
ggsave(filename=paste0(outdir,'/PCA_allSamples_clones.', plot_device), width=40, height=20,unit="cm", plot=grid)




# Comparisons 

for (i in 1:nrow(comparisons)){
  groups=as.character(unlist(comparisons[i,]))      
  print(groups)
  samples=as.character(sample_info$SampleID[sample_info$Group %in% groups])
  
  colData.sub <-sample_info[sample_info$SampleID%in%samples,]
  colData.sub=droplevels(colData.sub)
  rownames(colData.sub)=colData.sub$SampleID

  data.sub=data[,colnames(data)%in%samples]
  data.sub=data.sub[,match(as.character(colData.sub$SampleID),colnames(data.sub))]

  # Create DDS
  dds <- DESeqDataSetFromMatrix(countData= data.sub, colData=colData.sub, design= ~Group)

  # Prefiltering
  dds <- dds[rowSums(counts(dds)) >= 1, ]
  dds <- DESeq(dds)

  #dds.sub=dds[,samples]
  #colData(dds.sub)$Group=droplevels(colData(dds.sub)$Group)
  #print(colData(dds.sub))
  #print('Prefiltering: removing features with zero counts')
  # Prefiltering
  #dds.sub <- dds.sub[rowSums(counts(dds.sub)) >= 1, ]
  #dds.sub <- DESeq(dds.sub)
  
  print('Transforming values to rlog')
  rld <- rlog(dds, blind=FALSE)
 
  png(paste0(outdir,'/Heatmap_', paste(groups, collapse='_vs_'),".png"), height=850, width=900, res=120)
  heatmap.2(cor(assay(rld)), margins=c(14,14), scale=c('none'), density.info='density', trace='none')
  dev.off()
  
  print('Plotting PCA')
  pca_data=plotPCA_4pcs(rld, intgroup=c('Group'), returnData=TRUE)
  percentVar=round(100*attr(pca_data,'percentVar'))
  names(percentVar)=colnames(pca_data)[1:4]
 
  colnames(pca_data)[6]='SampleID'
  pca_data=left_join(pca_data, colData(dds)%>%as.data.frame, by='SampleID')
  rownames(pca_data)=pca_data$SampleID
  colnames(pca_data)[5]='Group'
  
  p1=pca_plot_single_label(pca_data, pcs=c('PC1','PC2'), col='Group', shape='Group', percentVar)
  p2=pca_plot_single_label(pca_data, pcs=c('PC1','PC3'), col='Group', shape='Group', percentVar)  
  p3=pca_plot_single_label(pca_data, pcs=c('PC1','PC4'), col='Group', shape='Group', percentVar)
  p4=pca_plot_single_label(pca_data, pcs=c('PC2','PC3'), col='Group', shape='Group', percentVar)
  p5=pca_plot_single_label(pca_data, pcs=c('PC2','PC4'), col='Group', shape='Group', percentVar)
  p6=pca_plot_single_label(pca_data, pcs=c('PC3','PC4'), col='Group', shape='Group', percentVar)
  
  #png(paste0(outdir,'/PCA_',paste(groups, collapse='_vs_'),".png"), height=880, width=1680, res=80 )
  grid<-grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
  ggsave(filename=paste0(outdir,'/PCA_',paste(groups, collapse='_vs_'),".", plot_device), width=40, height=20,unit="cm", plot=grid)
  #dev.off()
}

