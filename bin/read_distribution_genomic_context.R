#!/usr/local/bin/Rscript

suppressMessages(library(RColorBrewer))
suppressMessages(library(matrixStats))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
suppressMessages(require(grid))

in_dir=commandArgs(TRUE)[1]
out_dir= commandArgs(TRUE)[2]
patt='pct_dist_'

files=list.files(in_dir, pattern=patt)
sample_names=gsub('.txt','',gsub(patt,'',files))

# Read files into matrix m
m=data.frame(matrix(ncol=length(files),nrow=8),stringsAsFactors = FALSE)
m=read.table(paste0(in_dir,'/',files[1]),row.names=1, stringsAsFactors = FALSE)
colnames(m)=files[1]

for (i in 2:length(files)){
  temp=read.table(paste0(in_dir,'/',files[i]),row.names=1, stringsAsFactors = FALSE)
  m=cbind(m,temp[rownames(m),])
  colnames(m)[i]=files[i]
}

colnames(m)=gsub('.txt','',gsub(patt,'',colnames(m)))
rownames(m)=gsub('PCT_','',rownames(m))
rownames(m)=gsub('_BASES','',rownames(m))


# Plot including ribosomal, coding, utr, intronic and intergenic pcts, save strand in different vector
m=melt(cbind(cat=rownames(m),m))
strand=m[grep('STRAND',m[,1]),]
m=m[-grep('STRAND',m[,1]),]
m=m[-grep('USABLE',m[,1]),]
m=m[-grep('MRNA',m[,1]),]
png(paste0(out_dir,"/read_distribution_genomic_context.png"), height=700, width=900, res=150)
ggplot(m, aes(x=factor(variable),y=value, fill=cat)) + geom_bar(width=0.5, stat="identity") + coord_flip() + theme(axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'), axis.title.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank())  
dev.off()

# Table of strand accuracy
x=cbind(strand[,3])
rownames(x)=strand[,2]
colnames(x)=c('Percentage')
write.csv(x,paste0(out_dir,'/strand_mappingQC.csv'))
