#!/usr/local/bin/Rscript

suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))

indir = commandArgs(TRUE)[1]
outdir = commandArgs(TRUE)[2]
gtf.file= commandArgs(TRUE)[3]

dir.create(outdir, showWarnings = FALSE)

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

# Exclude reads not mapping to features, ambiguous, too low qual, multiple locations
data.notcounted=data[which(rownames(data)=='__no_feature'):dim(data)[1],]

# Extract reads mapping to features
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),]
total.data.counted=colSums(data)

# Read GTF File
GTF <- import.gff(gtf.file, format="gtf", feature.type="gene")
df=as.data.frame(GTF)
df=subset(df,select=c(seqnames,start,end,width,strand,gene_id,gene_name,gene_biotype))
rownames(df)=df$gene_id
df=df[rownames(data),]
df$gene_name[is.na(df$gene_name)]=df$gene_id[is.na(df$gene_name)]

# Set the shapes *this could vary depending on the number of samples*
shape_values=0:64
shape_values=shape_values[-c(27:33)]

# Library proportion of different biotypes #
#------------------------------------------
type=unique(df$gene_biotype)
for (i in 1:length(type)){ 
	name=type[i]
	assign(name,colSums(data[which(df$gene_biotype==type[i]),])/colSums(data[,]))
}

matrix=matrix(ncol=ncol(data),nrow=length(type))
rownames(matrix)=type
colnames(matrix)=colnames(data)
for (i in 1:length(type)){
	matrix[type[i],]=eval(as.name(type[i]))
}
matrix=matrix[-which(rowSums(matrix)==0),]
matrix.melt=melt(matrix)
colnames(matrix.melt)[2]='Sample'
png(paste0(outdir,'/biotype_proportion.png'), width=1400, height=700,res=150)
ggplot(matrix.melt, aes(x=X1, y=value, colour=Sample, shape=Sample)) + geom_point(size=1) + 
theme(axis.text.x=element_text(), axis.text=element_text(size=10, colour='black'), axis.title.x=element_blank(), legend.text=element_text(size=8), 
   panel.grid.major=element_line(colour='bisque3', linetype='dashed', size=.3), panel.background=element_rect(colour='white')) + 
scale_shape_manual(values=shape_values[1:ncol(data)]) + ylab('Library Proportion') + xlab('') +
coord_flip() 
dev.off()

# Library proportion of different genes #
#---------------------------------------#
prop=t(t(data)/total.data.counted)
prop.adjusted=prop/df$width

q=min(apply(prop.adjusted,2,function(x){quantile(x, .999)}))
genes_most_expressed=sapply(1:ncol(prop.adjusted),function(x){which(prop.adjusted[,x]>.0000002)})
genes_most_expressed=sapply(1:ncol(prop.adjusted),function(x){which(prop.adjusted[,x]>q)})

genes=unique(names(unlist(genes_most_expressed)))
prop.top=prop[genes,]
prop.top=cbind(prop.top,subset(df[genes,],select=c(seqnames,strand,gene_name,gene_biotype)))
prop.top=melt(prop.top)
prop.top$variable=factor(prop.top$variable)

# Top expressed genes all
png(paste0(outdir,'/top_expressed_genes.png'),width=1300, height=800, res=150)
ggplot(prop.top,aes(x=gene_name,y=value,group=variable,color=variable, shape=variable)) +
geom_point(size=1.5) + facet_wrap(~gene_biotype, scales='free_x')+ xlab('') + ylab('Library proportion') + theme(axis.text=element_text(size=6),axis.text.x=element_text(angle=90), axis.title=element_text(size=10), legend.text=element_text(size=8), legend.title=element_text(size=8)) + scale_shape_manual(values=shape_values[1:ncol(data)])
dev.off()

# Top expressed gene protein coding 
prop.top=prop.top[prop.top$gene_biotype=='protein_coding',]
prop.top$gene_name=factor(prop.top$gene_name,levels=sort(unique(prop.top$gene_name), decreasing=TRUE))
colnames(prop.top)[5]='Sample'

png(paste0(outdir,'/top_expressed_genes_protein_coding.png'),width=1100, height=1000, res=150)
ggplot(prop.top,aes(x=gene_name,y=value,group=Sample,color=Sample, shape=Sample)) +
geom_point(size=2.5) +
facet_wrap(~gene_biotype, scales='free_x')+ xlab('') + ylab('Library proportion') +
theme(axis.text=element_text(size=8, color='black'),axis.text.x=element_text(),
      axis.title=element_text(size=10),
      legend.text=element_text(size=8),  legend.title=element_text(size=8),
      panel.grid.major=element_line(colour='bisque3', linetype='dashed', size=.3),
      panel.background=element_rect(colour='white')) +
scale_shape_manual(values=shape_values[1:ncol(data)]) +
coord_flip()
dev.off()



# Do the same but Ignore ribosomal and mitochondrial protein genes
prop=t(t(data)/total.data.counted)
prop.adjusted=prop/df$width

rps=df[grep('Rps',df$gene_name, ignore.case=TRUE),]$gene_id
rpl=df[grep('Rpl',df$gene_name, ignore.case=TRUE),]$gene_id
mt=df[grep('mt-',df$gene_name, ignore.case=TRUE),]$gene_id
remove.genes=c(rps,rpl,mt)

prop.adjusted.noribosomalprotein=prop.adjusted[!rownames(prop.adjusted)%in%remove.genes,]

q=min(apply(prop.adjusted.noribosomalprotein,2,function(x){quantile(x, .999)}))
genes_most_expressed=sapply(1:ncol(prop.adjusted),function(x){which(prop.adjusted[,x]>.0000002)})
genes_most_expressed=sapply(1:ncol(prop.adjusted.noribosomalprotein),function(x){which(prop.adjusted.noribosomalprotein[,x]>q)})

genes=unique(names(unlist(genes_most_expressed)))
prop.top=prop[genes,]
prop.top=cbind(prop.top,subset(df[genes,],select=c(seqnames,strand,gene_name,gene_biotype)))
prop.top=melt(prop.top)
prop.top$variable=factor(prop.top$variable)


# Top expressed protein coding genes
prop.top=prop.top[prop.top$gene_biotype=='protein_coding',]
prop.top$gene_name=factor(prop.top$gene_name,levels=sort(unique(prop.top$gene_name), decreasing=TRUE))
colnames(prop.top)[5]='Sample'
shape_values=0:64
shape_values=shape_values[-c(27:33)]
#print(shape_values)

png(paste0(outdir,'/top_expressed_genes_protein_coding_ignore_ribosomal_and_mt.png'),width=1100, height=1100, res=150)
ggplot(prop.top,aes(x=gene_name,y=value,group=Sample,color=Sample, shape=Sample)) +
geom_point(size=2.5) + 
facet_wrap(~gene_biotype, scales='free_x')+ xlab('') + ylab('Library proportion') + 
theme(axis.text=element_text(size=8, color='black'),axis.text.x=element_text(), 
      axis.title=element_text(size=10), 
      legend.text=element_text(size=8),  legend.title=element_text(size=8), 
      panel.grid.major=element_line(colour='bisque3', linetype='dashed', size=.3), 
      panel.background=element_rect(colour='white')) + 
scale_shape_manual(values=shape_values[1:ncol(data)]) + 
coord_flip()
dev.off()
