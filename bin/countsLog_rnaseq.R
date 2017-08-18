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

indir = commandArgs(TRUE)[1]
outdir = commandArgs(TRUE)[2]
mapping_sum = commandArgs(TRUE)[3]

###############################
# Reading counted reads files #
###############################

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

# Exclude reads not mapping to features, ambiguous, too low qual, multiple locations
data.notcounted=data[which(rownames(data)=='__no_feature'):dim(data)[1],]

# Keep read mapping to features
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),]
total.data.counted=colSums(data)

# ***************************************************************************************#
# Distribution of counts: From the uniquely mapped reads, we get how many of these were in features, not in feature, ambiguous, and too low quality
n_reads_in_genes=colSums(data)
n_reads_no_feature=t(data.notcounted["__no_feature",])[,1]
n_reads_ambiguous=t(data.notcounted["__ambiguous",])[,1]
n_reads_lowqual=t(data.notcounted["__too_low_aQual",])[,1]

# Using the mapping summary, we get the uniquely mapped reads
n_reads_total=read.csv(mapping_sum,row.names=1)
sample_names=rownames(n_reads_total)
n_reads_total=n_reads_total$Mapped_num
names(n_reads_total)=sample_names

# Bind matrix
counting_summ=cbind(n_reads_total, n_reads_in_genes[sample_names], n_reads_no_feature[sample_names], n_reads_ambiguous[sample_names], n_reads_lowqual[sample_names])
colnames(counting_summ)=c('Uniquely_mapped','In_features', 'Not_in_features', 'Ambiguous', 'Low_aQual')
write.csv(counting_summ, paste0(outdir,'/counts_in_features.csv'))

# Plot
table=t(t(counting_summ[,-1]/counting_summ[,1]*100))
data.melt=melt(table)
png(paste0(outdir,"/counts_in_features.png"), width=1000, height=800, res=120)
ggplot(data.melt, aes(x=X1, y=value, fill=X2, group=X2 )) + geom_bar(stat='identity') + 
       theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) + labs(fill="") + xlab('Sample') + ylab('Percentage')
dev.off()

# Numbers of genes detected
n_genes=sapply(data,function(x){table(x==0)['FALSE']})
names(n_genes)=gsub('.FALSE','',names(n_genes))
counting_summ=cbind(N_features=n_genes)
write.csv(counting_summ, paste0(outdir,'/Number_of_features_detected.csv'))


# ************************************************ #
# Distribution of counts 
data.melt=suppressMessages(melt(data))
pdf(paste0(outdir,"/countsDistributionHist.pdf"), width=12, height=7)
colnames(data.melt)[1]='Sample'
ggplot(data.melt, aes(x=log(value), fill=Sample )) + geom_bar(binwidth=.2) + 
	theme(legend.text=element_text(size=10), legend.key.size = unit(.78, "cm"))
suppressMessages(dev.off())

png(paste0(outdir,"/countsDistributionHist_per_sample.png"), width=1400, height=800, res=150)
colnames(data.melt)[1]='Sample'
ggplot(data.melt, aes(x=log(value), fill=Sample )) + geom_bar(binwidth=.2) + facet_wrap(~Sample) + theme(legend.position='none', strip.text=element_text(size=9)) + xlab('Log2 raw counts')
suppressMessages(dev.off())

# Distribution lines
shape_values=0:64
shape_values=shape_values[-c(27:33)]

distribution=as.data.frame(matrix(0,ncol=length(files),nrow=14))
colnames(distribution)=colnames(data)
rownames(distribution)=c('0', '1-100', '101-200','200-500','501-1,000', '1,001-5,000', '5,001-10,000', '10,001-20,000', '20,001-30,000', '30,000-40,000', '40,001-50,000', '50,001-100,000', '100,001-50,000', '>500,000')
for (i in 1:length(colnames(data))) {
  sample_name=colnames(data)[i]
  x=data[,sample_name]
  zero=table(x==0)['TRUE']
  bin1=table(x>0 & x<=100)['TRUE']
  bin2=table(x>100 & x<=200)['TRUE']
  bin3=table(x>200 & x<=500)['TRUE']
  bin4=table(x>500 & x<=1000)['TRUE']
  bin5=table(x>1000 & x<=5000)['TRUE']
  bin6=table(x>5000 & x<=10000)['TRUE']
  bin7=table(x>10000 & x<=20000 )['TRUE']
  bin8=table(x>20000 & x<=30000 )['TRUE']
  bin9=table(x>30000 & x<=40000 )['TRUE']
  bin10=table(x>40000 & x<=50000 )['TRUE']
  bin11=table(x>50000 & x<=100000 )['TRUE']
  bin12=table(x>100000 & x<=500000 )['TRUE']
  bin13=table(x>500000)['TRUE']
  sample=c(zero,bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13)
  sample[is.na(sample)]=0
  distribution[,sample_name]=sample
}
distribution=distribution/colSums(distribution)*100
distribution.melt=melt(cbind(bins=factor(row.names(distribution), levels=row.names(distribution)), distribution))

png(paste0(outdir,'/countDistribution_lines.png'), width=1200, height=800, res=180)
ggplot(distribution.melt,aes(x=bins,y=value,color=variable,group=variable, shape=variable)) + geom_path() + geom_point(size=.8) +  
  xlab('Raw count value') + ylab('Percentage of genes') + 
  scale_y_continuous(breaks=seq(0,70,10), labels=c('0%','10%','20%','30%','40%','50%','60%','70%')) + 
  scale_colour_discrete(name='Sample') + 
  scale_shape_manual(values=shape_values[1:ncol(data)]) + 
  theme(axis.text.x=element_text(angle=90), axis.text.y=element_text(),legend.text=element_text(size=9), legend.key.size = unit(.48, "cm"))
dev.off()

# Distribution lines removing all zeros
data_filtered=data[!rowSums(data)==0,]
distribution=as.data.frame(matrix(0,ncol=length(files),nrow=14))
colnames(distribution)=colnames(data_filtered)
rownames(distribution)=c('0', '1-100', '101-200','200-500','501-1,000', '1,001-5,000', '5,001-10,000', '10,001-20,000', '20,001-30,000', '30,000-40,000', '40,001-50,000', '50,001-100,000', '100,001-50,000', '>500,000')
for (i in 1:length(colnames(data_filtered))) {
  sample_name=colnames(data_filtered)[i]
  x=data_filtered[,sample_name]
  zero=table(x==0)['TRUE']
  bin1=table(x>0 & x<=100)['TRUE']
  bin2=table(x>100 & x<=200)['TRUE']
  bin3=table(x>200 & x<=500)['TRUE']
  bin4=table(x>500 & x<=1000)['TRUE']
  bin5=table(x>1000 & x<=5000)['TRUE']
  bin6=table(x>5000 & x<=10000)['TRUE']
  bin7=table(x>10000 & x<=20000 )['TRUE']
  bin8=table(x>20000 & x<=30000 )['TRUE']
  bin9=table(x>30000 & x<=40000 )['TRUE']
  bin10=table(x>40000 & x<=50000 )['TRUE']
  bin11=table(x>50000 & x<=100000 )['TRUE']
  bin12=table(x>100000 & x<=500000 )['TRUE']
  bin13=table(x>500000)['TRUE']
  sample=c(zero,bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13)
  sample[is.na(sample)]=0
  distribution[,sample_name]=sample
}
distribution=distribution/colSums(distribution)*100
distribution.melt=melt(cbind(bins=factor(row.names(distribution), levels=row.names(distribution)), distribution))

png(paste0(outdir,'/countDistribution_lines_nozeros.png'),  width=1200, height=800, res=180)
ggplot(distribution.melt,aes(x=bins,y=value,color=variable,group=variable, shape=variable)) + geom_path() + 
  xlab('Raw count value (excluding undetected (i. e. 0 count) features in all samples )') + ylab('Percentage of genes') + 
  scale_y_continuous(breaks=seq(0,70,10), labels=c('0%','10%','20%','30%','40%','50%','60%','70%')) + 
  scale_colour_discrete(name='Sample') + 
  scale_shape_manual(values=shape_values[1:ncol(data)]) +
  theme(axis.text.x=element_text(angle=90), axis.text.y=element_text(),legend.text=element_text(size=9), legend.key.size = unit(.48, "cm"))
dev.off()


