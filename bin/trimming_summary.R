#!/usr/local/bin/Rscript
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))

indir_raw = commandArgs(TRUE)[1]
indir_trimmed = commandArgs(TRUE)[2]
outdir= commandArgs(TRUE)[3]


# Get table with number of reads raw #
#------------------------------------#
files=list.files(indir_raw, recursive=TRUE, pattern='fastqc_data.txt$')
info_matrix=matrix(ncol=2,nrow=length(files))
colnames(info_matrix)=c('Raw_reads','After_trimming')
for (i in 1:length(files)){
	info=readLines(paste0(indir_raw,'/',files[i]), n=7)[7]
	info=as.numeric(gsub('.*\\t','',info))
	info_matrix[i,1]=info
}
rownames(info_matrix)=gsub('.*\\/','',gsub('/fastqc_data.txt','',files))
rownames(info_matrix)=gsub('R1.*fastqc','R1',rownames(info_matrix))
rownames(info_matrix)=gsub('R2.*fastqc','R2',rownames(info_matrix))

# Get table with number of reads trimmed #
#----------------------------------------#
files=list.files(indir_trimmed, recursive=TRUE, pattern='fastqc_data.txt$')
info_matrix_trimmed=matrix(ncol=1,nrow=length(files))
colnames(info_matrix_trimmed)=c('After_trimming')
for (i in 1:length(files)){
	info=readLines(paste0(indir_trimmed,'/',files[i]), n=7)[7]
	info=as.numeric(gsub('.*\\t','',info))
	info_matrix_trimmed[i,1]=info
}
rownames(info_matrix_trimmed)=gsub('.*\\/','',gsub('/fastqc_data.txt','',files))
rownames(info_matrix_trimmed)=gsub('R1.*fastqc','R1',rownames(info_matrix_trimmed))
rownames(info_matrix_trimmed)=gsub('R2.*fastqc','R2',rownames(info_matrix_trimmed))


# Get table with number of reads trimmed and pct #
#------------------------------------------------#
info_matrix[,2]=info_matrix_trimmed[rownames(info_matrix),]
info_matrix=cbind(info_matrix,pct_removed=(1-info_matrix[,2]/info_matrix[,1])*100)

info_matrix[,3]=round(info_matrix[,3], digits=4)
info_matrix[,1]=format(info_matrix[,1], big.mark=',')
info_matrix[,2]=format(as.numeric(info_matrix[,2]), big.mark=',')

write.csv(info_matrix,paste0(outdir,"/nreads_trimming.csv"))
