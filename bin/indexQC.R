#!/usr/local/bin/Rscript
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))

indir = commandArgs(TRUE)[1]
outdir= commandArgs(TRUE)[2]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
files=list.files(indir, recursive=TRUE, pattern='fastqc_data.txt$')

# Get table with number of reads #
#--------------------------------#
info_matrix=matrix(ncol=1,nrow=length(files))
colnames(info_matrix)='NReads'
for (i in 1:length(files)){
	info=readLines(paste0(indir,'/',files[i]), n=7)[7]
	info=as.numeric(gsub('.*\\t','',info))
	info_matrix[i,1]=info
}
rownames(info_matrix)=gsub('.*\\/','',gsub('/fastqc_data.txt','',files))
rownames(info_matrix)=gsub('R1.*fastqc','R1',rownames(info_matrix))
rownames(info_matrix)=gsub('R2.*fastqc','R2',rownames(info_matrix))

# Make bar plot #
#---------------#
total=sum(info_matrix[,1])
info_matrix.melt=melt(info_matrix)
png(paste0(outdir,'/nreads.png'),width=900, height=700, res=100)
ggplot(info_matrix.melt, aes(x=X1, y=((value/total)*100))) + geom_bar(stat='identity') + theme(axis.text.x=element_text(face='bold',size=11,angle=90,hjust=0.5,vjust=0.5, colour='black'), axis.title.x=element_blank(), axis.text.y=element_text(size=11, colour='black'), axis.title.y=element_text(size=14, colour='black', vjust=1.5)) + ylab('Percentage')
dev.off()

info_matrix[,1]=format(info_matrix[,1], big.mark=',')
write.csv(info_matrix,paste0(outdir,"/nreads.csv"))
