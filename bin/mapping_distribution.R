#!/usr/local/bin/Rscript

dir = commandArgs(TRUE)[1]
name = commandArgs(TRUE)[2]

x=scan(paste0(dir,"/",name,'_metrics.txt'), what='character',skip=6,nlines=2, quiet=TRUE)
x=x[-grep("SAMPLE",x)]
x=x[-grep("LIBRARY",x)]
x=x[-grep("READ_GROUP",x)]
x=data.frame(cbind(matrix(x,ncol=2),NA),stringsAsFactors=FALSE)
x[,3]=format(round(as.numeric(x[,2]),3), big.mark=',', scientific=F)

# percentage
pct=x[grep("PCT_",x[,1]),]
pct[,2]=as.numeric(pct[,2])*100
write.table(pct[,1:2], paste0(dir,'/pct_dist_',name, '.txt'), quote=F, row.names=F, col.names=F)

