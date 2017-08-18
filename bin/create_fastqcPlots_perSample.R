#!/usr/local/bin/Rscript

suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape)))
suppressWarnings(suppressMessages(require(grid)))
suppressWarnings(suppressMessages(library(xtable)))

in_dir=commandArgs(TRUE)[1]
sample_name=commandArgs(TRUE)[2]
readType=commandArgs(TRUE)[3] 
outdir=commandArgs(TRUE)[4]
suffix=commandArgs(TRUE)[5]
plot_device=commandArgs(TRUE)[6]
if (is.na(suffix)){
  suffix=''
}

files=list.files(path = in_dir, full.names = TRUE, recursive = TRUE)
files=files[grep(sample_name,files)]

generate_qc_plot=function(files,sample_name,type){
  sample_temp=sample_name
  #files_temp=files[grep(sample_temp,files)]
  files_temp=files[grep(type,files)]
  r1=files_temp[grep('_R1_',files_temp)]
  if(file.info(r1)$size==0){
      print(paste0('File ', type, ', sample ', sample_name, ', R1 is empty'))
      dat_r1=data.frame()
  } else {
      dat_r1=read.table(r1, stringsAsFactors = FALSE)
  }
  if (readType=='pairedEnd'){
    r2=files_temp[grep('_R2_',files_temp)]
    if(file.info(r2)$size==0){
      print(paste0('File ', type, ', sample ', sample_name, ', R2 is empty'))
      dat_r2=data.frame()
    } else {
    dat_r2=read.table(r2, stringsAsFactors = FALSE)
  }
  }
  if (type=='per_base_sequence_quality.txt'){
    if(file.info(r1)$size==0){
        dat_r1=data.frame(matrix(ncol=7))
    }
    colnames(dat_r1)=c('Base','Mean','Median','Lower_quartile', 'Upper_quartile', 'x10th_percentile', 'x90th_percentile')
    dat_r1$Base <- factor(dat_r1$Base, as.character(dat_r1$Base))
    dat_r1$panel='Read 1' 
    if (readType=='pairedEnd'){
      if(file.info(r2)$size==0){
        dat_r2=data.frame(matrix(ncol=7))
      }
      colnames(dat_r2)=c('Base','Mean','Median','Lower_quartile', 'Upper_quartile', 'x10th_percentile', 'x90th_percentile')
      dat_r2$Base <- factor(dat_r2$Base, as.character(dat_r2$Base))
      dat_r2$panel='Read 2'
      dat=rbind(dat_r1,dat_r2)
      p=ggplot(dat,
             aes(x = Base, ymin = x10th_percentile, lower = Lower_quartile, middle = Median, upper = Upper_quartile, 
                 ymax = x90th_percentile)) + 
        geom_boxplot(stat='identity', fill='yellow') + facet_wrap(~panel) + ggtitle(sample_temp) + 
        theme(axis.text.x = element_text(size = 5, angle=90, hjust=1, vjust=1), axis.text.y=element_text(size=6),
              plot.title = element_text(lineheight=.8, face="bold", size=9), 
	      axis.title=element_text(size=10), strip.text=element_text(size=8)) + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=20, alpha=0.1, fill="red") + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=20, ymax=28, alpha=0.1, fill="yellow") + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=28, ymax=Inf, alpha=0.1, fill="green")  
        xlab("Position of base in read")
    } else {
      dat=dat_r1
      p=ggplot(dat, 
             aes(x = Base, ymin = x10th_percentile, lower = Lower_quartile, middle = Median, upper = Upper_quartile, 
                 ymax = x90th_percentile)) + 
        geom_boxplot(stat='identity', fill='yellow') + ggtitle(sample_temp) + 
        theme(axis.text.x = element_text(size = 8, angle=90, hjust=1, vjust=1), axis.text.y=element_text(size=6), 
              plot.title = element_text(lineheight=.8, face="bold"), axis.title=element_text(size=10)) + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=20, alpha=0.1, fill="red") + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=20, ymax=28, alpha=0.1, fill="yellow") + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=28, ymax=Inf, alpha=0.1, fill="green") + 
        xlab("Position of base in read")
    }
    ggsave(filename=paste0(outdir,'/',sample_temp, suffix, '_per_base_sequence_qual.', plot_device), plot=p, height=3.5, width=10)
  } else if (type=='per_base_sequence_content.txt'){
    if(file.info(r1)$size==0){
        dat_r1=data.frame(matrix(ncol=5))
    }
    colnames(dat_r1)=c('Base','%G','%A','%T','%C')
    dat_r1$Base <- factor(dat_r1$Base, as.character(dat_r1$Base))
    dat_r1=melt(dat_r1)
    dat_r1$panel='Read 1'
    if (readType=='pairedEnd'){
      if(file.info(r2)$size==0){
        dat_r2=data.frame(matrix(ncol=5))
      }
      colnames(dat_r2)=c('Base','%G','%A','%T','%C')
      dat_r2$Base <- factor(dat_r2$Base, as.character(dat_r2$Base))
      dat_r2=melt(dat_r2)
      dat_r2$panel='Read 2'
      dat=rbind(dat_r1,dat_r2)
      p=ggplot(dat, aes(x=Base, y=value, group=variable, colour=variable)) + 
        geom_line() + facet_wrap(~panel) + ggtitle(sample_temp) + 
        theme(axis.text.x = element_text(size = 5, angle=90), axis.title.x=element_text(size=10),axis.title.y=element_blank(), axis.text.y=element_text(size=6),strip.text=element_text(size=8), legend.position="top", plot.title = element_text(size=9,lineheight=.8, face="bold", vjust=-1.5))  + xlab("Position in read") +  ylim(0, 100)
    } else {
      dat=dat_r1
      p=ggplot(dat, aes(x=Base, y=value, group=variable, colour=variable)) + geom_line() + 
        ggtitle(sample_temp) + theme(axis.text.x = element_text(size = 8, angle=90), axis.title.y=element_blank(), 
                                     legend.position="top", plot.title = element_text(lineheight=.5, face="bold", vjust=-1.5)) + 
        xlab("Position in read") +  ylim(0, 100) 
    }
    ggsave(filename=paste0(outdir,'/',sample_temp, suffix, '_per_base_sequence_content.', plot_device), plot=p, height=3.5, width=10)
  } else if (type=='kmer_content.txt'){
    if(file.info(r1)$size==0){
        dat_r1=data.frame(matrix(ncol=5))
    }
    colnames(dat_r1)=c('Sequence','Content','PValue', 'Obs_Exp_Max', 'Max_Obs_Exp_Position' )
    nbases=max( dat_r1$Max_Obs_Exp_Position)
    m=c()
    for (i in 1:(min(length(dat_r1$Sequence),6))){
      m=rbind(m,cbind(dat_r1$Sequence[i],1:nbases,0))
    }
    for (i in 1:(min(length(dat_r1$Sequence),6))){
      l1 = m[,1]==dat_r1$Sequence[i] 
      l2 = m[,2]== dat_r1$Max_Obs_Exp_Position[i]
      rnum=intersect(which(l1==TRUE), which(l2==TRUE))
      m[rnum,3]=dat_r1$Obs_Exp_Max[i]
    }
    dat_r1=as.data.frame(m)
    dat_r1$panel='Read 1'
    if (readType=='pairedEnd'){
      if(file.info(r2)$size==0){
        dat_r2=data.frame(matrix(ncol=5))
      }
      colnames(dat_r2)=c('Sequence','Content','PValue', 'Obs_Exp_Max', 'Max_Obs_Exp_Position' )
      nbases=max(dat_r2$Max_Obs_Exp_Position,1)
      m=c()
      for (i in 1:(min(length(dat_r2$Sequence),6))){
        m=rbind(m,cbind(dat_r2$Sequence[i],1:nbases,0))
      }
      for (i in 1:(min(length(dat_r2$Sequence),6))){
        l1 = m[,1]==dat_r2$Sequence[i] 
        l2 = m[,2]==dat_r2$Max_Obs_Exp_Position[i]
        rnum=intersect(which(l1==TRUE), which(l2==TRUE))
        m[rnum,3]=dat_r2$Obs_Exp_Max[i]
      }
      dat_r2=as.data.frame(m)
      dat_r2$panel='Read 2'
      dat=rbind(dat_r1,dat_r2)
      colnames(dat)=c("Sequence","Position","value","panel")
      dat$Position=factor(dat$Position, 1:nbases,ordered = TRUE)
      dat$value=round(x = as.numeric(as.character(dat$value)),digits = 0)
      p=ggplot(dat, aes(x=Position, y=value, group=Sequence, colour=Sequence)) + geom_line() + 
        facet_wrap(~panel) + ggtitle(sample_temp) +  
        theme(axis.text.x = element_text(size = 6, angle=90), axis.title.y=element_blank(), 
              legend.position="top", plot.title = element_text(lineheight=.8, face="bold", vjust=-1.5))  + 
        xlab("Position in read")
    } else {
      dat = dat_r1
      colnames(dat)=c("Sequence","Position","value","panel")
      dat$Position=factor(dat$Position, 1:nbases,ordered = TRUE)
      dat$value=round(x = as.numeric(as.character(dat$value)),digits = 0)
      p=ggplot(dat, aes(x=Position, y=value, group=Sequence, colour=Sequence)) + geom_line() + 
        ggtitle(sample_temp) + theme(axis.text.x = element_text(size = 6, angle=90), 
                                     axis.title.y=element_blank(), legend.position="top", 
                                     plot.title = element_text(lineheight=.8, face="bold", vjust=-1.5))  + 
        xlab("Position in read")
    } 
    ggsave(filename=paste0(outdir,'/',sample_temp, suffix, '_kmer_content.', plot_device), plot=p, height=3.5, width=10)
  }
}


generate_qc_plot(files, sample_name, type='per_base_sequence_quality.txt')
#print(sample_name)
#generate_qc_plot(files, sample_name, type='per_base_sequence_content.txt')
generate_qc_plot(files, sample_name, type='kmer_content.txt')
#print(sample_name)

