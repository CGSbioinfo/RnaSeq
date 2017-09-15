#!/usr/local/bin/Rscript

suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape)))
suppressWarnings(suppressMessages(require(grid)))
suppressWarnings(suppressMessages(library(xtable)))
suppressWarnings(suppressMessages(library(dplyr)))

in_dir=commandArgs(TRUE)[1]
sample_names=commandArgs(TRUE)[2]
readType=commandArgs(TRUE)[3]
outdir=commandArgs(TRUE)[4]
suffix=commandArgs(TRUE)[5]
plot_device=commandArgs(TRUE)[6]
if (is.na(suffix)){
  suffix=''
}

files=list.files(path = in_dir, full.names = TRUE, recursive = TRUE)
sample_names=read.table(sample_names)[,1]

shape_values=0:64
shape_values=shape_values[-c(27:33,35,40,43,45,46,47)]
print(shape_values)

# Per base quality
#-----------------
qual_scores=files[grep('per_base_sequence_quality',files)]
mr1=data.frame(matrix(ncol=8))
colnames(mr1)=c('Base','Mean','Median', 'Lower_quartile', 'Upper_quartile', 'x10th_percentile', 'x90th_percentile', 'sample')
for (i in 1:length(sample_names)){
  x=qual_scores[grep(sample_names[i],qual_scores)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  x$sample=sample_names[i]
  colnames(x)=c('Base','Mean','Median', 'Lower_quartile', 'Upper_quartile', 'x10th_percentile', 'x90th_percentile', 'sample')
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
mr1=cbind(mr1,Read='Read 1')
d=mr1
p=d %>% filter(Read=='Read 1') %>% ggplot(aes(x = Base, ymin = x10th_percentile,
                                                lower = Lower_quartile, middle = Median, upper = Upper_quartile,
                                                ymax = x90th_percentile, color=sample, fill=sample, aes=0.5)) +
                                             geom_boxplot(stat='identity',position=position_dodge(0.7), alpha=0.3) +
                                             theme(axis.text.x = element_text(size = 8, angle=90, hjust=1, vjust=1),
                                                   axis.text.y=element_text(size=6),
                                                   plot.title = element_text(lineheight=.8, face="bold"),
                                                   axis.title=element_text(size=10)) +
                                             annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=20, alpha=0.1, fill="red") +
                                             annotate("rect", xmin=-Inf, xmax=Inf, ymin=20, ymax=28, alpha=0.1, fill="yellow") +
                                             annotate("rect", xmin=-Inf, xmax=Inf, ymin=28, ymax=Inf, alpha=0.1, fill="green") +
                                             xlab("Position of base in read") + ggtitle('Read 1') + 
                                             coord_cartesian(xlim=c(0,(max(d$Base)+1)))
ggsave(filename=paste0(outdir,'/per_base_sequence_quality_r1', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)

if (readType=='pairedEnd') {
  mr2=data.frame(matrix(ncol=8))
  colnames(mr2)=c('Base','Mean','Median', 'Lower_quartile', 'Upper_quartile', 'x10th_percentile', 'x90th_percentile', 'sample')
  for (i in 1:length(sample_names)){
    y=qual_scores[grep(sample_names[i],qual_scores)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    y$sample=sample_names[i]
    colnames(y)=c('Base','Mean','Median', 'Lower_quartile', 'Upper_quartile', 'x10th_percentile', 'x90th_percentile', 'sample')
    mr2=rbind(mr2,y)
  }
  mr2=mr2[-1,]
  mr2=cbind(mr2,Read='Read 2')

  d=rbind(mr1,mr2)
  p=d %>% filter(Read=='Read 2') %>% ggplot(aes(x = Base, ymin = x10th_percentile, 
                                                lower = Lower_quartile, middle = Median, upper = Upper_quartile,
                                                ymax = x90th_percentile, color=sample, fill=sample, aes=0.5)) +
                                             geom_boxplot(stat='identity',position=position_dodge(0.7), alpha=0.3) + 
                                             theme(axis.text.x = element_text(size = 8, angle=90, hjust=1, vjust=1), 
                                                   axis.text.y=element_text(size=6), 
                                                   plot.title = element_text(lineheight=.8, face="bold"), 
                                                   axis.title=element_text(size=10)) +
                                             annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=20, alpha=0.1, fill="red") +
                                             annotate("rect", xmin=-Inf, xmax=Inf, ymin=20, ymax=28, alpha=0.1, fill="yellow") +
                                             annotate("rect", xmin=-Inf, xmax=Inf, ymin=28, ymax=Inf, alpha=0.1, fill="green") + 
                                             xlab("Position of base in read") + ggtitle('Read 2') + 
                                             coord_cartesian(xlim=c(0,(max(d$Base)+1)))
  ggsave(filename=paste0(outdir,'/per_base_sequence_quality_r2', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)
}



# Per sequence quality scores
#----------------------------
qual_scores=files[grep('per_sequence_quality_scores.txt',files)]
mr1=matrix(ncol=3)
colnames(mr1)=c('x','y','Sample')
for (i in 1:length(sample_names)){
  x=qual_scores[grep(sample_names[i],qual_scores)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('x','y')
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
mr1=cbind(mr1,Read='Read 1')

if (readType=='pairedEnd') {
  mr2=matrix(ncol=3)
  colnames(mr2)=c('x','y','Sample')
  for (i in 1:length(sample_names)){
    y=qual_scores[grep(sample_names[i],qual_scores)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('x','y')
    y$Sample=sample_names[i]
    mr2=rbind(mr2,y)
  }
  mr2=mr2[-1,]
  mr2=cbind(mr2,Read='Read 2')
  
  d=rbind(mr1,mr2)
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=0.3) + geom_point(size=1) + facet_wrap(~Read) +
    theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), legend.text=element_text(size=8),  
           legend.key.height=unit(.8,"line"), axis.title.y=element_blank()) + ylab("") + xlab('Mean Quality score') + scale_shape_manual(values=shape_values[1:length(unique(d$Sample))]) 
} else {
  d=mr1
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=0.3) + geom_point(size=1) + facet_wrap(~Read) +
    theme( axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=7),  
           legend.key.height=unit(.8,"line"), axis.title.y=element_blank()) + 
           ylab("") + xlab('Mean Quality score') + scale_shape_manual(values=shape_values[1:length(unique(d$Sample))])
}

ggsave(filename=paste0(outdir,'/per_sequence_quality_scores', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)
print('Quality scores done... ')

# Per sequence gc content
#----------------------------
gc=files[grep('per_sequence_gc_content.txt',files)]
mr1=matrix(ncol=3)
colnames(mr1)=c('x','y','Sample')
for (i in 1:length(sample_names)){
  x=gc[grep(sample_names[i],gc)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('x','y')
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=cbind(mr1,Read='Read 1')
mr1=mr1[-1,]


if (readType=='pairedEnd') {
  mr2=matrix(ncol=3)
  colnames(mr2)=c('x','y','Sample')
  for (i in 1:length(sample_names)){
    y=gc[grep(sample_names[i],gc)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('x','y')
    y$Sample=sample_names[i]
    mr2=rbind(mr2,y)
  }
  mr2=cbind(mr2,Read='Read 2')
  mr2=mr2[-1,]
  
  d=rbind(mr1,mr2)
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) +  xlab('%GC') + ylab('') + theme(legend.key.height=unit(.8,"line")) + scale_shape_manual(values=shape_values[1:length(unique(d$Sample))])

} else {
  d=mr1
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) + 
    xlab('%GC') + ylab('') + theme(legend.key.height=unit(.8,"line")) +  scale_shape_manual(values=shape_values[1:length(unique(d$Sample))])

}
ggsave(filename=paste0(outdir,'/per_sequence_gc_content', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)
print('GC content done...')

# Per sequence length distribution
#---------------------------------
length_dist=files[grep('seq_length.txt',files)]
mr1=matrix(ncol=4)
colnames(mr1)=c('Position','Frequency','Sample', 'Read')
for (i in 1:length(sample_names)){
  x=length_dist[grep(sample_names[i],length_dist)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Position','Frequency')
  x$Position <- factor(x$Position, as.character(x$Position))
  x$Read='Read 1'
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
if (readType=='pairedEnd') {
  mr2=matrix(ncol=4)
  colnames(mr2)=c('Position','Frequency','Sample', 'Read')
  for (i in 1:length(sample_names)){
    y=length_dist[grep(sample_names[i],length_dist)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('Position','Frequency')
    y$Position <- factor(y$Position, as.character(y$Position))
    y$Read='Read 2'
    y$Sample=sample_names[i]
    mr2=rbind(mr2,y)
  }
  mr2=mr2[-1,]
  d=rbind(mr1,mr2)
  p=ggplot(d, aes(x = Position, y = Frequency, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) +
    theme( axis.title.x =element_text(size=12), axis.title.y =element_text(size=6), 
           axis.text.x=element_text(size=7,angle=90), axis.text.y=element_text(size=7))  + ylab("") + xlab('Length') + scale_shape_manual(values=scale_vaues[1:length(unique(d$Sample))])
} else {
  d=mr1
  p=ggplot(d, aes(x = Position, y = Frequency, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) +
    scale_shape_manual(values=shape_values[1:length(unique(d$Sample))]) + 
    theme( axis.title.x =element_text(size=12), axis.title.y =element_text(size=12), legend.key.height=unit(.8,"line"), 
           axis.text.x=element_text(size=7,angle=90), axis.text.y=element_text(size=12))  + ylab("") + xlab('length')
}
ggsave(filename=paste0(outdir,'/sequence_length_distribution', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)
print('Length distribution done')

# Per duplication levels
#-----------------------
dup_levels=files[grep('seq_dup_levels.txt',files)]
mr1=matrix(ncol=4)
colnames(mr1)=c('Duplication_Level','Percentage', 'Sample','Read')
for (i in 1:length(sample_names)){
  x=dup_levels[grep(sample_names[i],dup_levels)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Duplication_Level','Percentage of Deduplicated','Percentage')
  x$Duplication_Level <- factor(x$Duplication_Level, as.character(x$Duplication_Level))
  x=x[,-which(colnames(x)=="Percentage of Deduplicated")]
  x$Read='Read 1'
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
if (readType=='pairedEnd') {
  mr2=matrix(ncol=4)
  colnames(mr2)=c('Duplication_Level','Percentage', 'Sample','Read')
  for (i in 1:length(sample_names)){
    y=dup_levels[grep(sample_names[i],dup_levels)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('Duplication_Level','Percentage of Deduplicated','Percentage')
    y$Duplication_Level <- factor(y$Duplication_Level, as.character(y$Duplication_Level))
    y=y[,-which(colnames(y)=="Percentage of Deduplicated")]
    y$Read='Read 2'
    y$Sample=sample_names[i]
    mr2=rbind(mr2,y)
  }
  mr2=mr2[-1,]
  d=rbind(mr1,mr2)
  d$Duplication_Level=factor(d$Duplication_Level,levels=unique(d$Duplication_Level))
  p=ggplot(d, aes(x=Duplication_Level, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) + scale_shape_manual(values=shape_values[1:length(unique(d$Sample))]) +
    theme(axis.text.x = element_text(size = 9, angle=90), legend.key.height=unit(.8,"line"),
                        axis.title.y=element_blank())  + xlab("Number of copies per read")
} else {
  d=mr1
  d$Duplication_Level=factor(d$Duplication_Level,levels=unique(d$Duplication_Level))
  p=ggplot(d, aes(x=Duplication_Level, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) + 
    theme(axis.text.x = element_text(size = 9, angle=90), legend.key.height=unit(.8,"line"),
          axis.title.y=element_blank())  + xlab("Number of copies per read") + scale_shape_manual(values=shape_values[1:length(unique(d$Sample))])
}
ggsave(filename=paste0(outdir,'/sequence_dup_levels', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)
print('Duplication level done..')


# Base content
#-----------------------
base_content=files[grep('per_base_sequence_content.txt',files)]
mr1=matrix(ncol=5)
colnames(mr1)=c('Base', 'Percentage', 'Nucleotide','Sample','Read')
for (i in 1:length(sample_names)){
  x=base_content[grep(sample_names[i],base_content)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Base','G','A', 'T', 'C')
  all=matrix(ncol=3)
  colnames(all)=c('Base','Percentage','Nucleotide')
  nt=c('G','A','T','C')
  for (j in 1:length(nt)){
    all=rbind(all,cbind(x[,1], 'Percentage'=x[,which(colnames(x)==nt[j])], 'Nucleotide'=nt[j]))
  }
  all=as.data.frame(all[-1,])
  all$Percentage=as.numeric(as.character(all$Percentage))
  all$Read='Read 1'
  all$Sample=sample_names[i]
  mr1=rbind(mr1,all)
}
mr1=mr1[-1,]
mr1$Base=factor(mr1$Base, levels=unique(mr1$Base))
p=ggplot(mr1, aes(x=Base, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.2) + geom_point(size=0.8) + facet_wrap(~Nucleotide) + 
    theme(axis.text.x = element_text(size = 7, angle=90, v=0.5), axis.text.y = element_text(size = 7), legend.text=element_text(size=6), legend.title=element_text(size=7),
          axis.title=element_text(size=8), legend.key.height=unit(.8,"line"), strip.text=element_text(size=7, face='bold'), panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.2, linetype=2), panel.grid.minor=element_line(colour='grey',size=.2,linetype=2)) + 
    xlab("Position") + scale_shape_manual(values=shape_values[1:length(unique(mr1$Sample))]) + scale_x_discrete(breaks = c('1','25','50','75','100'))
ggsave(filename=paste0(outdir,'/per_base_sequence_content_r1', suffix, '.', plot_device), width=8, height=3.5, units='in', plot=p)

if (readType=='pairedEnd') {
    mr2=matrix(ncol=5)
    colnames(mr2)=c('Base', 'Percentage', 'Nucleotide','Sample','Read')
    for (i in 1:length(sample_names)){
      x=base_content[grep(sample_names[i],base_content)]
      x=x[grep('_R2_',x)]
      x=read.table(x, stringsAsFactors = FALSE)
      colnames(x)=c('Base','G','A', 'T', 'C')
      all=matrix(ncol=3)
      colnames(all)=c('Base','Percentage','Nucleotide')
      nt=c('G','A','T','C')
      for (j in 1:length(nt)){
        all=rbind(all,cbind(x[,1], 'Percentage'=x[,which(colnames(x)==nt[j])], 'Nucleotide'=nt[j]))
      }
      all=as.data.frame(all[-1,])
      all$Percentage=as.numeric(as.character(all$Percentage))
      all$Read='Read 2'
      all$Sample=sample_names[i]
      mr2=rbind(mr2,all)
    }
mr2=mr2[-1,]
mr2$Base=factor(mr2$Base, levels=unique(mr2$Base))
p=ggplot(mr2, aes(x=Base, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.2) + geom_point(size=0.8) + facet_wrap(~Nucleotide) + 
    theme(axis.text.x = element_text(size = 7, angle=90, v=0.5), axis.text.y = element_text(size = 7), legend.text=element_text(size=6), legend.title=element_text(size=7),
          axis.title=element_text(size=8), legend.key.height=unit(.8,"line"), strip.text=element_text(size=7, face='bold'), panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.2, linetype=2), panel.grid.minor=element_line(colour='grey',size=.2,linetype=2)) + 
    xlab("Position") + scale_shape_manual(values=shape_values[1:length(unique(mr2$Sample))]) + scale_x_discrete(breaks = c('1','25','50','75','100'))

}
ggsave(filename=paste0(outdir,'/per_base_sequence_content_r2', suffix, '.', plot_device), width=8, height=3.5, units='in', plot=p)
print('Base content done...')
