#!/usr/local/bin/Rscript

in_dir = commandArgs(TRUE)[1]
out_dir = commandArgs(TRUE)[2]

system(paste0("egrep -w \"Number of input reads | Uniquely mapped reads | Number of reads mapped to multiple loci | % of reads mapped to multiple loci |Number of reads mapped to too many loci | % of reads mapped to too many loci | % of reads unmapped: too many mismatches | % of reads unmapped: too short | % of reads unmapped: other | Uniquely mapped reads number | Number of splices:\" ", in_dir, "/*Log.final.out > star_logout.txt"))

system("sed 's/out: /out:\t/g' star_logout.txt  | sed 's/ /_/g' | sed 's/ /\t/g' | awk '{print $1,$2,$3}' > star_logout2.txt && mv star_logout2.txt star_logout.txt")

x=read.table("star_logout.txt", stringsAsFactors=FALSE)
names=gsub('.*/','',x[,1])
#names=strsplit(x[,1],"/")
#names=sapply(names,"[",2)
x[,1]=names
x[,1]=gsub("Log.final.out:","",x[,1])
x[,2]=gsub("__","",x[,2])
x[,2]=gsub("_ |","",x[,2])
#x[,2]=gsub("^_","",x[,2])
x[,3]=gsub("%","",x[,3])
star=x
names=unique(x[,1])
star=split(star,star$V2)
out=as.data.frame(cbind(Uniquely_Mapped_pct=as.numeric(star[[5]][match(names, star[[5]][,1]),3]),
                        Multiple_Loci_pct=as.numeric(star[[6]][match(names, star[[6]][,1]),3])+
                          as.numeric(star[[7]][match(names, star[[7]][,1]),3]),
                        Unmapped_pct=as.numeric(star[[8]][match(names, star[[8]][,1]),3])+
                          as.numeric(star[[9]][match(names, star[[9]][,1]),3])+
                          as.numeric(star[[10]][match(names, star[[10]][,1]),3]), 
			Input_num=as.numeric(star[[1]][match(names, star[[1]][,1]),3]),
			Mapped_num=as.numeric(star[[16]][match(names, star[[16]][,1]),3])))

rownames(out)=names
write.csv(out,paste0(out_dir,"/mapping_summary.csv"), quote=FALSE)
      

x=read.table("star_logout.txt", stringsAsFactors=FALSE)
names=gsub('.*/','',x[,1])
x[,1]=names
x[,1]=gsub("Log.final.out:","",x[,1])
x[,2]=gsub("__","",x[,2])
x[,2]=gsub("_ |","",x[,2])
x[,3]=gsub("%","",x[,3])
star=x
names=unique(x[,1])
star=split(star,star$V2)
out=as.data.frame(cbind(Splices_num=as.numeric(star[[9]][match(names, star[[9]][,1]),3]),
                        Annotated_splices_num=as.numeric(star[[4]][match(names, star[[4]][,1]),3])))
names=gsub('Log.final.out:','',names)
rownames(out)=names
write.csv(out,paste0(out_dir,"mapping_splices.csv"), quote=FALSE)
                        
system("rm star_logout.txt")
#Done
