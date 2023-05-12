#!/usr/bin/env Rscript
#start<-proc.time()
library(methylKit)

arg = commandArgs(T)
n = as.numeric(arg[1])
difference = as.numeric(arg[2])
data_dir = arg[3]
output_dir = arg[4]

temp1<-paste0(data_dir,"/test",1:n,".sorted.bam")
temp2<-paste0(data_dir,"/control",1:n,".sorted.bam")
temp<-c(temp1,temp2)
bam.list<-list()
n2<-2*n
for (i in 1:n2){
    bam.list[[i]]<-temp[i]
}

sample.id<-list()
for (i in 1:n){
    sample.id[[i]]<-paste0("test",i)
}
for (i in 1:n){
    sample.id[[i+n]]<-paste0("control",i)
}

treatment<-c(rep(1,n),rep(0,n))


#bam.list=list("test1.sorted.bam","test2.sorted.bam","test3.sorted.bam","test4.sorted.bam","test5.sorted.bam","test6.sorted.bam","test7.sorted.bam","test8.sorted.bam","control1.sorted.bam","control2.sorted.bam","control3.sorted.bam","control4.sorted.bam","control5.sorted.bam","control6.sorted.bam","control7.sorted.bam","control8.sorted.bam")
myobj=processBismarkAln(bam.list,sample.id=sample.id,assembly="hg19",treatment=treatment,read.context="CpG",save.folder=getwd())

#取并集
#meth.min=unite(myobj,min.per.group=1L) 

region<-tileMethylCounts(myobj)
meth<-unite(region,destrand=F)
myDiff<-calculateDiffMeth(meth)
differ = difference*100
all<-getMethylDiff(myDiff,qvalue=1,type="all",difference=differ)
#write.table (all, file=paste0(args[2],"/methylKit_DMR"), sep ="\t", row.names =F, col.names =T, quote =F)
write.table(all, file=output_dir, sep ="\t", row.names =F, col.names =T, quote =F)
#end<-proc.time() - start 

#output<-data.frame(name="methylKit",dir=basename(getwd()),user=end["user.self"],sys=end["sys.self"],elapsed=end["elapsed"])
#write.table(output,paste0("time_record_collection_",n,".",time),col.names =FALSE,row.names =FALSE,sep="\t",append=TRUE,quote=FALSE)