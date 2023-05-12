#!/usr/bin/env Rscript
#start<-proc.time()

library(BiSeq)

arg = commandArgs(T)
n = as.numeric(arg[1])
difference = as.numeric(arg[2])
data_dir = arg[3]
output_dir = arg[4]


for (i in 1:n){
    file_name = paste0(data_dir,"/test",i,".bed")
    x_name = paste0("test",i)
    assign(x_name,readBismark(file_name,colData=DataFrame(row.names=paste0("test_",i),group = c("test"))))
}


for (i in 1:n){
    file_name = paste0(data_dir,"/control",i,".bed")
    x_name = paste0("control",i)
    assign(x_name,readBismark(file_name,colData=DataFrame(row.names=paste0("control_",i),group = c("control"))))
}



if(n<=2){
temp<-combine(test1,test2)
for (i in 1:n){
   temp<-combine(temp,get(paste0("control",i)))
}
}else{
temp<-combine(test1,test2)
for (i in 3:n){
   temp<-combine(temp,get(paste0("test",i)))
}

for (i in 1:n){
   temp<-combine(temp,get(paste0("control",i)))
}
}

#covStatistics(temp)
rrbs=temp
#covBoxplots(rrbs, col = "cornflowerblue", las = 2)
rrbs.clust.unlim <- clusterSites(object = rrbs,perc.samples = 3/4,min.sites = 20,max.dist = 100)
#head(rowRanges(rrbs.clust.unlim))
clusterSitesToGR(rrbs.clust.unlim)
ind.cov <- totalReads(rrbs.clust.unlim) > 0
# quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov],0.9)
quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov])
rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)
predictedMeth <- predictMeth(object = rrbs.clust.lim)
#covBoxplots(rrbs.clust.lim, col = "cornflowerblue", las = 2)
test<- predictedMeth[, colData(predictedMeth)$group == "test"]
control <- predictedMeth[, colData(predictedMeth)$group == "control"]
mean.test <- rowMeans(methLevel(test))
mean.control <- rowMeans(methLevel(control))
#plot(mean.control, mean.test, col = "blue", xlab = "Methylation in controls", ylab = "Methylation in tests")
betaResults <- betaRegression(formula = ~group,link = "probit",object = predictedMeth,type = "BR")
vario <- makeVariogram(betaResults)
vario.sm <- smoothVariogram(vario, sill = 0.9)
# vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)
# vario.sm$pValsList <- vario.aux$pValsList
locCor <- estLocCor(vario.sm)
clusters.rej <- testClusters(locCor)
# clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = 0.1)
clusters.trimmed <- trimClusters(clusters.rej)
#head(clusters.trimmed)
DMRs <- findDMRs(clusters.trimmed,max.dist = 100,diff.dir = TRUE)
DMRss <- DMRs[abs(DMRs$median.meth.diff)>=difference,]
write.table(DMRss, output_dir, quote=F, row.names = F,col.names = T,sep="\t")

# end<-proc.time() - start 
# output<-data.frame(name="BiSeq",dir=basename(getwd()),user=end["user.self"],sys=end["sys.self"],elapsed=end["elapsed"])
# write.table(output,paste0("time_record_collection_",n,".",time),col.names =FALSE,row.names =FALSE,sep="\t",append=TRUE,quote=FALSE)
