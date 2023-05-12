library(dplyr)
arg = commandArgs(T)
num = as.numeric(arg[1])
dir = arg[2]
for(n in 1:num){
  if(n==1){
    data = read.table(paste0(dir,"/control",n,'_rate.bed'), header=F, sep='\t')
    data = data[,-2]
    names(data) = c('chr','pos',paste0('rate',n))
  }
  else{
    data2 = read.table(paste0(dir,"/control",n,'_rate.bed'), header=F, sep='\t')
    data2 = data2[,-2]
    names(data2) = c('chr','pos',paste0('rate',n))
    data = full_join(data, data2)
  }
}
for(n in 1:num){
  data2 = read.table(paste0(dir,"/test",n,'_rate.bed'), header=F, sep='\t')
  data2 = data2[,-2]
  names(data2) = c('chr','pos',paste0('rate',8+n))
  data = left_join(data, data2)
}
data[,3:(2+num*2)] = data[,3:(2+num*2)]/100
names(data) = c('chr','pos',rep('g1',num),rep('g2',num))
outputdir = paste0(dir,'/rate_combine_',num)
write.table(data, outputdir, row.names=F, col.names=T, quote=F, sep='\t')
