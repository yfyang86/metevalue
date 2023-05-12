##############################
library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

library("DESeq2")
colnames(cts)=paste0(colnames(cts),'fb')
cts = cts[,rownames(coldata)]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)


dat <- t(t(cts)/(dds$sizeFactor)) 
dat.out <- dat[rowSums(dat >5)>=0.8*ncol(dat),]

datasim <- log(dat.out)

######### add your e-value algorithm, row by row; we don't have regions in RNA-seq

library(stringr)
evalue_RNA = function(a_b, group1_name, group2_name){
  a_b = data.frame(a_b)
  innerf = function(x, innermu=0., innersig=1.){
    vector_temp = na.omit(as.numeric(x))
    n = length(vector_temp)
    value = mean(vector_temp)
    dnorm(x=value, mean = innermu, sd = innersig/sqrt(n))
  }
  
  innerlog = function(x){
    log(x[!is.na(x)])
  }
  
  site_1 = grep(paste0('^',group1_name),colnames(a_b),value=T)
  site_2 = grep(paste0('^',group2_name),colnames(a_b),value=T)
  
  evalue_all = c()
  for(i in rownames(a_b)){
  
    a_b_temp = a_b[i, ]
  
    vector_1 = unlist(a_b_temp[, site_1])
    vector_2 = unlist(a_b_temp[, site_2])
    vector = c(vector_1, vector_2)
    miu_1 = mean(vector_1, na.rm=T)
    sigma_1 = sd(vector_1, na.rm=T)
    miu_2 = mean(vector_2, na.rm=T)
    sigma_2 = sd(vector_2, na.rm=T)
    miu = mean(vector, na.rm=T)
    sigma = sd(vector, na.rm=T)
  
    norm_value_up_1 = apply(a_b_temp[, site_1], 2, function(x)innerf(x=x,innermu=miu_1, innersig=sigma_1))
    norm_value_down_1 = apply(a_b_temp[, site_1], 2, function(x)innerf(x=x, innermu=miu, innersig=sigma))
    norm_value_up_2 = apply(a_b_temp[, site_2], 2,  function(x)innerf(x=x, innermu=miu_2, innersig=sigma_2))
    norm_value_down_2 = apply(a_b_temp[, site_2], 2,  function(x)innerf(x=x, innermu=miu, innersig=sigma))
  
    e_value = exp(sum(innerlog(c(norm_value_up_1,norm_value_up_2))) -
                    sum(innerlog(c(norm_value_down_1,norm_value_down_2))))
    evalue_all = c(evalue_all,e_value)
  }
  return(cbind(a_b,evalue_all))
}
evalue = evalue_RNA(datasim, 'treated','untreated')
