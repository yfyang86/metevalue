## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(metevalue)
#  
#  ####Simulation Data ####
#  set.seed(1234)
#  
#  # methylKit is a Bioconductor package
#  library(methylKit)
#  file.list=list( system.file("extdata",
#                              "test1.myCpG.txt", package = "methylKit"),
#                  system.file("extdata",
#                              "test2.myCpG.txt", package = "methylKit"),
#                  system.file("extdata",
#                              "control1.myCpG.txt", package = "methylKit"),
#                  system.file("extdata",
#                              "control2.myCpG.txt", package = "methylKit") )
#  
#  
#  # read the files to a methylRawList object: myobj
#  myobj=methRead(file.list,
#                 sample.id=list("test1","test2","ctrl1","ctrl2"),
#                 assembly="hg18",
#                 treatment=c(1,1,0,0),
#                 context="CpG"
#  )
#  
#  meth=unite(myobj, destrand=FALSE)
#  meth.C <- getData(meth)[,seq(6,ncol(meth),3)]
#  meth.T <- getData(meth)[,seq(7,ncol(meth),3)]
#  mr <- meth.C/(meth.C + meth.T)
#  chr_pos = getData(meth)[,1:2]
#  methyrate = data.frame(chr_pos,mr)
#  names(methyrate) = c('chr', 'pos', rep('g1',2), rep('g2',2))
#  region<-tileMethylCounts(myobj)
#  meth<-unite(region,destrand=F)
#  myDiff<-calculateDiffMeth(meth)
#  met_all<-getMethylDiff(myDiff,type="all")
#  
#  example_tempfiles = tempfile(c("rate_combine", "methylKit_DMR_raw"))
#  tempdir()
#  write.table(methyrate, file=example_tempfiles[1], row.names=F, col.names=T, quote=F, sep='\t')
#  write.table (met_all, file=example_tempfiles[2], sep ="\t", row.names =F, col.names =T, quote =F)

## ----eval=FALSE---------------------------------------------------------------
#  result = metevalue.methylKit(example_tempfiles[1], example_tempfiles[2], bheader = T)
#  str(result)

## ----eval=FALSE---------------------------------------------------------------
#  result = evalue_buildin_var_fmt_nm(methyrate, met_all, method="methylKit")
#  result = list(a = result$a,
#                b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b, method="methylKit"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  str(result)

## ----eval=FALSE---------------------------------------------------------------
#  library(BiSeq)
#  library(dplyr)
#  data(rrbs)
#  rrbs.rel <- rawToRel(rrbs)
#  methyrate <- methLevel(rrbs.rel)
#  methyrate <- data.frame(methyrate)
#  methyrateq = cbind(rows = as.numeric(row.names(methyrate)), methyrate)
#  methypos = data.frame(rows = as.numeric(row.names(methyrate)), rowRanges(rrbs))
#  methyrate = left_join(methypos, methyrateq)
#  methyrate = methyrate[,c(2,3,7:16)]
#  names(methyrate) <- c('chr','pos',rep('g1',5),rep('g2',5))
#  
#  rrbs.clust.unlim <- clusterSites(object = rrbs,perc.samples = 3/4,min.sites = 20,max.dist = 100)
#  
#  clusterSitesToGR(rrbs.clust.unlim)
#  ind.cov <- totalReads(rrbs.clust.unlim) > 0
#  
#  quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov])
#  rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)
#  predictedMeth <- predictMeth(object = rrbs.clust.lim)
#  
#  test<- predictedMeth[, colData(predictedMeth)$group == "test"]
#  control <- predictedMeth[, colData(predictedMeth)$group == "control"]
#  mean.test <- rowMeans(methLevel(test))
#  mean.control <- rowMeans(methLevel(control))
#  
#  betaResults <- betaRegression(formula = ~group,link = "probit",object = predictedMeth,type = "BR")
#  vario <- makeVariogram(betaResults)
#  vario.sm <- smoothVariogram(vario, sill = 0.9)
#  
#  locCor <- estLocCor(vario.sm)
#  clusters.rej <- testClusters(locCor)
#  clusters.trimmed <- trimClusters(clusters.rej)
#  DMRs <- findDMRs(clusters.trimmed,max.dist = 100,diff.dir = TRUE)
#  
#  
#  example_tempfiles = tempfile(c('rate_combine', 'BiSeq_DMR'))
#  write.table(methyrate, example_tempfiles[1], row.names=F, col.names=T, quote=F, sep='\t')
#  write.table(DMRs, example_tempfiles[2], quote=F, row.names = F,col.names = F, sep = '\t')

## ----eval=FALSE---------------------------------------------------------------
#  hh <- data.frame(DMRs)
#  result = evalue_buildin_var_fmt_nm(rate_combine, hh, method="biseq")
#  result = list(a = result$a,  b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b,method="biseq"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  str(result)

## ----eval=FALSE---------------------------------------------------------------
#  rate_combine <- read.table("rate_combine_DMRfinder", header = T)
#  head(rate_combine)
#  
#  DMRs <- read.table("DMRfinder_DMR", header = T)
#  head(DMRs)
#  

## ----eval=FALSE---------------------------------------------------------------
#  result <- evalue.DMRfinder('rate_combine_DMRfinder', 'DMRfinder_DMR')
#  head(result)

## ----eval=FALSE---------------------------------------------------------------
#  result = evalue_buildin_var_fmt_nm(rate_combine, DMRs, method="DMRfinder")
#  result = list(a = result$a,
#                b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b, method="DMRfinder"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  head(result)

## ----eval=FALSE---------------------------------------------------------------
#  input <- read.table("metilene.input", header = T)
#  head(input)
#  
#  out <- read.table("metilene.out", header = F)
#  head(out)
#  

## ----eval=FALSE---------------------------------------------------------------
#  result <- evalue.metilene('metilene.input', 'metilene.out')
#  head(result)

## ----eval=FALSE---------------------------------------------------------------
#  result = evalue_buildin_var_fmt_nm(input, out, method="metilene")
#  result = list(a = result$a,
#                b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b, method="metilene"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  head(result)

## ----evalue-------------------------------------------------------------------
#### Initialization ####
n_seeds = 100  # how many seeds to consider
N = 10000      # the number of trials
delta = -0.1   # the parameter of the alternative hypothesis
# e_x are the e-values and iE are the cumulative e-values
e_x = rep(0, N)
# p_x are the p-values and FP are Fisher's overall p-values
p_x = rep(0, N)
iE_all   = matrix(1, nrow = N+1, ncol = n_seeds) # product
uni_all  = iE_all                                # universal
FP_all   = iE_all                                # Fisher
F_VS_all = iE_all                                # Fisher-VS

#### Calculation ####
for(seed in 1:n_seeds){
  set.seed(seed * 1e3 + 1)
  x = rnorm(N) + delta
  e_x = exp(delta * x - delta^2/2)
  iE_all[, seed] = cumprod(c(1, e_x))
  S = cumsum(x)
  nn = 1:N
  uni_all[-1,seed] = exp(S^2/2/nn) / sqrt(nn)
  p_x = pnorm(x)
  FF = -2 * cumsum(log(p_x))
  FP_all[-1, seed] = exp(pchisq(FF, df=2*nn, lower.tail = F, log.p = T))
  SELR_ = (FP_all[, seed] < exp(-1))
  F_VS_all[ SELR_, seed] = 1/(-exp(1)*FP_all[ SELR_, seed]*log(FP_all[ SELR_, seed]))   
}

iE = apply(iE_all, 1, median)
uni = apply(uni_all, 1, median)
FP = apply(FP_all, 1, median)
F_VS = apply(F_VS_all,1, median)

#### Plot ####
library(ggplot2)
library(tidyr)
library(dplyr)

sim_plots <- data.frame(
  x = 1:(N+1),
  product = iE,
  universal = uni,
  Fisher = 1/(FP),
  FisherVS = F_VS
) %>% 
  gather(key = "Method", value = "E_value", -x) %>%
  ggplot(aes(x = x, y = E_value)) + 
  geom_line(aes(color = Method), size = 1) + 
  scale_y_continuous(trans='log') +
  theme_grey() +  # Default
  theme(legend.position = "top") + 
  scale_color_brewer(palette="Dark2") +
  xlab("Numer of Observations") + 
  ylab("e-Value")
  
print(sim_plots)

