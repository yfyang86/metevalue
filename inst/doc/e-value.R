## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(evalue)

####Simulation Data ####
set.seed(1234)

simu_g_value <- function(n, r = 0.1){
  x = runif(n)
  x[runif(n) <= r] = 0
  return(x);
}


library(methylKit)
file.list=list( system.file("extdata", 
                            "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata",
                            "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", 
                            "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", 
                            "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
               sample.id=list("test1","test2","ctrl1","ctrl2"),
               assembly="hg18",
               treatment=c(1,1,0,0),
               context="CpG"
)

meth=unite(myobj, destrand=FALSE)
meth.C <- getData(meth)[,seq(6,ncol(meth),3)]
meth.T <- getData(meth)[,seq(7,ncol(meth),3)]
mr <- meth.C/(meth.C + meth.T)
chr_pos = getData(meth)[,1:2]
methyrate = data.frame(chr_pos,mr)
names(methyrate) = c('chr', 'pos', rep('g1',2), rep('g2',2))
region<-tileMethylCounts(myobj)
meth<-unite(region,destrand=F)
myDiff<-calculateDiffMeth(meth)
met_all<-getMethylDiff(myDiff,type="all")

example_tempfiles = tempfile(c("rate_combine", "methylKit_DMR_raw"))
tempdir()
write.table(methyrate, file=example_tempfiles[1], row.names=F, col.names=T, quote=F, sep='\t')
write.table (met_all, file=example_tempfiles[2], sep ="\t", row.names =F, col.names =T, quote =F)

## ----methylkit----------------------------------------------------------------
result = evalue.methylKit(example_tempfiles[1], example_tempfiles[2], bheader = T)
str(result)

## ----methylkit-function-------------------------------------------------------
result = evalue_buildin_var_fmt_nm(methyrate, met_all, method="methylKit")
result = list(a = result$a, 
              b = result$b, 
              a_b = evalue_buildin_sql(result$a, result$b, method="methylKit"))
result = varevalue.metilene(result$a, result$b, result$a_b)
str(result)

