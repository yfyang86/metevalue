
# Introductions

In this package, we provide e-value for four DMR (differentially methylated region) detection tools (MethylKit, Metilene, BiSeq and DMRfinder) and general purpose.

-   MethylKit
-   BiSeq
-   DMRfinder
-   Metilene
-   Other DNA methylation tools
-   RNA-seq data

For `DMR` (`methylKit`, `biseq`, `DMRfinder` or `metilene`), the met-evalue calculation is conducted by the `metevalue.[DMR]` function. 

| DMR | Method | Input.1 Example | Input.2 Example |
|:-----|:-----|:-----|:-----|
| MethylKit | `metevalue.methylKit` | `data(demo_methylkit_methyrate)` |  `data(demo_methylkit_met_all)` | 
| BiSeq | `metevalue.biseq` | `data(demo_biseq_methyrate)` | `data(demo_biseq_DMR)` |
| DMRfinder | `metevalue.DMRfinder`|  `data(demo_DMRfinder_rate_combine)` | `data(demo_DMRfinder_DMRs)` |
| Metilene | `metevalue.metilene` | `data(demo_metilene_input)` | `data(demo_metilene_out)` |
| Other DNA methylation tools | `varevalue.single_general` | `data(demo_metilene_input)` or any data above| | 
| RNA-seq data | `metevalue.RNA_general` | `data(demo_desq_out)` | | 

Two routines are supported to calculate the combined e-value:

- Call by **files**: Here the `files` include the outputs of given `DMR` packages and its corresponding e-value of each region;
- Call by **R data frames**: Here the `R data frames` are corresponding `data.frame` objects.


# Other Demos

Please vist the [metevalue-emo](https://github.com/yfyang86/metevalue-demos/tree/main/simulation) project for more demos.


## Call by files

We design the `metevalue.[DMR]` function to accept similar parameter patterns:

``` r
metevalue.[DMR](
  methyrate,                # methylation rates of each CpG site
  [DMR].output,             # Output file name of [DMR] with e-value of each region
  adjust.methods = "BH",    # Adjust methods of e-value
  sep = "\t",               # seperator, default is the TAB key
  bheader = FALSE           # A logical value indicating whether the [DMR].output file
                            # contains the names of the variables as its first line
)
```

Here  `[DMR]` could be one of `methylKit`, `biseq`, `DMRfinder` or `metilene`.

## Call by R data frames

We provide the `evalue_buildin_var_fmt_nm` and `varevalue.metilene` function to handle the general DMR e-value calculation in DNA methylation studies:

``` r
# Here  `[DMR]` coudle be one of `methylKit`, `biseq`, `DMRfinder` or `metilene`.
method_in_use = "[DMR]"
result = evalue_buildin_var_fmt_nm(
          methyrate,              # Data frame of the methylation rate
          DMR_evalue_output,      # Data frame of output data corresponding to the
                                  # "method" option
          method = method_in_use) # DMR: "metilene", "biseq", "DMRfinder" or "methylKit"
result = list(a = result$a,
              b = result$b,
              a_b = evalue_buildin_sql(result$a, result$b, method = method_in_use))
result = varevalue.metilene(result$a, result$b, result$a_b)
```

Replace `[DMR]` to one of `methylKit`, `biseq`, `DMRfinder` or `metilene` accordingly.

For `RNAseq` user, `metevalue.RNA_general` could be called directly. Example is:

```r
data("demo_desq_out")
evalue = metevalue.RNA_general(demo_desq_out, 'treated','untreated')
```

> Notice: for different `[DMR]`, the `data.frame` schemas are **different**!!! Check the R help document for details. Check the [Demo data](#demo-data) section for details.

## Example: MethylKit

methylKit is an R package for DNA methylation analysis and annotation
from high-throughput bisulfite sequencing. The package is designed to
deal with sequencing data from RRBS and its variants, but also
target-capture methods and whole genome bisulfite sequencing.

Currently, `metevalue` package supports the e-value calculation using the
`methylKit` output file.

``` r
library(metevalue)

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
#> two groups detected:
#>    will calculate methylation difference as the difference of
#>    treatment (group: 1) - control (group: 0)
met_all<-getMethylDiff(myDiff,type="all")

example_tempfiles = tempfile(c("rate_combine", "methylKit_DMR_raw"))
tempdir()
write.table(methyrate, file=example_tempfiles[1], row.names=F, col.names=T, quote=F, sep='\t')
write.table (met_all, file=example_tempfiles[2], sep ="\t", row.names =F, col.names =T, quote =F)
```

`evalue.methylKit` function could be used to tackle the problem.

``` r
result = metevalue.methylKit(example_tempfiles[1], example_tempfiles[2], bheader = T)
#> Joining, by = c("start", "end")
str(result)
#> 'data.frame':    24 obs. of  9 variables:
#>  $ chr      : chr  "chr21" "chr21" "chr21" "chr21" ...
#>  $ start    : int  9927001 9944001 9959001 9967001 10011001 10077001 10087001 10186001 13664001 13991001 ...
#>  $ end      : int  9928000 9945000 9960000 9968000 10012000 10078000 10088000 10187000 13665000 13992000 ...
#>  $ strand   : chr  "*" "*" "*" "*" ...
#>  $ p        : num  2.47e-10 2.57e-21 4.39e-23 3.08e-04 2.02e-65 ...
#>  $ qvalue   : num  3.24e-10 9.58e-21 2.36e-22 2.37e-04 3.27e-64 ...
#>  $ meth.diff: num  -34.1 -40.2 -25.4 -25.9 25.8 ...
#>  $ e_value  : num  1.65 1.65 1.65 1.65 1.65 ...
#>  $ e_adjust : num  1.65 1.65 1.65 1.65 1.65 ...
```

Alternatively, one could use the build-in functions to derive functions
which avoid the file operation:

``` r
result = evalue_buildin_var_fmt_nm(methyrate, met_all, method="methylKit")
result = list(a = result$a, 
              b = result$b, 
              a_b = evalue_buildin_sql(result$a, result$b, method="methylKit"))
result = varevalue.metilene(result$a, result$b, result$a_b)
#> Joining, by = c("start", "end")
str(result)
#> 'data.frame':    24 obs. of  9 variables:
#>  $ chr      : Factor w/ 1 level "chr21": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ start    : int  9927001 9944001 9959001 9967001 10011001 10077001 10087001 10186001 13664001 13991001 ...
#>  $ end      : int  9928000 9945000 9960000 9968000 10012000 10078000 10088000 10187000 13665000 13992000 ...
#>  $ strand   : Factor w/ 3 levels "+","-","*": 3 3 3 3 3 3 3 3 3 3 ...
#>  $ p        : num  2.47e-10 2.57e-21 4.39e-23 3.08e-04 2.02e-65 ...
#>  $ qvalue   : num  3.24e-10 9.58e-21 2.36e-22 2.37e-04 3.27e-64 ...
#>  $ meth.diff: num  -34.1 -40.2 -25.4 -25.9 25.8 ...
#>  $ e_value  : num  1.65 1.65 1.65 1.65 1.65 ...
#>  $ e_adjust : num  1.65 1.65 1.65 1.65 1.65 ...
```

## Example: BiSeq

First, we load the methylation data at CpG site levels from ‘BiSeq’
package. Then we cluster CpG sites into DMRs using ‘BiSeq’.

``` r
library(BiSeq)
library(dplyr)
data(rrbs)
rrbs.rel <- rawToRel(rrbs)
methyrate <- methLevel(rrbs.rel)
methyrate <- data.frame(methyrate)
methyrateq = cbind(rows = as.numeric(row.names(methyrate)), methyrate)
methypos = data.frame(rows = as.numeric(row.names(methyrate)), rowRanges(rrbs))
methyrate = left_join(methypos, methyrateq)
methyrate = methyrate[,c(2,3,7:16)]
names(methyrate) <- c('chr','pos',rep('g1',5),rep('g2',5))

rrbs.clust.unlim <- clusterSites(object = rrbs,perc.samples = 3/4,min.sites = 20,max.dist = 100)

clusterSitesToGR(rrbs.clust.unlim)
ind.cov <- totalReads(rrbs.clust.unlim) > 0

quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov])
rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)
predictedMeth <- predictMeth(object = rrbs.clust.lim)

test<- predictedMeth[, colData(predictedMeth)$group == "test"]
control <- predictedMeth[, colData(predictedMeth)$group == "control"]
mean.test <- rowMeans(methLevel(test))
mean.control <- rowMeans(methLevel(control))

betaResults <- betaRegression(formula = ~group,link = "probit",object = predictedMeth,type = "BR")
vario <- makeVariogram(betaResults)
vario.sm <- smoothVariogram(vario, sill = 0.9)

locCor <- estLocCor(vario.sm)
clusters.rej <- testClusters(locCor)
clusters.trimmed <- trimClusters(clusters.rej)
DMRs <- findDMRs(clusters.trimmed,max.dist = 100,diff.dir = TRUE)


example_tempfiles = tempfile(c('rate_combine', 'BiSeq_DMR'))
write.table(methyrate, example_tempfiles[1], row.names=F, col.names=T, quote=F, sep='\t')
write.table(DMRs, example_tempfiles[2], quote=F, row.names = F,col.names = F, sep = '\t')
```

Finally, we add E-values and adjusted E-values as additional columns
to the output file of ‘BiSeq’.`metevalue.biseq` function could be used to
tackle the problem.

``` r
result = metevalue.biseq(example_tempfiles[1], example_tempfiles[2])
str(result)
```

## Example: DMRfinder

Given the input file

-   `rate_combine_DMRfinder`: a file containing methylation rates at each CpG site

-   `DMRfinder_DMR`: the output file from ‘DMRfinder’

``` r
rate_combine <- read.table("rate_combine_DMRfinder", header = T)
head(rate_combine)

DMRs <- read.table("DMRfinder_DMR", header = T)
head(DMRs)
```

Adding E-values and adjusted E-values as additional columns to file
‘DMRfinder_DMR’

``` r
result <- metevalue.DMRfinder('rate_combine_DMRfinder', 'DMRfinder_DMR', bheader=T)
head(result)
```

Alternatively, function `varevalue.metilene` can also provide e-value
and adjusted e-value.

``` r
result = evalue_buildin_var_fmt_nm(rate_combine, DMRs, method="DMRfinder")
result = list(a = result$a, 
              b = result$b, 
              a_b = evalue_buildin_sql(result$a, result$b, method="DMRfinder"))
result = varevalue.metilene(result$a, result$b, result$a_b)
head(result)
```

## Example: Metilene

Given

-   `metilene.input`: the input file of `Metilene` containing methylation rates at each CpG site
-   `metilene.out`: the output file of `Metilene`

``` r
input <- read.table("metilene.input", header = T)
head(input)

out <- read.table("metilene.out", header = F)
head(out)
```

Adding E-values and adjusted E-values as additional columns to
`metilene.out`

``` r
result <- metevalue.metilene('metilene.input', 'metilene.out')
head(result)
```

Alternatively, function `varevalue.metilene` can also provide e-value
and adjusted e-value.

``` r
result = evalue_buildin_var_fmt_nm(input, out, method="metilene")
result = list(a = result$a, 
              b = result$b, 
              a_b = evalue_buildin_sql(result$a, result$b, method="metilene"))
result = varevalue.metilene(result$a, result$b, result$a_b)
head(result)
```

## Example: Other DNA methylation tools

In above examples, we have already provided examples to calculate E-values directly from DMR detection tools including BiSeq, DMRfinder, MethylKit and Metilene. All of these require users to prepare an output file of different tools.
However, users may wonder how to calculate the E-values directly from CpG sites or other DNA methylation tools not presented above.
We then facilitate the purpose in the following example.


-   `methyrate`: a file containing methylation rates at each CpG site of 2 different groups

By changing the group name, start site and end site, function `varevalue.single_general` can calculate e-value of any site or region using a general methylation rates data without using an output file of a specific tool.

```{r eval=FALSE}
input <- read.table("methyrate", header = T)
e_value <- varevalue.single_general(methyrate=input, group1_name='g1', group2_name='g2', chr='chr21', start=9439679, end=9439679)
head(e_value)
```

## Example: RNA-seq data
The framework of E-value calculation presented in this project is also able to be extended to other genomic data including RNA-seq. 
Here is an example to introduce the E-value calculation in RNA-seq. 


-   `desq_out`: the RNA data

function `metevalue.RNA_general` can provide e-values for each row of the normalized expression level of RNA-seq data.

```{r eval=FALSE}
input <- read.table("desq_out", header = T)
data_e <- metevalue.RNA_general(input, group1_name='treated', group2_name='untreated')
head(data_e)
```



# Misc

## Demo data {#demo-data}

Demo data for different `metevalue.[DMR]` functions are listed in the section.

### Input Data Examples: MethylKit

**methyrate Example**

|chr   |     pos|        g1|        g1|        g2|        g2|
|:-----|-------:|---------:|---------:|---------:|---------:|
|chr21 | 9853296| 0.5882353| 0.8048048| 0.8888889| 0.8632911|
|chr21 | 9853326| 0.7058824| 0.7591463| 0.8750000| 0.7493404|


**methylKit.output Example**

|chr| start|     end| strand | pvalue| qvalue|meth.diff|
|:-----|-------:|-------:|:------|------:|------:|---------:|
|chr21 | 9927001| 9928000|*      |      0|      0| -34.07557|
|chr21 | 9944001| 9945000|*      |      0|      0| -40.19089|

### Input Data Examples: BiSeq

**methyrate Example**

|chr  |    pos|        g1| g1|  g1|  g1|  g1|     g2|        g2|        g2|        g2|        g2|
|:----|------:|---------:|--:|---:|---:|---:|------:|---------:|---------:|---------:|---------:|
|chr1 | 870425| 0.8205128|  1| 0.7| NaN| NaN| 0.3125| 0.7419355| 0.2461538| 0.1794872| 0.2413793|
|chr1 | 870443| 0.8461538|  1| 0.7| NaN| NaN| 0.3750| 0.3225806| 0.2923077| 0.0512821| 0.2413793|


**biseq.output Example**

|seqnames |  start|    end| width|strand |  median.p| median.meth.group1| median.meth.group2| median.meth.diff|
|:--------|------:|------:|-----:|:------|---------:|------------------:|------------------:|----------------:|
|chr1     | 872369| 872616|   248|*      | 0.0753559|          0.9385462|          0.8666990|        0.0710524|
|chr1     | 875227| 875470|   244|*      | 0.0000026|          0.5136315|          0.1991452|        0.2942668|

### Input Data Examples: DMRfinder

**methyrate Example**

|chr  |       pos| g1|      g1.1| g2| g2.1|
|:----|---------:|--:|---------:|--:|----:|
|chr1 | 202833315|  0| 0.0000000|  0|    0|
|chr1 | 202833323|  1| 0.8095238|  1|    1|

**DMRfinder.output Example**

|chr   |    start|      end| CpG| Control.mu|  Exptl.mu| Control..Exptl.diff| Control..Exptl.pval|
|:-----|--------:|--------:|---:|----------:|---------:|-------------------:|-------------------:|
|chr8  | 25164078| 25164102|   3|  0.9241646| 0.7803819|          -0.1437827|           0.0333849|
|chr21 |  9437432|  9437538|  14|  0.7216685| 0.1215506|          -0.6001179|           0.0000000|

### Input Data Examples: DMRfinder

**methyrate Example**

|chr |     pos|        g1| g1.1|      g1.2| g1.3| g1.4|      g1.5| g1.6|      g1.7|        g2| g2.1| g2.2| g2.3|      g2.4| g2.5| g2.6| g2.7|
|:-----|-------:|---------:|----:|---------:|----:|----:|---------:|----:|---------:|---------:|----:|----:|----:|---------:|----:|----:|----:|
|chr21 | 9437433| 0.9285714|   NA| 0.7222222| 0.75|    1| 0.6666667|    1| 0.8695652| 0.0000000|    0|    0|    0| 0.0000000|  0.0|   NA| 0.00|
|chr21 | 9437445| 1.0000000|   NA| 0.9444444| 0.75|    1| 0.6666667|    0| 0.8695652| 0.6111111|    0|    0|    0| 0.7333333|  0.6|   NA| 0.75|

**metilene.output Example**

| chr  | start  |  end  |  q-value | methyl.diff  |  CpGs  |  p | p2  | m1 | m2 |
|:-----|-------:|-------:|--:|--------:|--:|--:|--:|-------:|-------:|
|chr21 | 9437432| 9437540|  0| 0.610989| 26|  0|  0| 0.73705| 0.12606|
|chr21 | 9708982| 9709189|  0| 0.475630| 28|  0|  0| 0.58862| 0.11299|


### Input Data Examples: Metilene

**metilene.input Example**

|chr |     pos|        g1| g1.1|      g1.2| g1.3| g1.4|      g1.5| g1.6|      g1.7|        g2| g2.1| g2.2| g2.3|      g2.4| g2.5| g2.6| g2.7|
|:-----|-------:|---------:|----:|---------:|----:|----:|---------:|----:|---------:|---------:|----:|----:|----:|---------:|----:|----:|----:|
|chr21 | 9437433| 0.9285714|   NA| 0.7222222| 0.75|    1| 0.6666667|    1| 0.8695652| 0.0000000|    0|    0|    0| 0.0000000|  0.0|   NA| 0.00|
|chr21 | 9437445| 1.0000000|   NA| 0.9444444| 0.75|    1| 0.6666667|    0| 0.8695652| 0.6111111|    0|    0|    0| 0.7333333|  0.6|   NA| 0.75|

**metilene.output Example**

| chr  | start  |  end  |  q-value | methyl.diff  |  CpGs  |  p | p2  | m1 | m2 |
|:-----|-------:|-------:|--:|--------:|--:|--:|--:|-------:|-------:|
|chr21 | 9437432| 9437540|  0| 0.610989| 26|  0|  0| 0.73705| 0.12606|
|chr21 | 9708982| 9709189|  0| 0.475630| 28|  0|  0| 0.58862| 0.11299|

### Input Data Examples: Other DNA methylation tools

**methyrate Example**

|chr |     pos|        g1| g1.1|      g1.2| g1.3| g1.4|      g1.5| g1.6|      g1.7|        g2| g2.1| g2.2| g2.3|      g2.4| g2.5| g2.6| g2.7|
|:-----|-------:|---------:|----:|---------:|----:|----:|---------:|----:|---------:|---------:|----:|----:|----:|---------:|----:|----:|----:|
|chr21 | 9437433| 0.9285714|   NA| 0.7222222| 0.75|    1| 0.6666667|    1| 0.8695652| 0.0000000|    0|    0|    0| 0.0000000|  0.0|   NA| 0.00|
|chr21 | 9437445| 1.0000000|   NA| 0.9444444| 0.75|    1| 0.6666667|    0| 0.8695652| 0.6111111|    0|    0|    0| 0.7333333|  0.6|   NA| 0.75|


### Input Data Examples: RNA-seq data

**desq_out Example**

| treated1fb| treated2fb| treated3fb| untreated1fb| untreated2fb| untreated3fb| untreated4fb|
|:----------|:----------|:----------|:------------|:------------|:------------|:------------|
|   4.449648|   4.750104|   4.431634|     4.392285|     4.497514|     4.762213|     4.533928|
|   6.090031|   5.973211|   5.913239|     6.238684|     6.050743|     5.932738|     6.022005|

