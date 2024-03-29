#' Calculate E-value of the Metilene data
#'
#' The data file could be pre-handled by the evalue.metilene.chk function.
#' @param a A data.frame object:
#' \tabular{rrrrrrrr}{
#'  chr	\tab  pos	 \tab   g1	\tab ...  \tab  g1 \tab  g2 \tab ... \tab g2 \cr
#' chr1 \tab  1    \tab  0.1 \tab  ... \tab   0.1\tab  0.2\tab ... \tab 0.2\cr
#' }
#' i.e two key columns (chrom, pos) with several value columns in groups.
#' @param b A data.frame object stores the data, the columns are (in order):
#'
#'     - chr:   Chromosome
#'
#'     - start: The positions of the start sites of the corresponding region
#'
#'     - end: The positions of the end sites of the corresponding region
#'
#'     - q-value: The adjusted p-value based on BH method in MWU-test
#'
#'     - methyl.diff: The difference between the group means of methylation level
#'
#'     - CpGs:  The number of CpG sites within the corresponding region
#'
#'     - p : p-value based on MWU-test
#'
#'     - p2: p-value based on 2D KS-test
#'
#'     - m1:  The absolute mean methylation level for the corresponding segment of group 1
#'
#'     - m2:  The absolute mean methylation level for the corresponding segment of group 2
#'
#' @param a_b A data.frame object of a join b with particular data clean processes. Check the function [evalue.methylKit.chk()] for more details.
#' @param group1_name charactor: The name of the first group. For example, "g1" in the above example.
#' @param group2_name charactor: The name of the second group. For example, "g2" in the above example.
#' @param adjust.methods is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'. The default value is 'BH'.
#' @return a dataframe, the columns are (in order):
#'
#'     - chr:   Chromosome
#'
#'     - start: The positions of the start sites of the corresponding region
#'
#'     - end: The positions of the end sites of the corresponding region
#'
#'     - q-value: The adjusted p-value based on BH method in MWU-test
#'
#'     - methyl.diff: The difference between the group means of methylation level
#'
#'     - CpGs:  The number of CpG sites within the corresponding region
#'
#'     - p : p-value based on MWU-test
#'
#'     - p2: p-value based on 2D KS-test
#'
#'     - m1:  The absolute mean methylation level for the corresponding segment of group 1
#'
#'     - m2:  The absolute mean methylation level for the corresponding segment of group 2
#'
#'     - e_value: The e-value of the corresponding region
#'
#' @examples
#' #data(demo_metilene_input)
#' #data(demo_metilene_out)
#' #result = evalue_buildin_var_fmt_nm(demo_metilene_input, demo_metilene_out, method="metilene")
#' #result = list(a = result$a, 
#' #              b = result$b, 
#' #              a_b = evalue_buildin_sql(result$a, result$b, method="metilene"))
#' #result = varevalue.metilene(result$a, result$b, result$a_b)
varevalue.metilene <- function(a, b, a_b, group1_name = 'g1', group2_name = 'g2', adjust.methods='BH'){
  innerf = function(x, innermu=0., innersig=1.){
    vector_temp = na.omit(as.numeric(x))
    n = length(vector_temp)
    value = mean(vector_temp)
    dnorm(x=value, mean = innermu, sd = innersig/sqrt(n))
  }

  innerlog = function(x){
    log(x[!is.na(x)])
  }

  uid = unique(data.frame(start=a_b$start, end=a_b$end))
  g_count = (ncol(a) - 2) / 2

  site_1 = grep(paste0('^',group1_name), names(a_b), value=T)
  site_2 = grep(paste0('^',group2_name), names(a_b), value=T)
  
  withe = cbind(uid,e_value=0)

  for(i in seq_len(nrow(uid))){
    uid_temp = uid[i,]
    start_temp = uid_temp[1,1]
    end_temp = uid_temp[1,2]
    a_b_temp = a_b[(a_b$start==start_temp & a_b$end==end_temp), ]
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
    withe[i,3] = e_value
  }
  data_withe = left_join(b, withe)
  e_adjust= 1 / p.adjust(1/(data_withe$e_value), method=adjust.methods)
  data_withe = cbind(data_withe, e_adjust)
  return(data_withe);
}



#' A general method to calculate the e-value for other DNA methylation tools not described above. The input data is the DNA methylation rates using the similar format with Metilene.
#'
#' The input data file is just the DNA methylation rates using the similar format above, with no need for another data file output by different tools.
#' The Chromosome name, start and end sites shoule be specified in the function.
#' @param methyrate data.frame: A data.frame object of methylation rates, the columns should be (name of groups can be self-defined)
#' \tabular{rrrrrrrr}{
#'  chr	\tab  pos	 \tab   g1	\tab ...  \tab  g1 \tab  g2 \tab ... \tab g2 \cr
#' chr1 \tab  1    \tab  0.1 \tab  ... \tab   0.1\tab  0.2\tab ... \tab 0.2\cr
#' }
#' @param group1_name charactor: The name (pattern) of the first group. For example, "g1" in the above example.
#' For example `g1_abc` and `g1` will be considered as the same group if `group1_name = "g1"`.  Use this with care in practice.
#' @param group2_name charactor: The name (pattern) of the second group. For example, "g2" in the above example.
#' For example `g2_abc` and `g2` will be considered as the same group if `group2_name = "g2"`.  Use this with care in practice.
#' @param chr charactor: The Chromosome name. Typically, it is a string like "chr21" and so on.
#' @param start integer:  The position of the start site of the corresponding region
#' @param end integer: The position of the end site of the corresponding region
#' @return evalue
#' @examples
#' #data("demo_metilene_input")
#' #varevalue.single_general(demo_metilene_input, chr = "chr21", start = 9437432, end = 9437540)
#' # [1] 2.626126e+43
#' 
#' #### Compare to `metevalue.metilene`  ####
#' data(demo_metilene_out)
#' #example_tempfiles = tempfile(c("metilene_input", "metilene_out"))
#' #tempdir()
#' #write.table(demo_metilene_input, file=example_tempfiles[1],
#' #      row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
#' #write.table (demo_metilene_out, file=example_tempfiles[2],
#' #      sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
#' #result = metevalue.metilene(example_tempfiles[1], example_tempfiles[2],
#' #      bheader = TRUE)
#' # result[with(result, chr == 'chr21' & start == '9437432' & end == '9437540'), ncol(result)]
#' # [1] 2.626126e+43
varevalue.single_general = function(methyrate, group1_name='g1', group2_name='g2', chr, start, end){
  innerf = function(x, innermu=0., innersig=1.){
    vector_temp = na.omit(as.numeric(x))
    n = length(vector_temp)
    value = mean(vector_temp)
    dnorm(x=value, mean = innermu, sd = innersig/sqrt(n))
  }
  
  innerlog = function(x){
    log(x[!is.na(x)])
  }
  a_b = methyrate
  site_1 = grep(paste0('^',group1_name), names(a_b), value=T)
  site_2 = grep(paste0('^',group2_name), names(a_b), value=T)
  
  start_temp = start
  end_temp = end
  chr_temp = chr
  a_b_temp = a_b[(a_b$pos>=start_temp & a_b$pos<=end_temp & a_b$'chr'==chr_temp), ]
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
  # print(chr_temp)
  return(e_value)
}

#' A general method to calculate the e-value for RNA-seq data.
#'
#' @param rna data.frame: A data.frame object of RNAseq data. For example:
#' 
#' \tabular{rrrrr}{
#' TAG   \tab  treated1fb \tab treated2fb  \tab untreated1fb \tab untreated2fb \cr
#' TAG1  \tab 4.449648    \tab 4.750104    \tab   4.392285   \tab   4.497514 \cr
#' TAG2  \tab 8.241116    \tab 8.302852    \tab  8.318125    \tab  8.488796 \cr
#' ...   \tab ...         \tab ...         \tab  ...         \tab  ...      \cr
#' }
#' 
#' Row names (TAG1 and TAG2 in the above example) is also suggested.
#' @param group1_name charactor: The name (pattern) of the first group. For example, "treated" in the above example.
#' For example `treated_abc` and `treated` will be considered as the same group if `group1_name = "treated"`.  Use this with care in practice.
#' @param group2_name charactor: The name (pattern) of the second group. For example, "untreated" in the above example.
#' For example `untreated_abc` and `untreated` will be considered as the same group if `group2_name = "untreated"`.  Use this with care in practice.
#' @return evalue 
#' @examples
#' data("demo_desq_out")
#' evalue = metevalue.RNA_general(demo_desq_out, 'treated','untreated')
metevalue.RNA_general = function(rna, group1_name, group2_name){
  a_b = data.frame(rna)
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



#' DESeq Output Dataset
#'
#' The output dummy data for "RNA" meythod illustrating purpose.
#'
#' The data includes 10 columns.
#'
#' - treated1fb:
#' 
#' - treated2fb:
#' 
#' - treated3fb:
#' 
#' - untreated1fb:
#' 
#' - untreated2fb:
#' 
#' - untreated3fb:
#'  
#' - untreated4fb:
#' 
#' This data contains 8166 rows and 7 columns.
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_desq_out
#' @docType data
#' @keywords metevalue
#' The data is a simulation data:
#' @examples
#' # library("pasilla")
#' # pasCts <- system.file("extdata",
#' #                       "pasilla_gene_counts.tsv",
#' #                       package="pasilla", mustWork=TRUE)
#' # pasAnno <- system.file("extdata",
#' #                        "pasilla_sample_annotation.csv",
#' #                        package="pasilla", mustWork=TRUE)
#' # cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
#' # coldata <- read.csv(pasAnno, row.names=1)
#' # coldata <- coldata[,c("condition","type")]
#' # coldata$condition <- factor(coldata$condition)
#' # coldata$type <- factor(coldata$type)
#' # 
#' # library("DESeq2")
#' # colnames(cts)=paste0(colnames(cts),'fb')
#' # cts = cts[,rownames(coldata)]
#' # dds <- DESeqDataSetFromMatrix(countData = cts,
#' #                               colData = coldata,
#' #                               design = ~ condition)
#' # dds <- DESeq(dds)
#' # 
#' # 
#' # dat <- t(t(cts)/(dds$sizeFactor)) 
#' # dat.out <- dat[rowSums(dat >5)>=0.8*ncol(dat),]
#' # 
#' # demo_desq_out <- log(dat.out)
NULL


#' Methyrate Dataset
#'
#' The methyrate dataset samples "myCpG" data from the methylKit (a bioconductor package) for illustrating purpose.
#'
#' The data includes 6 columns.
#'
#' - chr: string Chromosome
#'
#' - pos: int Position
#'
#' - g1~g2: methylation rate data in groups (4 columns)
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_methylkit_methyrate
#' @docType data
#' @references Akalin, Altuna, et al. "methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles." Genome biology 13.10 (2012): 1-9. \doi{10.1186/gb-2012-13-10-r87}
#' @keywords metevalue
NULL

#' Methyrate output dataset from methylKit
#'
#' The methyrate dataset samples "myCpG" data from the methylKit (a bioconductor package) for illustrating purpose.
#'
#' The data includes 7 columns:
#'
#'  - chr:   Chromosome
#'
#'  - start: The positions of the start sites of the corresponding region
#'
#'  - end: The positions of the end sites of the corresponding region
#'
#'  - strand: Strand
#'
#'  - pvalue: The adjusted p-value based on BH method in MWU-test
#'
#'  - qvalue: cutoff for qvalue of differential methylation statistic
#'
#'  - methyl.diff: The difference between the group means of methylation level
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_methylkit_met_all
#' @docType data
#' @references Akalin, Altuna, et al. "methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles." Genome biology 13.10 (2012): 1-9. \doi{10.1186/gb-2012-13-10-r87}
#' @keywords metevalue
NULL

# #' DMR BiSeq Demo Dataset
# #'
# #' The BiSeq dataset for demo purpose. The data are dummy data. It includes 8 columns:
# #'
# #'
# #' - seqnames: Chromosome
# #'
# #'  - start: The positions of the start sites of the corresponding region
# #'
# #'  - end: The positions of the end sites of the corresponding region
# #'
# #'  - strand: Strand
# #'
# #'  - median.p
# #'
# #'  - median.meth.group1
# #'
# #'  - median.meth.group2
# #'
# #'  - median.meth.diff
# #'
# #' @name demo_biseq_DMR
# #' @docType data
# #' @keywords metevalue
# #NULL

#' BiSeq Methyrate Demo Dataset
#'
#' The methyrate for BiSeq illustrating purpose. It is dummy.
#'
#' The data includes 12 columns.
#'
#' - chr: string Chromosome
#'
#' - pos: int Position
#'
#' - g1~g2: methylation rate data in groups, repeat 5 times.
#' Notice that there are "NaN" within the feature columns.
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_biseq_methyrate
#' @docType data
#' @keywords metevalue
NULL

#' BiSeq Output Demo Dataset
#'
#' The dummy output for BiSeq illustrating purpose. It is dummy.
#'
#' - seqnames: Chromosome
#'
#'  - start: The positions of the start sites of the corresponding region
#'
#'  - end: The positions of the end sites of the corresponding region
#'
#' - width
#'
#'  - strand: Strand
#'
#'  - median.p
#'
#' - median.meth.group1
#'
#' - median.meth.group2
#'
#' - median.meth.diff
#'
#' Notice that there are "NaN" within the feature columns.
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_biseq_DMR
#' @docType data
#' @keywords metevalue
NULL

#' DMRfinder Methyrate Demo Dataset
#'
#' The methyrate for BiSeq illustrating purpose. It is dummy.
#'
#' The data includes 6 columns.
#'
#' - chr: string Chromosome
#'
#' - pos: int Position
#'
#' - g1~g2: methylation rate data in groups, repeat 2 times.
#' Notice that there are "NaN" within the feature columns.
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_DMRfinder_rate_combine
#' @docType data
#' @keywords metevalue
NULL

#' DMRfinder Output Demo Dataset
#'
#' The output dummy dataset for DMRfinder illustrating purpose.
#'
#' The data includes 6 columns.
#'
#' - chr: string Chromosome
#'
#' - pos: int Position
#'
#' - g1~g2: methylation rate data in groups, repeat 2 times.
#' Notice that there are "NaN" within the feature columns.
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_DMRfinder_DMRs
#' @docType data
#' @keywords metevalue
NULL

#' Metilene Methyrate Demo Dataset
#'
#' The methyrate for metilene illustrating purpose. It is dummy.
#'
#' The data includes 18 columns.
#'
#' - chr: string Chromosome
#'
#' - pos: int Position
#'
#' - g1~g2: methylation rate data in groups.
#'
#' Notice that there are "NaN" within the feature columns.
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_metilene_input
#' @docType data
#' @keywords metevalue
NULL

#' Metilene Demo Output Dataset
#'
#' The output dummy data for "metilene" meythod illustrating purpose.
#'
#' The data includes 10 columns.
#'
#' - V1: string Chromosome
#'
#' - V2: The positions of the start sites of the corresponding region
#'
#' - V3: The positions of the end sites of the corresponding region
#'
#' - V4- V10: data value.
#'
#' Please check the vignette "metevalue" for details.
#' @name demo_metilene_out
#' @docType data
#' @keywords metevalue
NULL

## #' General Purpose E-value for Metilene Data with 3 dimensions
## #' @param value vector (dimension must be 3)
## #' @param level integer
## #' @return A numeric vector of corrected e-values
## #' @export
## metevalue.adjust = function(value, level = 1){
##   if (length(dim(value)) != 3){
##     stop("Dimensions of value should be 3!");
##   }
##   rijbar = apply(value, c(1,2), mean, na.rm = TRUE)
##   mui = apply(value, 1, mean, na.rm = TRUE)
##   sigmai = apply(value, 1, function(x) var(as.vector(x), na.rm = TRUE))
##   mu_overall = mean(as.vector(value), na.rm = TRUE)
##   sigma_overall = var(as.vector(value))
## }
