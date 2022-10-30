#' Evalue of the Metilene data
#'
#' Perform the Evaluation for the Metilene data. The data file could be pre-handled by the evalue.metilene.chk function.
#' @param a A data.frame object, the columns should be (in order):
#'
#' chrom	pos	g1	g1	g1	g1	g1	g1	g1	g1	g2	g2	g2	g2	g2	g2	g2	g2
#'
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
#' @param a_b A data.frame object of a join b with particular data clean processes. Check the function [evalue.methylKit.chk()] for more details.
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
#' @examples
#' #### methylKit example ####
#' data(demo_methylkit_methyrate)
#' data(demo_methylkit_met_all)
#' example_tempfiles = tempfile(c("rate_combine", "methylKit_DMR_raw"))
#' tempdir()
#' write.table(demo_methylkit_methyrate, file=example_tempfiles[1],
#'       row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
#' write.table (demo_methylkit_met_all, file=example_tempfiles[2],
#'       sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
#' result = metevalue.methylKit(example_tempfiles[1], example_tempfiles[2],
#'       bheader = TRUE)
#' str(result)
varevalue.metilene <- function(a, b, a_b, adjust.methods='BH'){
  innerf = Vectorize(function(x, innermu=0., innersig=1.){
    vector_temp = na.omit(x)
    n = length(vector_temp)
    value = mean(vector_temp)
    dnorm((value-innermu)*sqrt(n)/innersig)
  }, vectorize.args = 'x'
  )

  innerlog = function(x){
    log(x[!is.na(x)])
  }

  uid = unique(data.frame(start=a_b$start, end=a_b$end))
  g_count = (ncol(a) - 2) / 2

  site_1 <- 'g1'
  site_2 <- 'g2'

  if (g_count > 1){
    g_count = g_count - 1
    site_1 = c(site_1, paste("g1.", 1:g_count, sep = ''))
    site_2 = c(site_2, paste("g2.", 1:g_count, sep = ''))
  }

  withe = cbind(uid,e_value=0)

  for(i in 1:nrow(uid)){
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
  e_adjust= 1 / p.adjust(1/data_withe$e_value, method=adjust.methods)
  data_withe = cbind(data_withe, e_adjust)
  return(data_withe);
}

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

#' DMR BiSeq Demo Dataset
#'
#' The BiSeq dataset for demo purpose. The data are dummy data. It includes 9 columns:
#'
#'
#' - seqnames: Chromosome
#'
#'  - start: The positions of the start sites of the corresponding region
#'
#'  - end: The positions of the end sites of the corresponding region
#'
#'  - strand: Strand
#'
#'  - median.p
#'
#'  - median.meth.group1
#'
#'  - median.meth.group2
#'
#'  - median.meth.diff
#' @name demo_biseq_DMR
#' @docType data
#' @keywords metevalue
NULL

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
#' - seqnames
#'
#' - start
#'
#' - end
#'
#' - width
#'
#' - strand
#'
#' - median.p
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
#' - g1~g2: methylation rate data in groups, repeat 8 times.
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
