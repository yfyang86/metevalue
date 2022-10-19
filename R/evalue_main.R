#' Evalue of the Metilene data
#'
#' Perform the Evaluation for the Metilene data. The data file could be pre-handled by the evalue.metilene.chk function.
#' @param a A data.frame object, the columns should be (in order):
#'
#' chrom	pos	g1	g1	g1	g1	g1	g1	g1	g1	g2	g2	g2	g2	g2	g2	g2	g2
#'
#' i.e two key columns (chrom, pos) with several value columns in pairs.
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
#'@examples
#'\dontrun{
#' s = evalue.metilene(input_filename_a = "metilene.input",
#' input_filename_b = "metilene.out")
#' ## > str(s)
#' ## data.frame':	723 obs. of  11 variables:
#' ## $ chr        : chr  "chr21" "chr21" "chr21" "chr21" ...
#' ## $ start      : int  9437432 9708982 9825467 ...
#' ## $ end        : int  9437540 9709189 9825508 ...
#' ## $ q-value    : num  2.49e-25 4.62e-29 6.00e-02 3.40e-01 2.82e-07 ...
#' ## $ methyl.diff: num  0.611 0.476 -0.274 -0.164 -0.261 ...
#' ## $ CpGs       : int  26 28 12 26 26 31 13 10 73 10 ...
#' ## $ p         : num  3.86e-14 4.34e-14 2.60e-07 2.55e-05 1.23e-11 ...
#' ## $ p2        : num  2.23e-29 4.13e-33 5.37e-06 3.04e-05 2.52e-11 ...
#' ## $ m1         : num  0.737 0.589 0.298 0.374 0.353 ...
#' ## $ m2         : num  0.126 0.113 0.573 0.538 0.615 ...
#' ## $ e_value    : num  8.50e+40 2.71e+38 7.20e+05 1.89e+06 5.36e+12 ...
#'}
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
