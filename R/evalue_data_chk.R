#' Buildin data process function
#' @param a filename
#' @param b filenmae
#' @param method "metilene" or "biseq" or "methylKit"
evalue_buildin_sql <- function(a, b, method="metilene"){
  result = NA
  if (method == "metilene"){
    result = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON (f.pos <= g.end AND f.pos>=g.start AND f.chr = g.chr)")
  }
  if (method == "biseq"){
    result = a_b = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON (f.pos <= g.end AND f.pos>=g.start AND f.chr = g.chr)")
  }
  if (method == "methylKit"){
    result = a_b = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON (f.pos <= g.end AND f.pos>=g.start AND f.chr = g.chr)")
  }
  if (method == "DMRfinder"){
    result = a_b = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON (f.pos <= g.end AND f.pos>=g.start AND f.chr = g.chr)")
  }
  return(result);
}

#' Buildin check file format function
#' @param a filename
#' @param b filenmae
#' @param method "metilene" or "biseq" or "methylKit"
evalue_buildin_var_fmt_nm <- function(a, b, method="metilene"){
  result = NA
  a = data.frame(a)
  b = data.frame(b)
  if (method == "metilene"){
    if ((ncol(a) - 2) %% 2 != 0) warning("Format error!!! ")
    site_1 <- 'g1'
    site_2 <- 'g2'
    g_count = (ncol(a) - 2) / 2
    if (g_count > 1){
      g_count = g_count - 1
      site_1 = c(site_1, paste("g1.", 1:g_count, sep = ''))
      site_2 = c(site_2, paste("g2.", 1:g_count, sep = ''))
    }
    names(a) <- c("chr", "pos", site_1, site_2)
    b_tab_names <- c("chr", "start", "end", "q-value",
                     "methyl.diff", "CpGs", "p", "p2",
                     "m1", "m2")
    names(b) = b_tab_names
  }
  if (method == "biseq"){
    if ((ncol(a) - 2) %% 2 != 0) warning("Format error!!! ")

    site_1 <- 'g1'
    site_2 <- 'g2'

    g_count = (ncol(a) - 2) / 2

    if (g_count > 1){
      g_count = g_count - 1
      site_1 = c(site_1, paste("g1.", 1:g_count, sep = ''))
      site_2 = c(site_2, paste("g2.", 1:g_count, sep = ''))
    }

    names(a) <- c("chr", "pos", site_1, site_2)


    b_tab_names <- c("chr","start","end","range","strand","median.p",
                     "median.meth.group1","median.meth.group2",
                     "median.meth.diff")
    names(b)[1:9] = b_tab_names
  }
  if (method == "methylKit"){
    if ((ncol(a) - 2) %% 2 != 0) warning("Format error!!! ")

    site_1 <- 'g1'
    site_2 <- 'g2'

    g_count = (ncol(a) - 2) / 2

    if (g_count > 1){
      g_count = g_count - 1
      site_1 = c(site_1, paste("g1.", 1:g_count, sep = ''))
      site_2 = c(site_2, paste("g2.", 1:g_count, sep = ''))
    }

    names(a) <- c("chr", "pos", site_1, site_2)

    b_tab_names <- c("chr", "start", "end", "strand", "p", "qvalue", "meth.diff")
    names(b) = b_tab_names
  }

  if (method == "DMRfinder"){
    if ((ncol(a) - 2) %% 2 != 0) warning("Format error!!! ")

    site_1 <- 'g1'
    site_2 <- 'g2'

    g_count = (ncol(a) - 2) / 2

    if (g_count > 1){
      g_count = g_count - 1
      site_1 = c(site_1, paste("g1.", 1:g_count, sep = ''))
      site_2 = c(site_2, paste("g2.", 1:g_count, sep = ''))
    }

    names(a) <- c("chr", "pos", site_1, site_2)
    b_tab_names <- c("chr", "start", "end", "CpG", "Control_mu",
    "Exptl_mu",	"Control_Exptl_diff", "p")
    names(b) = b_tab_names
  }

  return(result=list(a = a, b = b));
}



#' Check the Metilene data format
#' @param input_filename_a metilene input file path. This file is a sep (e.g. TAB) separated file with two key columns and several value columns in pairs:
#' For exampe:
#'
#' chrom	pos	g1	g1	g1	g1	g1	g1	g1	g1	g2	g2	g2	g2	g2	g2	g2	g2
#'
#' chrom and pos are keys;
#' g1 g1 g2 g2 must be stored in pairs.

#' @param input_filename_b  metilene input file path. This file should stored as a sep(e.g. TAB) separated file with two key columns and several value columns:
#' The columns are (in order):
#'
#'     - chr:   Chromosome
#'
#'     - start: The position of the start sites of the corresponding region
#'
#'     - end: The position of the end sites of the corresponding region
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
#' @param sep separator, default is the TAB key.
#' @param bheader a logical value indicating whether the input_filename_b file contains the names of the variables as its first line. By default, bheader = FALSE.
#' @return list(file_a, file_b, file_a_b) returns a list with three pr-handled data.frames corresponding to the input_filename_a, input_filename_b file and a A JOIN B file.
#'
#'@examples
#'\dontrun{
#' s = evalue.metilene(input_filename_a = "../metilene/metilene_8.input",
#' input_filename_b = "../metilene/metilene_8.out")
#' ## > str(s)
#' ## data.frame':	723 obs. of  11 variables:
#' ## $ chr        : chr  "chr21" "chr21" "chr21" "chr21" ...
#' ## $ start      : int  9437432 9708982 9825467 9825514 9825794 9825882 9826220 9827335 9926563 9963396 ...
#' ## $ end        : int  9437540 9709189 9825508 9825788 9825876 9826191 9826387 9827356 9927097 9963753 ...
#' ## $ q-value    : num  2.49e-25 4.62e-29 6.00e-02 3.40e-01 2.82e-07 ...
#' ## $ methyl.diff: num  0.611 0.476 -0.274 -0.164 -0.261 ...
#' ## $ CpGs       : int  26 28 12 26 26 31 13 10 73 10 ...
#' ## $ p         : num  3.86e-14 4.34e-14 2.60e-07 2.55e-05 1.23e-11 ...
#' ## $ p2        : num  2.23e-29 4.13e-33 5.37e-06 3.04e-05 2.52e-11 ...
#' ## $ m1         : num  0.737 0.589 0.298 0.374 0.353 ...
#' ## $ m2         : num  0.126 0.113 0.573 0.538 0.615 ...
#' ## $ e_value    : num  8.50e+40 2.71e+38 7.20e+05 1.89e+06 5.36e+12 ...
#'}
evalue.metilene.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){

    a <- read.table(input_filename_a, header=T, sep=sep)
    b <- read.table(input_filename_b, header=bheader, sep=sep)

    re <- evalue_buildin_var_fmt_nm(a, b, method = "metilene")


return(list(file_a = re$a, file_b = re$b, file_a_b = evalue_buildin_sql(re$a, re$b, method = "metilene")));
}

#' Check the BiSeq data format
#' @param input_filename_a metilene input file path. This file is a sep (e.g. TAB) separated file with two key columns and several value columns in pairs:
#' For exampe:
#'
#' chrom	pos	g1	g1	g1	g1	g1	g1	g1	g1	g2	g2	g2	g2	g2	g2	g2	g2
#'
#' chrom and pos are keys;
#' g1 g1 g2 g2 must be stored in pairs.

#' @param input_filename_b  metilene input file path. This file should stored as a sep(e.g. TAB) separated file with two key columns and several value columns:
#' The columns are (in order):
#'
#'     - chr:   Chromosome
#'
#'     - start: The position of the start site of the corresponding region
#'
#'     - end: The position of the end site of the corresponding region
#'
#'     - range: The range of the corresponding region
#'
#'     - strand: Strand
#'
#'     - median.p:  The median of p-values in the corresponding region
#'
#'     - median.meth.group1 : The median of methylation level for the corresponding segment of group 1
#'
#'     - median.meth.group2 : The median of methylation level for the corresponding segment of group 2
#'
#'     - median.meth.diff:  The median of the difference between the methylation level
#'
#'
#' @param sep separator, default is the TAB key.
#' @param bheader a logical value indicating whether the input_filename_b file contains the names of the variables as its first line. By default, bheader = FALSE.
#' @return list(file_a, file_b, file_a_b) returns a list with three pr-handled data.frames corresponding to the input_filename_a, input_filename_b file and a A JOIN B file.
#' @examples
#' \dontrun{
#' s = evalue.metilene(input_filename_a = "../metilene/metilene_8.input",
#' input_filename_b = "../metilene/metilene_8.out")
#' ## > str(s)
#' ## data.frame':	723 obs. of  11 variables:
#' ## $ chr        : chr  "chr21" "chr21" "chr21" "chr21" ...
#' ## $ start      : int  9437432 9708982 9825467 9825514 9825794 9825882 9826220 9827335 9926563 9963396 ...
#' ## $ end        : int  9437540 9709189 9825508 9825788 9825876 9826191 9826387 9827356 9927097 9963753 ...
#' ## $ q-value    : num  2.49e-25 4.62e-29 6.00e-02 3.40e-01 2.82e-07 ...
#' ## $ methyl.diff: num  0.611 0.476 -0.274 -0.164 -0.261 ...
#' ## $ CpGs       : int  26 28 12 26 26 31 13 10 73 10 ...
#' ## $ p         : num  3.86e-14 4.34e-14 2.60e-07 2.55e-05 1.23e-11 ...
#' ## $ p2        : num  2.23e-29 4.13e-33 5.37e-06 3.04e-05 2.52e-11 ...
#' ## $ m1         : num  0.737 0.589 0.298 0.374 0.353 ...
#' ## $ m2         : num  0.126 0.113 0.573 0.538 0.615 ...
#' ## $ e_value    : num  8.50e+40 2.71e+38 7.20e+05 1.89e+06 5.36e+12 ...
#' }
evalue.biseq.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){
  a <- read.table(input_filename_a, header=T, sep=sep)
  b <- read.table(input_filename_b, header=bheader, sep=sep)
  re <- evalue_buildin_var_fmt_nm(a, b, method = "biseq")
  return(list(file_a = re$a, file_b = re$b, file_a_b = evalue_buildin_sql(re$a, re$b, method = "biseq")));
}

#' Check the methylKit data format
#' @param input_filename_a the combined data of methylation rate file. This file is a sep (e.g. TAB) separated file with two key columns and several value columns in pairs:
#' For exampe:
#'
#' chrom	pos	g1	g1	g1	g1	g1	g1	g1	g1	g2	g2	g2	g2	g2	g2	g2	g2
#'
#' chrom and pos are keys;
#' g1 g1 g2 g2 must be stored in pairs.

#' @param input_filename_b  the output file of methylKit. a methylDiff or methylDiffDB object containing the differential methylated locations satisfying the criteria.
#' The columns are (in order):
#'
#'     - chr:   Chromosome
#'
#'     - start: The position of the start sites of the corresponding region
#'
#'     - end: The position of the end sites of the corresponding region
#'
#'     - strand: Strand
#'
#'     - p: p-value
#'
#'     - qvalue:  The adjusted p-value based on BH method
#'
#'     - meth.diff : The difference between the group means of methylation level
#'
#' @param sep separator, default is the TAB key.
#' @param bheader a logical value indicating whether the input_filename_b file contains the names of the variables as its first line. By default, bheader = FALSE.
#' @return list(file_a, file_b, file_a_b) returns a list with three pr-handled data.frames corresponding to the input_filename_a, input_filename_b file and a A JOIN B file.
#' @examples
#' \dontrun{
#' s = evalue.metilene(input_filename_a = "../metilene/metilene_8.input",
#' input_filename_b = "../metilene/metilene_8.out")
#' ## > str(s)
#' ## data.frame':	723 obs. of  11 variables:
#' ## $ chr        : chr  "chr21" "chr21" "chr21" "chr21" ...
#' ## $ start      : int  9437432 9708982 9825467 9825514 9825794 9825882 9826220 9827335 9926563 9963396 ...
#' ## $ end        : int  9437540 9709189 9825508 9825788 9825876 9826191 9826387 9827356 9927097 9963753 ...
#' ## $ q-value    : num  2.49e-25 4.62e-29 6.00e-02 3.40e-01 2.82e-07 ...
#' ## $ methyl.diff: num  0.611 0.476 -0.274 -0.164 -0.261 ...
#' ## $ CpGs       : int  26 28 12 26 26 31 13 10 73 10 ...
#' ## $ p         : num  3.86e-14 4.34e-14 2.60e-07 2.55e-05 1.23e-11 ...
#' ## $ p2        : num  2.23e-29 4.13e-33 5.37e-06 3.04e-05 2.52e-11 ...
#' ## $ m1         : num  0.737 0.589 0.298 0.374 0.353 ...
#' ## $ m2         : num  0.126 0.113 0.573 0.538 0.615 ...
#' ## $ e_value    : num  8.50e+40 2.71e+38 7.20e+05 1.89e+06 5.36e+12 ...
#' }
evalue.methylKit.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){
  a <- read.table(input_filename_a, header=T, sep=sep)
  b <- read.table(input_filename_b, header=bheader, sep=sep)
  re <- evalue_buildin_var_fmt_nm(a, b, method = "methylKit")
  return(list(file_a = re$a, file_b = re$b, file_a_b = evalue_buildin_sql(re$a, re$b, method = "methylKit")));
}


#' Check the DMRfinder data format
#' @param input_filename_a the combined data of methylation rate file. This file is a sep (e.g. TAB) separated file with two key columns and several value columns in pairs:
#' For exampe:
#'
#' chrom	pos	g1	g1	g1	g1	g1	g1	g1	g1	g2	g2	g2	g2	g2	g2	g2	g2
#'
#' chrom and pos are keys;
#' g1 g1 g2 g2 must be stored in pairs.

#' @param input_filename_b  the output file of DMRfinder. xxxxxxxxxxxxxx
#' The columns are (in order):
#'
#'     - chr:   Chromosome
#'
#'     - start: The position of the start sites of the corresponding region
#'
#'     - end: The position of the end sites of the corresponding region
#'
#'     - CpG: The number of CpG sites within the corresponding region
#'
#'     - `Control:mu`: The absolute mean methylation level for the corresponding segment of the control group
#'
#'     - `Exptl:mu`: The absolute mean methylation level for the corresponding segment of the experimental group
#'
#'     - `Control->Exptl:diff`: The difference between the group means of methylation level
#'
#'     - p: p-value
#'
#' @param sep separator, default is the TAB key.
#' @param bheader a logical value indicating whether the input_filename_b file contains the names of the variables as its first line. By default, bheader = FALSE.
#' @return list(file_a, file_b, file_a_b) returns a list with three pr-handled data.frames corresponding to the input_filename_a, input_filename_b file and a A JOIN B file.
#' @examples
#' \dontrun{
#' s = evalue.metilene(input_filename_a = "../metilene/metilene_8.input",
#' input_filename_b = "../metilene/metilene_8.out")
#' ## > str(s)
#' ## data.frame':	723 obs. of  11 variables:
#' ##   $ chr        : chr  "chr21" "chr21" "chr21" "chr21" ...
#' ## $ start      : int  9437432 9708982 9825467 9825514 9825794 9825882 9826220 9827335 9926563 9963396 ...
#' ## $ end        : int  9437540 9709189 9825508 9825788 9825876 9826191 9826387 9827356 9927097 9963753 ...
#' ## $ q-value    : num  2.49e-25 4.62e-29 6.00e-02 3.40e-01 2.82e-07 ...
#' ## $ methyl.diff: num  0.611 0.476 -0.274 -0.164 -0.261 ...
#' ## $ CpGs       : int  26 28 12 26 26 31 13 10 73 10 ...
#' ## $ p         : num  3.86e-14 4.34e-14 2.60e-07 2.55e-05 1.23e-11 ...
#' ## $ p2        : num  2.23e-29 4.13e-33 5.37e-06 3.04e-05 2.52e-11 ...
#' ## $ m1         : num  0.737 0.589 0.298 0.374 0.353 ...
#' ## $ m2         : num  0.126 0.113 0.573 0.538 0.615 ...
#' ## $ e_value    : num  8.50e+40 2.71e+38 7.20e+05 1.89e+06 5.36e+12 ...
#' }
evalue.DMRfinder.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){
  a <- read.table(input_filename_a, header=T, sep=sep)
  b <- read.table(input_filename_b, header=bheader, sep=sep)
  re <- evalue_buildin_var_fmt_nm(a, b, method = "DMRfinder")
  return(list(file_a = re$a, file_b = re$b, file_a_b = evalue_buildin_sql(re$a, re$b, method = "DMRfinder")));
}