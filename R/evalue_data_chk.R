#' Build-in data process function
#' @param a data frame of the methylation rate
#' @param b data frame of output data corresponding to the "method" option
#' @param method "metilene" or "biseq", "DMRfinder" or "methylKit"
#' @return a data frame combines data frame a and b corresponding to the "method" option
#'
#' @examples
#' \donttest{
#' data("demo_metilene_out")
#' data("demo_metilene_input")
#' result = evalue_buildin_var_fmt_nm(demo_metilene_input,
#'                                    demo_metilene_out, method="metilene")
#' result_sql = evalue_buildin_sql(result$a, result$b, method="metilene")
#' }
evalue_buildin_sql <- function(a, b, method="metilene"){
  result = NA
  if (method == "metilene"){
    result = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON f.chr = g.chr where f.pos<=g.end AND f.pos>=g.start")
  }
  if (method == "biseq"){
    result = a_b = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON f.chr = g.chr where f.pos<=g.end AND f.pos>=g.start")
  }
  if (method == "methylKit"){
    result = a_b = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON f.chr = g.chr where f.pos<=g.end AND f.pos>=g.start")
  }
  if (method == "DMRfinder"){
    result = a_b = sqldf("SELECT * FROM b AS g
                 LEFT JOIN a AS f
                 ON f.chr = g.chr where f.pos<=g.end AND f.pos>=g.start")
  }
  return(result);
}

#' Build-in check file format function
#' Perform the format check and data clean for the  "metilene" or "biseq", "DMRfinder" or "methylKit" method correspondingly.
#' @param a data frame of the methylation rate
#' @param b data frame of output data corresponding to the "method" option
#' @param method "metilene" or "biseq", "DMRfinder" or "methylKit"
#' @return list(a, b) which contains the cleaned data correspondingly
#'
#' @examples
#' \donttest{
#' data("demo_metilene_out")
#' data("demo_metilene_input")
#' evalue_buildin_var_fmt_nm(demo_metilene_input,
#'                           demo_metilene_out, method="metilene")
#' }
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
    if(ncol(b) < length(b_tab_names)){
      stop("File Column Mismatch. Please Check the file format.")
    }
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
    if(ncol(b) < length(b_tab_names)){
      stop("File Column Mismatch. Please Check the file format.")
    }
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
    if(ncol(b) < length(b_tab_names)){
      stop("File Column Mismatch. Please Check the file format.")
    }
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
    if(ncol(b) < length(b_tab_names)){
      stop("File Column Mismatch. Please Check the file format.")
    }
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
#' @examples
#' data("demo_metilene_out")
#' data("demo_metilene_input")
metevalue.metilene.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){

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
#' data("demo_biseq_DMR")
#' data("demo_biseq_methyrate")
metevalue.biseq.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){
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
#' #### methylKit example ####
#' data(demo_methylkit_methyrate)
#' data(demo_methylkit_met_all)
metevalue.methylKit.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){
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

#' @param input_filename_b  the output file of DMRfinder.
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
#'
#' @examples
#' data("demo_DMRfinder_rate_combine")
#' data("demo_DMRfinder_DMRs")
metevalue.DMRfinder.chk <- function(input_filename_a, input_filename_b, sep = "\t", bheader = FALSE){
  a <- read.table(input_filename_a, header=T, sep=sep)
  b <- read.table(input_filename_b, header=bheader, sep=sep)
  re <- evalue_buildin_var_fmt_nm(a, b, method = "DMRfinder")
  return(list(file_a = re$a, file_b = re$b, file_a_b = evalue_buildin_sql(re$a, re$b, method = "DMRfinder")));
}
