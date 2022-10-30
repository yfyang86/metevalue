#' Evalue of the methylKit data format
#'
#' Perform the Evaluation for the BiSeq data.
#' @param methyrate is the output file of methylKit, the columns are (in order):
#'     - chr: Chromosome
#'
#'     - pos: int Position
#'
#'     - g1~g2: methylation rate data in groups
#'
#' @param methylKit.output is the output data with e-value of each region
#'     - chr: Chromosome
#'
#'     - start: The positions of the start sites of the corresponding region
#'
#'     - end: The positions of the end sites of the corresponding region
#'
#'     - strand: Strand
#'
#'     - pvalue: The adjusted p-value based on BH method in MWU-test
#'
#'     - qvalue: cutoff for qvalue of differential methylation statistic
#'
#'     - methyl.diff: The difference between the group means of methylation level
#'
#' @param adjust.methods is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'
#' @param sep seperator, default is the TAB key.
#' @param bheader a logical value indicating whether the input_filename_b file contains the names of the variables as its first line. By default, bheader = FALSE.
#' @return a dataframe, the columns are (in order):
#'
#'     - chr: Chromosome
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
metevalue.methylKit <- function(methyrate, methylKit.output, adjust.methods='BH', sep = "\t", bheader = FALSE){
  re = metevalue.methylKit.chk(methyrate, methylKit.output, sep, bheader)
  return(varevalue.metilene(re$file_a, re$file_b, re$file_a_b, adjust.methods=adjust.methods));
}
