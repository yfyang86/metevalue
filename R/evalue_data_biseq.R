#' Evalue of the BiSeq data format
#'
#' Perform the Evaluation for the BiSeq data. Please check vignette "metevalue" for details.
#' @param methyrate is the methyrate file. The columns are (in order):
#'     - chr: Chromosome
#'
#'     - pos: int Position
#'
#'     - g1~g2: methylation rate data in groups
#' @param BiSeq.output is the output file of BiSeq. The columns are (in order):
#'     - seqnames: Chromosome
#'
#'     - start: The positions of the start sites of the corresponding region
#'
#'     - end: The positions of the end sites of the corresponding region
#'
#'     - width: The number of CpG sites within the corresponding region
#'
#'     - strand: Strand
#'
#'     - median.p: The median p-value among CpG sites within the corresponding region
#'
#'     - median.meth.group1: The median methylation rate in the first group among CpG sites within the corresponding region
#'
#'     - median.meth.group2: The median methylation rate in the second group among CpG sites within the corresponding region
#'
#'     - median.meth.diff: The median methylation difference between groups among CpG sites within the corresponding region
#' @param adjust.methods is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'
#' @param sep seperator, default is the TAB key.
#' @param bheader a logical value indicating whether the BiSeq.output file contains the names of the variables as its first line. By default, bheader = FALSE.
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
#' \donttest{
#' data("demo_biseq_methyrate")
#' data("demo_biseq_DMR")
#' example_tempfiles = tempfile(c("demo_biseq_methyrate", "demo_biseq_DMR"))
#' tempdir()
#' #### write to temp file ####
#' write.table(demo_biseq_methyrate, file=example_tempfiles[1],row.names=FALSE,
#'             col.names=TRUE, quote=FALSE, sep='\t')
#' write.table (demo_biseq_DMR, file=example_tempfiles[2],
#'              sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
#' #### compute e-value and its adjustment ####
#' result = metevalue.biseq(example_tempfiles[1],
#'                          example_tempfiles[2], bheader = TRUE)
#' }
metevalue.biseq <- function(methyrate, BiSeq.output, adjust.methods='BH', sep = "\t", bheader = FALSE){
  re = metevalue.biseq.chk (methyrate, BiSeq.output, sep, bheader)
  return(varevalue.metilene(re$file_a, re$file_b, re$file_a_b, adjust.methods=adjust.methods));
}
