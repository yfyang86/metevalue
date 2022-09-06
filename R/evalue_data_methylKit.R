#' Evalue of the methylKit data format
#'
#' Perform the Evaluation for the BiSeq data.
#' @param methylKit.output' is the output file of methylKit
#' @param outputfile' is the output data with e-value of each region
#' @param adjust.methods' is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'
#' @param sep seperator, default is the TAB key.
#' @param bheader a logical value indicating whether the input_filename_b file contains the names of the variables as its first line. By default, bheader = FALSE.
evalue.methylKit <- function(methyrate, methylKit.output, adjust.methods='BH', sep = "\t", bheader = FALSE){
  re = evalue.methylKit.chk(methyrate, methylKit.output, sep, bheader)
  return(varevalue.metilene(re$file_a, re$file_b, re$file_a_b, adjust.methods=adjust.methods));
}
