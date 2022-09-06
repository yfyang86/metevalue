#' Evalue of the DMRfinder data format
#'
#' Perform the Evaluation for the BiSeq data.
#' @param methyrate is the methyrate file.
#' @param DMRfinder.output' is the output file of DMRfinder.
#' @param adjust.methods' is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'
#' @param sep seperator, default is the TAB key.
#' @param bheader a logical value indicating whether the DMRfinder.output file contains the names of the variables as its first line. By default, bheader = FALSE.
evalue.DMRfinder <- function(methyrate, DMRfinder.output, adjust.methods='BH', sep = "\t", bheader = FALSE){
  re = evalue.DMRfinder.chk (methyrate, DMRfinder.output, sep, bheader)
  return(varevalue.metilene(re$file_a, re$file_b, re$file_a_b, adjust.methods=adjust.methods));
}
