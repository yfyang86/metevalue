#' Evalue of the BiSeq data format
#'
#' Perform the Evaluation for the BiSeq data.
#' #' @param methyrate is the methyrate file.
#' @param BiSeq.output' is the output file of BiSeq
#' @param adjust.methods' is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'
#' @param sep seperator, default is the TAB key.
#' @param bheader a logical value indicating whether the BiSeq.output file contains the names of the variables as its first line. By default, bheader = FALSE.
evalue.biseq <- function(methyrate, BiSeq.output, adjust.methods='BH', sep = "\t", bheader = FALSE){
  re = evalue.biseq.chk (methyrate, BiSeq.output, sep, bheader)
  return(varevalue.metilene(re$file_a, re$file_b, re$file_a_b, adjust.methods=adjust.methods));
}
