#' Calculate E-value of the DMRfinder data format
#'
#' @param methyrate is the methyrate file. For example:
#' \tabular{rrrrrrrr}{
#' chr	\tab  pos	 \tab   g1	\tab ...  \tab  g1 \tab  g2 \tab ... \tab g2 \cr
#' chr1 \tab  1    \tab  0.1 \tab  ... \tab   0.1\tab  0.2\tab ... \tab 0.2\cr
#' }
#' The columns are (in order):
#' 
#'     - chr: Chromosome
#'
#'     - pos: int Position
#'
#'     - g1~g2: methylation rate data in groups
#'
#' @param DMRfinder.output is the output file of DMRfinder.
#'
#'     - chr: Chromosome
#'
#'     - start: The positions of the start sites of the corresponding region
#'
#'     - end: The positions of the end sites of the corresponding region
#'
#'     - CpG: The number of CpG sites within the corresponding region
#'
#'     - Control.mu: The average methylation rate in control group
#'
#'     - Expt1.mu: The average methylation rate in experiment group
#'
#'     - Control.Expt1.diff: The methylation difference between control and experiment groups
#'
#'     - Control.Expt1.pval: P-value based on Wald-test.
#'
#'
#' @param adjust.methods is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'
#' @param sep seperator, default is the TAB key.
#' @param bheader a logical value indicating whether the DMRfinder.output file contains the names of the variables as its first line. By default, bheader = FALSE.
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
#'
#' @examples
#' data(demo_DMRfinder_rate_combine)
#' data(demo_DMRfinder_DMRs)
#' #example_tempfiles = tempfile(c("rate_combine", "DMRfinder_out"))
#' #tempdir()
#' #write.table(demo_DMRfinder_rate_combine, file=example_tempfiles[1],
#' #      row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
#' #write.table(demo_DMRfinder_DMRs, file=example_tempfiles[2],
#' #      sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
#' #result = metevalue._DMRfinder(example_tempfiles[1], example_tempfiles[2],
#' #      bheader = TRUE)
#' #head(result)
metevalue.DMRfinder <- function(methyrate, DMRfinder.output, adjust.methods='BH', sep = "\t", bheader = FALSE){
  re = metevalue.DMRfinder.chk (methyrate, DMRfinder.output, sep, bheader)
  return(varevalue.metilene(re$file_a, re$file_b, re$file_a_b, adjust.methods=adjust.methods));
}
