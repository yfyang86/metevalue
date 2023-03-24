#' Evalue of the Metilene data format
#' @param methyrate metilene input file path. This file is a sep (e.g. TAB) separated file with two key columns and several value columns in pairs:
#' For exampe:
#'
#' chrom	pos	g1	g1	g1	g1	g1	g1	g1	g1	g2	g2	g2	g2	g2	g2	g2	g2
#'
#' chrom and pos are keys;
#' g1 g1 g2 g2 must be stored in groups.
#'
#' @param metilene.output  metilene input file path. This file should stored as a sep(e.g. TAB) separated file with two key columns and several value columns:
#' The columns are (in order):
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
#' @param adjust.methods is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'
#' @param sep seperator, default is the TAB key.
#' @param bheader a logical value indicating whether the metilene.output file contains the names of the variables as its first line. By default, bheader = FALSE.
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
#' #### metilene example ####'
#' data(demo_metilene_input)
#' data(demo_metilene_out)
metevalue.metilene <- function(methyrate, metilene.output, adjust.methods='BH', sep = "\t", bheader = FALSE){
    re = metevalue.metilene.chk (methyrate, metilene.output, sep, bheader)
return(varevalue.metilene(re$file_a, re$file_b, re$file_a_b, adjust.methods=adjust.methods));
}


