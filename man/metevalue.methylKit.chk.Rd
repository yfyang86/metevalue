% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalue_data_chk.R
\name{metevalue.methylKit.chk}
\alias{metevalue.methylKit.chk}
\title{Check the methylKit data format}
\usage{
metevalue.methylKit.chk(
  input_filename_a,
  input_filename_b,
  sep = "\\t",
  bheader = FALSE
)
}
\arguments{
\item{input_filename_a}{the combined data of methylation rate file. This file is a sep (e.g. TAB) separated file with two key columns and several value columns:
For exampe:

\tabular{rrrrrrrr}{
chr    \tab  pos     \tab   g1    \tab ...  \tab  g1 \tab  g2 \tab ... \tab g2 \cr
chr1 \tab  1    \tab  0.1 \tab  ... \tab   0.1\tab  0.2\tab ... \tab 0.2\cr
}

     - chr and pos are keys;

     - g1~g2: methylation rate data in groups.}

\item{input_filename_b}{the output file of methylKit. a methylDiff or methylDiffDB object containing the differential methylated locations satisfying the criteria.
The columns are (in order):

    - chr:   Chromosome

    - start: The position of the start sites of the corresponding region

    - end: The position of the end sites of the corresponding region

    - strand: Strand

    - p: p-value

    - qvalue:  The adjusted p-value based on BH method

    - meth.diff : The difference between the group means of methylation level}

\item{sep}{separator, default is the TAB key.}

\item{bheader}{a logical value indicating whether the input_filename_b file contains the names of the variables as its first line. By default, bheader = FALSE.}
}
\value{
list(file_a, file_b, file_a_b) returns a list with three pr-handled data.frames corresponding to the input_filename_a, input_filename_b file and a A JOIN B file.
}
\description{
Check the methylKit data format
}
\examples{
data(demo_methylkit_methyrate)
data(demo_methylkit_met_all)
example_tempfiles = tempfile(c("rate_combine", "methylKit_DMR_raw"))
tempdir()
write.table(demo_methylkit_methyrate, file=example_tempfiles[1],
      row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.table(demo_methylkit_met_all, file=example_tempfiles[2],
      sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
result = metevalue.methylKit.chk(example_tempfiles[1], example_tempfiles[2],
      bheader = TRUE)
}
