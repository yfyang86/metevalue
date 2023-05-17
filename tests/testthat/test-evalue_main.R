UNITEST_evalue_main_methylkit <- function() {
  data(demo_methylkit_methyrate)
  data(demo_methylkit_met_all)
  example_tempfiles <- tempfile(c("rate_combine", "methylKit_DMR_raw"))
  tempdir()
  #### write to temp file ####
  write.table(demo_methylkit_methyrate,
    file = example_tempfiles[1], row.names = FALSE,
    col.names = TRUE, quote = FALSE, sep = "\t"
  )
  write.table(demo_methylkit_met_all,
    file = example_tempfiles[2],
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  #### compute e-value and its adjustment ####
  result <- metevalue.methylKit(example_tempfiles[1],
    example_tempfiles[2],
    bheader = TRUE
  )
  return(sprintf("%0.5f", result[1, "e_value"]))
}


UNITEST_metevalue_RNA_general <- function() {
  data("demo_desq_out")
  result = metevalue.RNA_general(demo_desq_out, 'treated','untreated')[1,]
  return(sprintf("%0.5f", result[1, "evalue_all"]))
}

UNITEST_evalue_main_biseq <- function() {
  data("demo_biseq_methyrate")
  data("demo_biseq_DMR")
  example_tempfiles <- tempfile(c("demo_biseq_methyrate", "demo_biseq_DMR"))
  tempdir()
  #### write to temp file ####
  write.table(demo_biseq_methyrate,
    file = example_tempfiles[1], row.names = FALSE,
    col.names = TRUE, quote = FALSE, sep = "\t"
  )
  write.table(demo_biseq_DMR,
    file = example_tempfiles[2],
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  #### compute e-value and its adjustment ####
  result <- metevalue.biseq(example_tempfiles[1],
    example_tempfiles[2],
    bheader = TRUE
  )
  return(sprintf("%0.5f", result[1, "e_adjust"]))
}



test_that("Metevalue works", {
  expect_equal(UNITEST_evalue_main_methylkit(), "18.74770")
  expect_equal(UNITEST_evalue_main_biseq(), "28077.12788")
  expect_equal(UNITEST_metevalue_RNA_general(), "0.87551")
})
