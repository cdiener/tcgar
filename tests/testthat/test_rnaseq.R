#  test_rnaseq.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

context("RNASeqV2")

gbm <- system.file("extdata", "GBM", package="tcgar")
other <- system.file("extdata", "", package="tcgar")

test_that("RNASeq data can be loaded", {
    expect_output(rna <- read_rnaseq(gbm), "Reading experiment")
    expect_error(read_rnaseq(other))
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 1)
    expect_equal(names(rna$features), c("transcript_id", "symbol", "entrez"))
    expect_equal(nrow(rna$samples), 1)
})
