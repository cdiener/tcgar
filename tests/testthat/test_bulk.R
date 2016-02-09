#  test_clinical.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

context("Bulk reading")

gbm <- system.file("extdata", "GBM", package="tcgar")
other <- system.file("extdata", "", package="tcgar")

test_that("bulk reading works", {
    out <- capture.output(bulk <- read_bulk(gbm))
    out2 <- capture.output(bulk2 <- read_bulk(other))
    expect_true(is.null(bulk2))
    expect_equal(names(bulk), c("HuEx", "RNASeqV2", "clinical", "summary"))
    expect_true(!is.null(grep("Samples", out)))
})
