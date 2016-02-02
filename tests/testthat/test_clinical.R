#  test_clinical.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

context("Clinical")

gbm <- system.file("extdata", "GBM", package="tcgar")
other <- system.file("extdata", "", package="tcgar")

test_that("Exon expression data can be loaded", {
    clin <- read_clinical(gbm)
    expect_error(read_clinical(other))
    expect_equal(names(clin), c("patients", "samples"))
    expect_true("barcode" %in% names(clin$patients))
    expect_true(all(clin$samples$cancer_panel == "GBM"))
    expect_true(is.logical(clin$samples$tumor))
})
