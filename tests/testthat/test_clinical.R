#  test_clinical.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("Clinical")

manifest <- system.file("extdata", "manifest.tsv", package="tcgar")
d <- tempdir()

test_that("Clinical data can be loaded", {
    clin <- read_clinical(manifest, d, progress=FALSE)
    expect_equal(ncol(clin), 17)
    expect_equal(nrow(clin), 1)
    expect_type(clin$follow_ups, "integer")
})
