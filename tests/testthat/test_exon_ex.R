#  test_exon_ex.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

context("HuEx 1.0 ST v2")

gbm <- system.file("extdata", "GBM", package="tcgar")
other <- system.file("extdata", "", package="tcgar")

test_that("Exon expression data can be loaded", {
    huex <- read_huex(gbm)
    expect_error(read_huex(other))
    expect_equal(names(huex), c("assay", "samples", "features"))
    expect_equal(class(huex$assay), "matrix")
    expect_equal(dim(huex$assay)[2], 29)
    expect_equal(names(huex$samples), c("name", "barcode", "label", "id", "panel", "tumor"))
})
