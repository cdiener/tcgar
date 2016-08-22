#  test_exon_ex.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("HuEx 1.0 ST v2")

manifest <- system.file("extdata", "manifest_legacy.tsv", package="tcgar")
d <- tempdir()

test_that("Summarized exon expression data can be loaded", {
    huex <- read_huex(manifest, d, progress=FALSE)
    expect_equal(names(huex), c("assay", "samples", "features"))
    expect_equal(class(huex$assay), "matrix")
    expect_equal(ncol(huex$assay), 50)
    expect_equal(names(huex$samples), c("name", "barcode", "label", "id", "panel", "tumor"))
    expect_equal(names(huex$features), "symbol")
})
