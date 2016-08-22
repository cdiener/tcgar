## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

context("panel stats")

manifest <- system.file("extdata", "manifest_legacy.tsv", package="tcgar")
d <- tempdir()

test_that("panel stats can be obtained", {
    ps <- panel_samples(rnaseq=read_rnaseq_legacy(manifest, d, progress=FALSE),
        huex=read_huex(manifest, d, progress=FALSE),
        clinical=read_clinical(manifest, d, progress=FALSE))
    expect_true(all(ps[1,-1] == c(2, 50, 1)))
})
