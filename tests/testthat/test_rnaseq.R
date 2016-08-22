#  test_rnaseq.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("RNA-Seq")

manifest <- system.file("extdata", "manifest.tsv", package="tcgar")
d <- tempdir()

test_that("RNA-Seq counts can be read", {
    rna <- read_rnaseq(manifest, d, progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(ncol(rna$counts), 1)
    expect_equal(names(rna$features),
        c("ensgene", "description", "symbol", "entrez"))
    expect_equal(nrow(rna$samples), 1)
})

test_that("RNA-Seq FPKM can be read", {
    rna <- read_rnaseq(manifest, d, normalization="FPKM", progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(ncol(rna$counts), 1)
    expect_equal(names(rna$features),
        c("ensgene", "description", "symbol", "entrez"))
    expect_equal(nrow(rna$samples), 1)
})

test_that("RNA-Seq FPKM-UQ can be read", {
    rna <- read_rnaseq(manifest, d, normalization="FPKM-UQ", progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(ncol(rna$counts), 1)
    expect_equal(names(rna$features),
        c("ensgene", "description", "symbol", "entrez"))
    expect_equal(nrow(rna$samples), 1)
})
