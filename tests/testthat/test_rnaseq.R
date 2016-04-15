#  test_rnaseq.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("RNASeqV2")

gbm <- system.file("extdata", "GBM", package="tcgar")
other <- system.file("extdata", "", package="tcgar")

test_that("RNASeq per-gene data can be loaded", {
    expect_output(rna <- read_rnaseq(gbm, progress=T), "Reading experiment")
    expect_error(read_rnaseq(other))
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 1)
    expect_equal(names(rna$features), c("transcript_id", "symbol", "entrez"))
    expect_equal(nrow(rna$samples), 1)
})

test_that("RNASeq junction data can be loaded", {
    expect_output(rna <- read_rnaseq(gbm, features="junctions", progress=T),
        "Reading experiment")
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 1)
    expect_equal(names(rna$features), "junction")
    expect_equal(nrow(rna$samples), 1)
})

test_that("RNASeq exon data can be loaded", {
    expect_output(rna <- read_rnaseq(gbm, features="exons", progress=T),
        "Reading experiment")
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 1)
    expect_equal(names(rna$features), c("exon", "median_length_normalized"))
    expect_equal(nrow(rna$samples), 1)
})

test_that("RNASeq isoform data can be loaded", {
    expect_output(rna <- read_rnaseq(gbm, features="isoforms", progress=T),
        "Reading experiment")
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 1)
    expect_equal(names(rna$features), "isoform_id")
    expect_equal(nrow(rna$samples), 1)
})
