#  test_rnaseq.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("RNA-Seq")

manifest <- system.file("extdata", "manifest.tsv", package="tcgar")
d <- tempdir()

test_that("RNASeq per-gene data can be read", {
    rna <- read_rnaseq(manifest, d, progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(ncol(rna$counts), 2)
    expect_equal(names(rna$features), c("symbol", "entrez"))
    expect_equal(nrow(rna$samples), 2)
})

test_that("normalized RNASeq per-gene data can be read", {
    rna <- read_rnaseq(manifest, d, normalization="Q75", progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(ncol(rna$counts), 2)
    expect_equal(names(rna$features), c("symbol", "entrez"))
    expect_equal(nrow(rna$samples), 2)
})

# test_that("RNASeq junction data can be loaded", {
#     expect_output(rna <- read_rnaseq(gbm, features="junctions", progress=T),
#         "Reading experiment")
#     expect_equal(names(rna), c("counts", "features", "samples"))
#     expect_equal(class(rna$counts), "matrix")
#     expect_equal(dim(rna$counts)[2], 1)
#     expect_equal(names(rna$features), "junction")
#     expect_equal(nrow(rna$samples), 1)
# })
#
# test_that("RNASeq exon data can be loaded", {
#     expect_output(rna <- read_rnaseq(gbm, features="exons", progress=T),
#         "Reading experiment")
#     expect_equal(names(rna), c("counts", "features", "samples"))
#     expect_equal(class(rna$counts), "matrix")
#     expect_equal(dim(rna$counts)[2], 1)
#     expect_equal(names(rna$features), c("exon", "median_length_normalized"))
#     expect_equal(nrow(rna$samples), 1)
# })
#
# test_that("RNASeq isoform data can be loaded", {
#     expect_output(rna <- read_rnaseq(gbm, features="isoforms", progress=T),
#         "Reading experiment")
#     expect_equal(names(rna), c("counts", "features", "samples"))
#     expect_equal(class(rna$counts), "matrix")
#     expect_equal(dim(rna$counts)[2], 1)
#     expect_equal(names(rna$features), "isoform_id")
#     expect_equal(nrow(rna$samples), 1)
# })
