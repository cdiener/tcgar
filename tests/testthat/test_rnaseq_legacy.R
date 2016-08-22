#  test_rnaseq.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("Legacy RNA-Seq")

manifest <- system.file("extdata", "manifest_legacy.tsv", package="tcgar")
d <- tempdir()

test_that("RNASeq per-gene data can be read", {
    rna <- read_rnaseq_legacy(manifest, d, progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(ncol(rna$counts), 2)
    expect_equal(names(rna$features),
        c("entrez", "ensgene", "description", "symbol"))
    expect_equal(nrow(rna$samples), 2)
})

test_that("normalized RNASeq per-gene data can be read", {
    rna <- read_rnaseq_legacy(manifest, d, normalization="Q75", progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(ncol(rna$counts), 2)
    expect_equal(names(rna$features),
        c("entrez", "ensgene", "description", "symbol"))
    expect_equal(nrow(rna$samples), 2)
})

test_that("RNASeq junction data can be loaded", {
    rna <- read_rnaseq_legacy(manifest, d, features="junctions", progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 2)
    expect_equal(names(rna$features), "junction")
    expect_equal(nrow(rna$samples), 2)
})

test_that("RNASeq exon data can be loaded", {
    rna <- read_rnaseq_legacy(manifest, d, features="exons", progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 2)
    expect_equal(names(rna$features), c("exon", "median_length_normalized"))
    expect_equal(nrow(rna$samples), 2)
})

test_that("RNASeq isoform data can be loaded", {
    rna <- read_rnaseq_legacy(manifest, d, features="isoforms", progress=FALSE)
    expect_equal(names(rna), c("counts", "features", "samples"))
    expect_equal(class(rna$counts), "matrix")
    expect_equal(dim(rna$counts)[2], 2)
    expect_equal(names(rna$features), "isoform_id")
    expect_equal(nrow(rna$samples), 2)
})
