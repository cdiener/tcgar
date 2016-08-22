#  test_download.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("Data download")

manifest <- system.file("extdata", "manifest.tsv", package="tcgar")
legacy <- system.file("extdata", "manifest_legacy.tsv", package="tcgar")
d <- tempdir()

test_that("Downloading data works", {
    man <- fread(manifest)
    get_data(manifest, d, quiet=TRUE)
    folders <- list.dirs(path=d, full.names=FALSE, recursive=FALSE)
    expect_true(all(man$id %in% folders))
})

test_that("Downloading legacy data works", {
    man <- fread(legacy)
    get_data(legacy, d, quiet=TRUE)
    folders <- list.dirs(path=d, full.names=FALSE, recursive=FALSE)
    expect_true(all(man$id %in% folders))
})
