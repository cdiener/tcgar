#  test_download.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

context("Data download")

manifest <- system.file("extdata", "manifest.tsv", package="tcgar")
d <- tempdir()

test_that("Downloading data works", {
    man <- data.table::fread(manifest)
    get_data(manifest, d, quiet=TRUE)
    folders <- list.dirs(path=d, full.names=FALSE, recursive=FALSE)
    expect_true(all(man$id %in% folders))
})
