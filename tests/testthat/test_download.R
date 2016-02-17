#  test_download.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

context("Data download")

test_that("Panel list can be obtained", {
    panels <- get_panels()
    expect_true("GBM" %in% panels)
    expect_true("LUSC" %in% panels)
    expect_true(all(grepl("[A-Z_]{2,8}", panels)))
})

test_that("Downloading data works", {
    d <- file.path(tempdir(), "GBM")
    get_data("GBM", tech="clinical", out=d, quiet=TRUE)
    expect_equal(length(list.files(d)), 2)
})
