## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

context("file info")

GDC_BASE <- "https://gdc-api.nci.nih.gov/files?"
fields <- c("file_id", "file_name", "md5sum",
    "cases.case_id", "cases.submitter_id", "cases.samples.sample_id",
    "cases.samples.submitter_id", "data_category", "data_type",
    "experimental_strategy", "updated_datetime")

query <- sprintf("%sfields=%s&pretty=true",
    GDC_BASE, paste(fields, collapse=","))

test_that("we can get full info", {
    files <- list_files(query=query, unique=FALSE, chunk=32, max_size=96)
    expect_equal(nrow(files), 96)
    expect_true("cases" %in% names(files))
    expect_equal(class(files$cases[[1]]), "data.frame")
})

test_that("we can get unique sample mappings", {
    files <- list_files(query=query, unique=TRUE, chunk=32, max_size=96)
    expect_true(nrow(files) <= 96)
    expect_false("cases" %in% names(files))
    expect_true("patient_uuid" %in% names(files))
    expect_true("patient_barcode" %in% names(files))
    expect_true("sample_uuid" %in% names(files))
    expect_true("sample_barcode" %in% names(files))
})
