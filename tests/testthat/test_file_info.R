## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

context("file info")

GDC_BASE <- "https://gdc-api.nci.nih.gov/files?"
fields <- c("file_id", "file_name", "md5sum",
    "cases.case_id", "cases.submitter_id", "data_category", "data_type",
    "experimental_strategy", "updated_datetime")

query <- sprintf("%sfields=%s&pretty=true&size=32",
    GDC_BASE, paste(fields, collapse=","))

test_that("we can get full info", {
    files <- list_files(query=query, unique_patient=FALSE)
    expect_equal(nrow(files), 32)
    expect_true("cases" %in% names(files))
    expect_equal(class(files$cases[[1]]), "data.frame")
})

test_that("we can get unique patient mappings", {
    files <- list_files(query=query, unique_patient=TRUE)
    expect_true(nrow(files) <= 32)
    expect_false("cases" %in% names(files))
    expect_true("patient_uuid" %in% names(files))
    expect_true("barcode" %in% names(files))
})
