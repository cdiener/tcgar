## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

GDC_BASE <- "https://gdc-api.nci.nih.gov/files?"
fields <- c("file_id", "file_name", "md5sum",
    "cases.case_id", "cases.submitter_id", "data_category", "data_type",
    "experimental_strategy", "updated_datetime")

# To fix stupid CRAN notes
utils::globalVariables("cases")

#' Obtains infromation about the files in GDC.
#'
#' This function downloads a list of all files in GDC together with additional
#' metadata. Here the most important metadata is the mapping of individual
#' experimental data to the corresponding patient.
#'
#' @param query Optional. The GDC query to be used. By default downloads
#' detailed information for all files in the GDC.
#' @param unique_patient optional boolean. Whether only to return file info
#' for file that can be mapped to a single patient.
#' @return A data table mapping files to various ids and the the patient
#' @examples
#'  NULL
#'
#' @export
#' @importFrom jsonlite fromJSON
#' @importFrom data.table as.data.table
#' @importFrom curl curl
list_files <- function(query="default", unique_patient=TRUE) {
    if (query == "default")
        query <- sprintf("%sfields=%s&pretty=true&size=1000000",
            GDC_BASE, paste(fields, collapse=","))

    json <- fromJSON(query)
    fi <- as.data.table(json$data$hits)
    if (unique_patient) {
        fi <- fi[sapply(cases, nrow) == 1]
        fi[, c("patient_uuid", "patient_barcode") := rbindlist(cases)]
        fi[, cases := NULL]
    }

    return(fi)
}

#' File infos for the GDC release v2 from 08-2016
#'
#' This data set contains detailed information about all files in the GDC
#' release that have a corresponding read method in `tcgar`
#'
#' @format A data frame with the following columns.
#' \describe{
#'   \item{data_type}{The type of the file.}
#'   \item{updated_datetime}{Last update of the file.}
#'   \item{file_name}{The filename of the file.}
#'   \item{md5sum}{md5 hash for the file.}
#'   \item{file_id}{The GDC UUID for the file.}
#'   \item{data_category}{Category for the file.}
#'   \item{experimental_strategy}{Type of experiment used to obtain the data.}
#'   \item{patient_uuid}{Unique ID for the patient the data came from.}
#'   \item{barcode}{Barcode for the patient the data came from.}
#' }
"gdc_files"
