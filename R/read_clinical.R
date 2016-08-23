#  template.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.


COLS <- c("patient_uuid" = "//shared:bcr_patient_uuid",
    "patient_barcode" = "//shared:bcr_patient_barcode",
    "gender" = "//shared:gender", "vital" = "//clin_shared:vital_status",
    "days_to_contact" = "//clin_shared:days_to_last_followup",
    "days_to_death" = "//clin_shared:days_to_death",
    "days_to_birth" = "//clin_shared:days_to_birth",
    "panel" = "//admin:disease_code",
    "histology" = "//shared:histological_type",
    "tissue_site" = "//clin_shared:tumor_tissue_site",
    "stage" = "//shared_stage:pathologic_stage",
    "T" = "//shared_stage:pathologic_T", "N" = "//shared_stage:pathologic_N",
    "M" = "//shared_stage:pathologic_M",
    "residual_tumor" = "//clin_shared:residual_tumor")
COUNTS <- c("new_tumor_events" = "//nte:new_tumor_event_after_initial_treatment",
    "follow_ups" = "//clin_shared:bcr_followup_uuid")
TO_NUM <- c("days_to_contact", "days_to_birth", "days_to_birth")
TO_LOWER <- c("patient_uuid", "gender")

XML_RE <- "\\.xml"

#' Reads clinical data from TCGA for multiple samples.
#'
#' This function assembles a data set with the major clinical indicators for each
#' patient as well as sample procurement data for all the samples available for
#' each patient from TCGA Biotab data.
#'
#' @export
#' @keywords clinical TCGA read GDC
#' @param manifest Path to the GDC file manifest.
#' @param folder The folder that contains the data.
#' @param progress Whether to show progress information.
#' @return A data table containing the information for patients on its rows.
#' @examples
#' # Not run due to large download...
#' # gbm <- system.file("extdata", "manifest.tsv", package = "tcgar")
#' # d <- tempdir()
#' # clin <- read_clinical(gbm, d)
#'
#' @importFrom magrittr '%>%'
#' @importFrom xml2 read_xml xml_find_all xml_text xml_ns
#' @importFrom data.table rbindlist
#' @importFrom pbapply pbapply pblapply
read_clinical <- function(manifest, folder, progress=TRUE) {
    man <- fread(manifest)
    files <- man[grep(XML_RE, filename)]
    afun <- ifelse(progress, pbapply, apply)

    if (progress) cat("Reading clinical data:\n")
    patients <- afun(files, 1, function(fi) {
        xml_doc <- read_xml(file.path(folder, fi["id"], fi["filename"]))
        ns <- fix_ns(xml_ns(xml_doc))
        vals <- sapply(COLS, function(co)
            xml_doc %>% xml_find_all(co, ns=ns) %>% xml_text() %>% last_or_na()
        )
        names(vals) <- names(COLS)
        counts <- sapply(COUNTS, function(co)
            xml_doc %>% xml_find_all(co, ns=ns) %>% xml_text() %>% length()
        )
        names(counts) <- names(COUNTS)
        data.table(t(vals), t(counts))
    })

    patients <- rbindlist(patients)

    for (co in TO_NUM) set(patients, j=co, value=as.numeric(patients[[co]]))
    for (co in TO_LOWER) set(patients, j=co, value=tolower(patients[[co]]))

    return(patients)
}
