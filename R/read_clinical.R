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
    "histology" = "//shared:histological_type",
    "tissue_site" = "//clin_shared:tumor_tissue_site",
    "stage" = "//shared_stage:pathologic_stage",
    "T" = "//shared_stage:pathologic_T", "N" = "//shared_stage:pathologic_N",
    "M" = "//shared_stage:pathologic_M",
    "residual_tumor" = "//clin_shared:residual_tumor")
COUNTS <- c("new_tumor_events" = "//nte:new_tumor_event_after_initial_treatment",
    "follow_ups" = "//clin_shared:bcr_followup_uuid")
TO_NUM <- c("days_to_contact", "days_to_birth", "days_to_birth")

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
#' @return A data table containing the information for patients on its rows.
#' @examples
#' gbm <- system.file("extdata", "GBM", package = "tcgar")
#' clin <- read_clinical(gbm)
#'
#' @importFrom magrittr '%>%'
#' @importFrom xml2 read_xml xml_find_all xml_text
#' @importFrom data.table rbindlist
#' @importFrom utils txtProgressBar setTxtProgressBar
read_clinical <- function(manifest, folder, progress=FALSE) {
    man <- fread(manifest)
    files <- man[grep(XML_RE, filename)]

    if (progress) {
        pb <- txtProgressBar(min=0, max=nrow(files), style=3)
        i <- 0
        env <- environment()
    }

    patients <- apply(files, 1, function(fi) {
        xml_doc <- read_xml(file.path(folder, fi["id"], fi["filename"]))
        ns <- fix_ns(xml_ns(xml_doc))
        vals <- lapply(COLS, function(co)
            xml_doc %>% xml_find_all(co, ns=ns) %>% xml_text() %>% last_or_na()
        )
        names(vals) <- names(COLS)
        counts <- lapply(COUNTS, function(co)
            xml_doc %>% xml_find_all(co, ns=ns) %>% xml_text() %>% length()
        )
        names(counts) <- names(COUNTS)
        if(progress) {
            env$i <- env$i + 1
            setTxtProgressBar(env$pb, env$i)
        }

        data.table(t(vals), t(counts))
    })
    if (progress) close(pb)

    patients <- rbindlist(patients)
    for (co in TO_NUM) set(patients, j=co, as.numeric(patients[[co]]))

    return(patients)
}
