#  template.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

SAMPLE_COLS <- c(bcr_patient_uuid="patient_uuid", bcr_sample_barcode="barcode", 
    bcr_sample_uuid="sample_uuid", days_to_collection="days_to_collection", 
    sample_type="sample_type", sample_type_id="tumor")  

PATIENT_COLS <- c(bcr_patient_uuid="patient_uuid", bcr_patient_barcode="barcode", 
    gender="gender", vital_status="vital", last_contact_days_to="days_to_contact",
    death_days_to="days_to_death", birth_days_to="days_to_birth", 
    histological_type="histology", tumor_tissue_site="tissue_site", 
    clinical_M="M", clinical_N="N", clinical_T="T", clinical_stage="stage")
    
SAMPLE_RE <- "biospecimen_sample_([a-z]{4})\\.txt"

#' Reads clinical data from TCGA for multiple samples.
#'
#' This function assembles a data set with the major clinical indicators for each 
#' patient as well as sample procurement data for all the samples available for 
#' each patient from TCGA Biotab data.
#' 
#' @export
#' @keywords clinical TCGA read
#' @param folder The folder that contains the data.
#' @param add_sample_cols Names of additional columns read from normal and tumor
#'  sample data or NULL for no additional columns. If the vector is named the 
#'  names are expected to be the column names in the TCGA Biotab files and the
#'  vector elements will be the column names in the returned data table. 
#' @param add_patient_cols Names of additional columns read from patient data 
#'  or NULL for no additional columns. If the vector is named the 
#'  names are expected to be the column names in the TCGA Biotab files and the
#'  vector elements will be the column names in the returned data table.
#' @return A list with two data tables, samples and patient, containing the 
#'  information for patients and samples.
#' @importFrom stringr str_match
#' @importFrom data.table fread tstrsplit ':=' rbindlist
read_clinical <- function(folder, add_sample_cols=NULL, add_patient_cols=NULL) {
    files <- fread(file.path(folder, "file_manifest.txt"))
    files <- setNames(files, gsub("\\s+","_",names(files)))
    colnames(files) <- tolower(colnames(files))
    
    sample_files = files[grep("nationwidechildrens.org_biospecimen_sample", file_name)]
    patient_files = files[grep("nationwidechildrens.org_clinical_patient", file_name)]
    
    # Add additional columns, if no names are given orginal column names are
    # conserved
    sc <- SAMPLE_COLS
    pc <- PATIENT_COLS
    
    if (!is.null(add_sample_cols)) {
        if (is.null(names(add_sample_cols))) 
            names(add_sample_cols) <- add_sample_cols
        sc <- c(sc, add_sample_cols)
    }
    
    if (!is.null(add_patient_cols)) {
        if (is.null(names(add_patient_cols))) 
            names(add_patient_cols) <- add_patient_cols
        pc <- c(pc, add_patient_cols)
    }
    
    samples <- lapply(sample_files$file_name, function(fi) {
        fi <- file.path(folder, fi)
        cols <- colnames(fread(fi, nrows=0))
        s <- fread(fi, skip=2, select=which(cols %in% names(SAMPLE_COLS)), 
            header=F, na.strings=c("[Not Available]", "[Unknown]"))
        colnames(s) <- cols[cols %in% names(SAMPLE_COLS)]
        colnames(s) <- SAMPLE_COLS[colnames(s)]
        s$tumor <- s$tumor < 10
        s$cancer_panel <- str_match(fi, SAMPLE_RE)[1,2]
        s$sample_uuid <- tolower(s$sample_uuid)
        s$patient_uuid <- tolower(s$patient_uuid)
        s
        })
    
    patients <- lapply(patient_files$file_name, function(fi) {
        fi <- file.path(folder, fi)
        cols <- colnames(fread(fi, nrows=0))
        p <- fread(fi, skip=3, select=which(cols %in% names(PATIENT_COLS)), 
            header=F, na.strings=c("[Not Available]", "[Unknown]"))
        colnames(p) <- cols[cols %in% names(PATIENT_COLS)]
        colnames(p) <- PATIENT_COLS[colnames(p)]
        p$patient_uuid <- tolower(p$patient_uuid)
        p
        })
    
    return(list(patients=rbindlist(patients), samples=rbindlist(samples)))
}
