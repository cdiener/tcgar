#  read_rnaseq.R
#
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

RNA_SDRF_COLS <- c(1, 2, 15, 26)
RNA_SDRF_NAMES <- c("name", "barcode", "genome", "panel")

# To fix stupid CRAN notes
utils::globalVariables(c("gene_id", "name"))

hiseq_or_last <- function(files) {
    hiseq <- grep("HiSeq", files, value=T)
    if (length(hiseq) > 0) return(hiseq[1]) else return(files[length(files)])
}

remove_dupes <- function(files, pref=hiseq_or_last) {
    ids <- str_match(files, "unc\\.edu\\.([0-9a-z-]+)\\.")[,2]
    dupes <- unique(ids[duplicated(ids)])
    dupe_files <- files[ids %in% dupes]
    clean_files <- tapply(dupe_files, ids[ids %in% dupes], pref)
    files <- files[!(ids %in% dupes)]
    ids <- ids[!(ids %in% dupes)]
    files <- c(files, clean_files)
    ids <- c(ids, dupes)
    files <- files[order(ids)]

    return(list(files=files, ids=ids))
}

#' Reads RNA sequencing data from TCGA for several samples.
#'
#' Reads RNASeqV2 Level 3 data downloaded from TCGA. The functions does some
#' conservative checking in order to validate that all samples and features are
#' annotated correctly. Also note that duplicate data (files with the same
#' extract name) are dealt with implicitly. This will for instance occur if you
#' have data for several Illumina technologies in the same directory (for TCGA
#' mostly IlluminaGA or IlluminaHiSeq). The tie breaker used is to prefer HiSeq
#' or use the last of the sorted file names if HiSeq is not available.
#'
#' @seealso \code{\link{read_huex}} to read exon expression data.
#' @export
#' @keywords TCGA read RNASeq
#' @param folder Folder where the data files reside.
#' @param features The feature type. Must be one of "genes", "isoforms",
#'  "junctions" or "exons".
#' @param normalization The normalization method. Must be one of
#'  \describe{
#'      \item{"raw"}{The raw counts.}
#'      \item{"Q75"}{Normalization by dividing through the 75% quantile of the
#'      raw counts (TCGA default). This only available for transcripts and
#'      isoforms.}
#'      \item{"XPM"}{X per million. This will use transcripts per million (TPM)
#'      for genes and isoforms and reads per kilobase of transcript per million
#'      (RPKM) for junctions and exons.}
#'  }
#' @param progress Logical. Show progress info?
#' @return A data table containing the features as rows and the samples in the
#'  columns.
#' @examples
#' gbm <- system.file("extdata", "GBM", package = "tcgar")
#' rna <- read_rnaseq(gbm)
#'
#' @importFrom stringr str_match
#' @importFrom data.table fread set tstrsplit ':=' rbindlist
read_rnaseq <- function(folder, features="genes", normalization="raw",
    progress=FALSE) {
    files <- list.files(folder, recursive=TRUE, full.names=TRUE)

    mult <- 1
    if (features == "genes") {
        ann_idx <- c(1,4)
        if (normalization == "raw") {
            files <- grep(".genes.results", files, value=TRUE)
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- grep(".genes.results", files, value=TRUE)
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- grep(".genes.normalized_results", files, value=TRUE)
           data_idx <- 2
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "isoforms") {
        ann_idx <- 1
        if (normalization == "raw") {
            files <- grep(".isoforms.results", files, value=TRUE)
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- grep(".isoforms.results", files, value=TRUE)
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- grep(".isoforms.normalized_results", files, value=TRUE)
            data_idx <- 2
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "junctions") {
        ann_idx <- 1
        files <- grep(".junction_quantification.txt", files, value=TRUE)
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    } else if (features == "exons") {
        ann_idx <- c(1,3)
        files <- grep(".exon_quantification.txt", files, value=TRUE)
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    }
    else stop("Not a valid feature type.")

    files <- remove_dupes(files)
    ids <- files$ids
    files <- files$files

    sdrf_path <- list.files(path=folder, pattern="RNASeqV2.+\\.sdrf\\.txt",
        full.names=TRUE)

    sdrf <- lapply(sdrf_path, function(p) {
        sdrf <- fread(p, na.strings="->", select=RNA_SDRF_COLS,
            colClasses="character")
        names(sdrf) <- RNA_SDRF_NAMES
        sdrf <- unique(sdrf, by="name")
        is_tumor <- as.numeric(sapply(sdrf$barcode, substr, 14, 15)) < 10
        sdrf$tumor <- is_tumor
        sdrf$panel <- toupper(str_match(sdrf$panel, "unc\\.edu_(\\w+)\\.")[,2])
        sdrf[, "barcode" := sapply(barcode, substr, 0, 16)]
    })
    sdrf <- rbindlist(sdrf)
    sdrf <- unique(sdrf, by="name")

    # Some of the public downloads are missing files, also order the names
    sdrf <- sdrf[name %in% ids]
    sdrf <- sdrf[order(name)]

    feat <- fread(files[1], select=ann_idx)
    if (features == "genes") {
        feat[, c("symbol", "entrez") := tstrsplit(gene_id, "|", fixed=TRUE)]
        feat[, "gene_id" := NULL]
    }

    if (any(sdrf$name != sort(ids))) stop("Some samples not annotated!")

    counts <- matrix(0, nrow=nrow(feat), ncol=length(files))  # pre-allocate
    for(i in 1:length(files)) {
        if (progress) {
            cat("                                                           \r")
            cat(sprintf("Reading experiment %d/%d...", i, length(files)))
        }
        counts[, i] <- mult * fread(files[i], select=data_idx)[[1]]
    }
    if (progress) cat("\n")
    dimnames(counts) <- list(feat[[1]], sdrf$barcode)

    return(list(counts=counts, features=feat, samples=sdrf))
}

#' Gene annotations for the TCGA RNASeqV2 data.
#'
#' This data set contains additional annotations for the RNASeq genes
#' obtained from BioMart.
#'
#' @format A data frame with 20531 rows and 4 variables:
#' \describe{
#'   \item{ensgene}{Ensembl Gene ID}
#'   \item{description}{Description of the Gene}
#'   \item{symbol}{The Gene symbol}
#'   \item{entrez}{Entrez Gene ID}
#' }
"rnaseq_bm"
