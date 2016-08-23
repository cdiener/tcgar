#  read_rnaseq.R
#
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

RNA_SDRF_COLS <- c(1, 2, 15, 26)
RNA_SDRF_NAMES <- c("name", "barcode", "genome", "panel")
RNA_PANEL_RE <- "unc\\.edu_(\\w+)\\."
RNA_ID_RE <- "unc\\.edu\\.([a-z0-9-]+)\\."

# To fix stupid CRAN notes
utils::globalVariables(c("gene_id", "name", "sample_uuid"))

#' Reads RNA sequencing data from the GDC legacy archive.
#'
#' Reads RNASeq Level 3 data downloaded from the legacy GDC. The functions does some
#' conservative checking in order to validate that all samples and features are
#' annotated correctly. Also note that duplicate data (files with the same
#' extract name) will be treated as independent duplicates for the same sample.
#'
#' @seealso \code{\link{read_huex}} to read exon expression data.
#' @export
#' @keywords TCGA read RNASeq GDC
#' @param manifest Path to the GDC file manifest.
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
#' # Not run due to large download...
#' # gbm <- system.file("extdata", "manifest.tsv", package = "tcgar")
#' # d <- tempdir()
#' # rna <- read_rnaseq(gbm, d)
#'
#' @importFrom stringr str_match
#' @importFrom data.table fread set tstrsplit ':=' rbindlist setkey
#' @importFrom pbapply startpb setpb closepb
read_rnaseq_legacy <- function(manifest, folder, features="genes", normalization="raw",
    progress=TRUE) {
    man <- fread(manifest)
    afun <- ifelse(progress, pbapply, apply)

    mult <- 1
    if (features == "genes") {
        ann_idx <- 1
        nam <- "entrez"
        if (normalization == "raw") {
            files <- man[grep("\\.genes\\.results", filename)]
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- man[grep("\\.genes\\.results", filename)]
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- man[grep("\\.genes\\.normalized_results", filename)]
            data_idx <- 2
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "isoforms") {
        ann_idx <- 1
        nam <- "isoform_id"
        if (normalization == "raw") {
            files <- man[grep("\\.isoforms\\.results", filename)]
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- man[grep("\\.isoforms\\.results", filename)]
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- man[grep("\\.isoforms\\.normalized_results", filename)]
            data_idx <- 2
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "junctions") {
        ann_idx <- 1
        nam <- "junction"
        files <- man[grep("\\.junction_quantification\\.txt", filename)]
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    } else if (features == "exons") {
        ann_idx <- c(1,3)
        nam <- "exon"
        files <- man[grep("\\.exon_quantification\\.txt", filename)]
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    }
    else stop("Not a valid feature type.")

    files[, sample_uuid := str_match(filename, RNA_ID_RE)[,2]]

    if(progress) cat("Reading annotations:\n")
    sdrf <- afun(files, 1, function(fi) {
        path <- file.path(folder, fi["id"])
        mage <- find_dir(path, "mage-tab")
        if (length(mage) == 0) {
            mage <- list.files(path=path, pattern="mage-tab.+\\.tar\\.gz",
                full.names=TRUE)
            if (length(mage) != 1) {
                warning(sprintf("ID %s has no or several mage tab annotations!",
                    fi["id"]))
                return(NULL)
            }
            err <- untar(mage[1], exdir=path.expand(path))
            if (err != 0) stop(sprintf("Error unpacking file %s!", mage))
            mage <- find_dir(path, "mage-tab")
        }
        path <- mage
        sdrf_path <- list.files(path=path, pattern="\\.sdrf\\.txt", full.names=TRUE)
        fread(sdrf_path, na.strings="->", select=RNA_SDRF_COLS,
            colClasses="character")
    })
    sdrf <- rbindlist(sdrf)
    names(sdrf) <- RNA_SDRF_NAMES
    sdrf <- unique(sdrf, by="name")
    is_tumor <- as.numeric(sapply(sdrf$barcode, substr, 14, 15)) < 10
    sdrf$tumor <- is_tumor
    sdrf$panel <- toupper(str_match(sdrf$panel, RNA_PANEL_RE)[,2])
    sdrf[, "sample_barcode" := sapply(barcode, substr, 0, 16)]
    sdrf[, "patient_barcode" := sapply(barcode, substr, 0, 12)]
    sdrf[, barcode := NULL]

    # Some of the public downloads are missing files, also order the names
    setkey(sdrf, "name")
    files <- files[sample_uuid %in% sdrf$name]
    sdrf <- sdrf[files$sample_uuid]
    setkey(sdrf, NULL)

    if (any(sdrf$name != files$sample_uuid))
        stop("Missing sample annotations for some IDs.")

    feat <- fread(files[1, file.path(folder, id, filename)], select=ann_idx)
    if (features == "genes") {
        feat[, "entrez" := sub(".+\\|", "", gene_id)]
        feat[, "gene_id" := NULL]
        feat <- merge(feat, unique(genemap, by="entrez"), by="entrez",
            all.x=TRUE)
    }

    if(progress) {
        cat("Reading assays:\n")
        pb <- startpb(min=0, max=nrow(files))
        on.exit(closepb(pb))
    }
    counts <- matrix(0, nrow=nrow(feat), ncol=nrow(files))  # pre-allocate
    for(i in 1:nrow(files)) {
        path <- files[i, file.path(folder, id, filename)]
        if (progress) setpb(pb, i)
        counts[, i] <- mult * fread(path, select=data_idx)[[1]]
    }
    if (normalization == "raw") counts <- round(counts)
    dimnames(counts) <- list(feat[[nam]], files$sample_uuid)

    return(list(counts=counts, features=feat, samples=sdrf))
}
