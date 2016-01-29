#  read_rnaseq.R
#  
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

#UUID_RE = "\\.([0-9a-f]+-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+)\\."

#' Reads RNA sequencing data from TCGA for several samples.
#'
#' More description...
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
#' @return A data table containing the features as rows and the samples in the
#'  columns.
#' @examples
#' NULL
#'
#' @importFrom stringr str_match
#' @importFrom data.table fread set tstrsplit ':='
read_rnaseq <- function(folder, features="genes", normalization="raw") {
    if (!file.exists(file.path(folder, "file_manifest.txt")))
        stop("Not a valid TCGA download folder!")
    
    files <- fread(file.path(folder, "file_manifest.txt"))
    files <- setNames(files, gsub("\\s+","_",tolower(names(files))))
    if (length(grep("RNASeqV2", files$platform)) == 0) {
        stop("No RNASeq data found (need version 2)!")
    }
    
    mult <- 1
    if (features == "genes") {
        ann_idx <- c(1,4)
        if (normalization == "raw") {
            files <- files[grep(".genes.results", file_name)]
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- files[grep(".genes.results", file_name)]
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- files[grep(".genes.normalized_results", file_name)]
           data_idx <- 2 
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "isoforms") {
        ann_idx <- 1
        if (normalization == "raw") {
            files <- files[grep(".isoforms.results", file_name)]
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- files[grep(".isoforms.results", file_name)]
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- files[grep(".isoforms.normalized_results", file_name)]
            data_idx <- 2 
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "junctions") {
        ann_idx <- 1
        files <- files[grep(".junctions_quantification.txt", file_name)]
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    } else if (features == "exons") {
        ann_idx <- c(1,3)
        files <- files[grep(".exon_quantification.txt", file_name)]
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    }
    else stop("Not a valid feature type.")
    
    #tcga_uuid <- str_match(files$file_name, UUID_RE)[,2]
    #files[, "uuid" := tcga_uuid]
    files[, "barcode" := substr(barcode, 1, 16)]
    is_tumor <- as.numeric(tstrsplit(files$sample, "-", fixed=TRUE)[[4]]) < 10
    files[, "tumor" := is_tumor]
    files[, "sample" := NULL]
    
    feat <- fread(file.path(folder, files$file_name[1]), select=ann_idx)
    if (features == "genes") {
        feat[, c("symbol", "entrez") := tstrsplit(gene_id, "|", fixed=TRUE)]
        feat[, "gene_id" := NULL]
    }
    
    counts <- matrix(0, nrow=nrow(feat), ncol=nrow(files))  # pre-allocate
    for(i in 1:nrow(files)) {
        cat("                                                               \r")
        cat(sprintf("Reading experiment %d/%d...", i, nrow(files)))
        counts[, i] <- mult * fread(file.path(folder, files$file_name[i]), 
            select=data_idx)[[1]]
    }
    cat("\n")
    dimnames(counts) <- list(feat$entrez, files$barcode)
    
    return(list(counts=counts, features=feat, samples=files))
}
