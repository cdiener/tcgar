#  read_rnaseq.R
#  
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

#' Reads RNA sequencing data from TCGA for several samples.
#'
#' More description...
#' 
#' @seealso \code{\link{bla}} for doing it slower.
#' @export
#' @keywords TCGA read
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
#' @importFrom string str_match
#' @importFrom data.table fread set tstrsplit ':='
read_rnaseq <- function(folder, features="genes", normalization="XPM") {
    if (features == "genes") {
        ann_idx <- c(1,4)
        if (normalization == "raw") {
            files <- list.files(folder, pattern=".genes.results", full.names=TRUE)
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- list.files(folder, pattern=".genes.results", full.names=TRUE)
            data_idx <- 3
        }
        else if (normalization == "Q75") {
            files <- list.files(folder, pattern=".genes.normalized_results", full.names=TRUE)
           data_idx <- 2 
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "isoforms") {
        ann_idx <- c(1,4)
        if (normalization == "raw") {
            files <- list.files(folder, pattern=".isoforms.results", full.names=TRUE)
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- list.files(folder, pattern=".isoforms.results", full.names=TRUE)
            data_idx <- 3
        }
        else if (normalization == "Q75") {
            files <- list.files(folder, pattern=".isoforms.normalized_results", full.names=TRUE)
           data_idx <- 2 
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "junctions") {
        ann_idx <- 1
        if (normalization == "raw") {
            files <- list.files(folder, pattern=".junction_quantification.txt", full.names=TRUE)
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- list.files(folder, pattern=".junction_quantification.txt", full.names=TRUE)
            data_idx <- 4
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "exons") {
        ann_idx <- c(1,3)
        if (normalization == "raw") {
            files <- list.files(folder, pattern=".genes.results", full.names=TRUE)
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- list.files(folder, pattern=".genes.results", full.names=TRUE)
            data_idx <- 4
        }
        else stop("Not a valid features/normalization combination.")
    }
    else stop("Not a valid feature type.")
    
    tcga_uuid <- str_match(files, "\\.([0-9a-f]+-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+)\\.")[,2]
    
    res <- fread(files[1], select=ann_idx)
    
    cols <- lapply(1:length(tcga_uuid), function(i) {
        co <- fread(files[i], select=data_idx)[[1]]
        res[, (tcga_uuid[i]) := co]
        })
    
    return(res)
}
