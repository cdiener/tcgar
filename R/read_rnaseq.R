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
#' @importFrom data.table fread set tstrsplit ':=' setNames
read_rnaseq <- function(folder, features="genes", normalization="XPM") {
    files <- fread(file.path(folder, "file_manifest.txt"))
    files <- setNames(files, gsub("\\s+","_",names(files)))
    
    mult <- 1
    if (features == "genes") {
        ann_idx <- c(1,4)
        if (normalization == "raw") {
            files <- files[grep(".genes.results", File_Name)]
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- files[grep(".genes.results", File_Name)]
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- files[grep(".genes.normalized_results", File_Name)]
           data_idx <- 2 
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "isoforms") {
        ann_idx <- c(1,4)
        if (normalization == "raw") {
            files <- files[grep(".isoforms.results", File_Name)]
            data_idx <- 2
        } else  if (normalization == "XPM") {
            files <- files[grep(".isoforms.results", File_Name)]
            data_idx <- 3
            mult <- 1e6
        }
        else if (normalization == "Q75") {
            files <- files[grep(".isoforms.normalized_results", File_Name)]
           data_idx <- 2 
        }
        else stop("Not a valid features/normalization combination.")
    } else if (features == "junctions") {
        ann_idx <- 1
        files <- files[grep(".junctions_quantification.txt", File_Name)]
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    } else if (features == "exons") {
        ann_idx <- c(1,3)
        files <- files[grep(".exon_quantification.txt", File_Name)]
        if (normalization == "raw") data_idx <- 2
        else if (normalization == "XPM") data_idx <- 4
        else stop("Not a valid features/normalization combination.")
    }
    else stop("Not a valid feature type.")
    
    tcga_uuid <- str_match(files$File_Name, "\\.([0-9a-f]+-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+)\\.")[,2]
    files[, "uuid" := tcga_uuid]
    is_tumor <- as.numeric(tstrsplit(files$Sample, "-")[[4]]) < 10
    files[, "Tissue_Type" := factor(c("normal", "tumor")[is_tumor + 1])]
    
    features <- fread(file.path(folder, files$File_Name[1]), select=ann_idx)
    
    counts <- matrix(0, nrow=nrow(features), ncol=nrow(files))
    for(i in 1:nrow(files)) {
        counts[, i] <- mult*fread(file.path(folder, files$File_Name[i]), 
        select=data_idx)[[1]]
    }
    
    return(list(counts=round(counts), features=features, samples=files))
}
