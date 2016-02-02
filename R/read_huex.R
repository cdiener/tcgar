#  template.R
#  
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

SDRF_COLS <- c(13, 14, 17, 19)
SDRF_NAMES <- c("name", "barcode", "label", "id")

#' Reads HuEx 1.0 st v2 exon expression data from the TCGA project.
#'
#' This reads level 3 data for exon expression data for the TCGA project.
#' The function is currently limited to level 3 data since lower level data
#' are not publically available.
#' 
#' @seealso \code{\link{read_rnaseq}} to read RNASeq data.
#' @export
#' @keywords HuEx microarray TCGA
#' @param folder Location of the directory containing downloaded TCGA data.
#' @param features Which features to use. Can be either "genes" or "exons".
#' @return A list with two elements, a matrix of log-expression values and a
#'  data table giving additional information for samples.
#' @examples
#' gbm <- system.file("extdata", "GBM", package = "tcgar")
#' huex <- read_huex(gbm)
#'
#' @importFrom data.table fread set ':='
read_huex <- function(folder, features="genes") {
    if (!file.exists(file.path(folder, "file_manifest.txt")))
        stop("Not a valid TCGA download folder (missing file manifest)!")
    
    files <- fread(file.path(folder, "file_manifest.txt"))
    files <- setNames(files, gsub("\\s+","_",tolower(names(files))))
    if (length(grep("HuEx-1_0-st-v2", files$platform)) == 0)
        stop("No HuEx exon transcription data found!")
    
    sdrf_path <- files[grep(".HuEx-1_0-st-v2.sdrf.txt", file_name), file_name]
    sdrf_path <- file.path(folder, sdrf_path)
    sdrf <- fread(sdrf_path, na.strings="->", select=SDRF_COLS, 
        colClasses="character")
    names(sdrf) <- SDRF_NAMES
    is_tumor <- as.numeric(sapply(sdrf$barcode, substr, 14, 15)) < 10
    sdrf$tumor <- is_tumor
    sdrf[, "barcode" := sapply(barcode, substr, 0, 15)]
    
    if (features == "genes")
        files <- files[grep(".HuEx-1_0-st-v2.\\d+.gene.txt", file_name)]
    else 
        files <- files[grep(".HuEx-1_0-st-v2.\\d+.FIRMA.txt", file_name)]
        
    arrays <- lapply(files$file_name, function(fi) {
        fi <- file.path(folder, fi)
        cols <- colnames(fread(fi, nrows=0, header=T))
        arr <- fread(fi, skip=2)
        rows <- arr[[1]]
        set(arr, j=1L, value=NULL)
        arr <- as.matrix(arr)
        rownames(arr) <- rows
        colnames(arr) <- cols[-1]
        arr
    })
    arrays <- do.call("cbind", arrays)
    
    return(list(assay=arrays, samples=sdrf))
}
