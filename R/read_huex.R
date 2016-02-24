#  template.R
#  
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

HUEX_RE <- "lbl\\.gov_(\\w+)\\.HuEx-1_0-st-v2\\.\\d+.+\\.txt"
PANEL_RE <- "lbl\\.gov_(\\w+)\\."
HUEX_SDRF_COLS <- c(13, 14, 17, 19, 24)
HUEX_SDRF_NAMES <- c("name", "barcode", "label", "id", "panel")

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
    files <- list.files(path=folder, pattern=HUEX_RE, recursive=T, full.names=T) 
    if (length(files) == 0) stop("No HuEx exon transcription data found!")
    
    sdrf_path <- list.files(path=folder, pattern="\\.HuEx-1_0-st-v2\\.sdrf\\.txt", 
        full.names=TRUE)
    sdrf <- fread(sdrf_path, na.strings="->", select=HUEX_SDRF_COLS, 
        colClasses="character")
    names(sdrf) <- HUEX_SDRF_NAMES
    is_tumor <- as.numeric(sapply(sdrf$barcode, substr, 14, 15)) < 10
    sdrf$tumor <- is_tumor
    sdrf$panel <- toupper(str_match(sdrf$panel, PANEL_RE)[,2])
    sdrf[, "barcode" := sapply(barcode, substr, 0, 16)]
    
    if (features == "genes")
        files <- grep("HuEx-1_0-st-v2.\\d+\\.gene\\.txt", files, value=TRUE)
    else 
        files <- grep("HuEx-1_0-st-v2\\.\\d+\\.FIRMA\\.txt", files, value=TRUE)
        
    arrays <- lapply(files, function(fi) {
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
    arrays <- arrays[, sort(colnames(arrays))]

    sdrf <- sdrf[order(id)]
    
    if (any(as.numeric(sdrf$id) != as.numeric(colnames(arrays)))) 
        stop("Sample IDs do not match between assay and annotation!")
        
    colnames(arrays) <- sdrf$barcode
    
    if (any(huex_bm$symbol != rownames(arrays)))
        stop("Something is wrong with the features!")
    
    return(list(assay=arrays, samples=sdrf, features=huex_bm))
}

#' Gene annotations for the TCGA HuEx data.
#'
#' This data set contains additional annotations for the HuEx features obtained
#' from BioMart.
#'
#' @format A data frame with 18632 rows and 4 variables:
#' \describe{
#'   \item{ensgene}{Ensembl Gene ID}
#'   \item{description}{Description of the Gene}
#'   \item{symbol}{The Gene symbol}
#'   \item{entrez}{Entrez Gene ID}
#' }
"huex_bm"
