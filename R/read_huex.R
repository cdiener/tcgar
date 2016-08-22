#  template.R
#
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

HUEX_GENE_RE <- "HuEx-1_0-st-v2.\\d+\\.gene\\.txt"
HUEX_EXON_RE <- "HuEx-1_0-st-v2\\.\\d+\\.FIRMA\\.txt"
HUEX_SDRF_RE <- "\\.HuEx-1_0-st-v2\\.sdrf\\.txt"
PANEL_RE <- "lbl\\.gov_(\\w+)\\."
HUEX_SDRF_COLS <- c(13, 14, 17, 19, 24)
HUEX_SDRF_NAMES <- c("name", "barcode", "label", "id", "panel")

# To fix stupid CRAN notes
utils::globalVariables(c("barcode", "id", "filename", "genemap", "gdc_files"))

#' Reads HuEx 1.0 st v2 exon expression data.
#'
#' This reads level 3 data for exon expression data.
#'
#' @seealso \code{\link{read_rnaseq}} to read RNASeq data.
#' @export
#' @keywords HuEx microarray TCGA GDC
#' @param manifest Path to the GDC file manifest.
#' @param folder Location of the directory containing downloaded TCGA data.
#' @param features Which features to use. Can be either "genes" or "exons".
#' @param progress Whether to show progress information.
#' @return A list with three elements, a matrix of log-expression values, a
#'  data table giving additional information for samples and a data table
#'  describing the features (genes or exon IDs).
#' @examples
#' # Not run due to large download...
#' # gbm <- system.file("extdata", "manifest.tsv", package = "tcgar")
#' # d <- tempdir()
#' # huex <- read_huex(gbm, d)
#'
#' @importFrom data.table data.table fread set ':='
read_huex <- function(manifest, folder, features="genes", progress=TRUE) {
    man <- fread(manifest)
    if (features == "genes") {
        files <- man[grep(HUEX_GENE_RE, filename)]
    }
    else
        files <- man[grep(HUEX_EXON_RE, filename)]
    if (nrow(files) == 0) stop("No HuEx exon transcription data found!")

    afun <- ifelse(progress, pbapply, apply)

    if(progress) cat("Getting annotations: \n")
    sdrfs <- afun(files, 1, function(fi) {
        path <- file.path(folder, fi["id"])
        mage <- find_dir(path, "mage-tab")
        if (length(mage) == 0) {
            untar(list.files(path=path, pattern="\\.tar\\.gz", full.names=TRUE),
                exdir=path.expand(path))
            mage <- find_dir(path, "mage-tab")
        }
        path <- mage
        sdrf_path <- list.files(path=path, pattern=HUEX_SDRF_RE, full.names=TRUE)
        fread(sdrf_path, na.strings="->", select=HUEX_SDRF_COLS,
            colClasses="character")
    })
    sdrf <- unique(rbindlist(sdrfs), by=NULL)
    names(sdrf) <- HUEX_SDRF_NAMES
    is_tumor <- as.numeric(sapply(sdrf$barcode, substr, 14, 15)) < 10
    sdrf$tumor <- is_tumor
    sdrf$panel <- toupper(str_match(sdrf$panel, PANEL_RE)[,2])
    sdrf[, "barcode" := sapply(barcode, substr, 0, 16)]

    afun <- ifelse(progress, pblapply, lapply)
    if(progress) cat("Reading arrays: \n")
    arrays <- afun(1:nrow(files), function(i) {
        fi <- files[i]
        path <- file.path(folder, fi[, id], fi[, filename])
        cols <- colnames(fread(path, nrows=0, header=T))
        arr <- fread(path, skip=2)
        rows <- arr[[1]]
        set(arr, j=1L, value=NULL)
        arr <- as.matrix(arr)
        rownames(arr) <- rows
        colnames(arr) <- cols[-1]
        arr
    })
    arrays <- do.call("cbind", arrays)
    arrays <- arrays[, sort(colnames(arrays))]
    if (features == "genes") feat <- data.table(symbol=rownames(arrays))
    else feat <- data.table(ID=rownames(arrays))

    matched <- as.numeric(sdrf$id) %in% as.numeric(colnames(arrays))
    ids <- sort(sdrf[matched, id])
    sdrf <- sdrf[matched]
    sdrf <- sdrf[order(id)]
    arrays <- arrays[, ids]

    colnames(arrays) <- sdrf$barcode

    return(list(assay=arrays, samples=sdrf, features=feat))
}
