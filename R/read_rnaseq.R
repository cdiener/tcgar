#  read_rnaseq.R
#
#  Copyright 2015 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

RNA_FILE_EXT <- c("counts"="htseq\\.counts\\.gz", "FPKM"="FPKM\\.txt\\.gz",
    "FPKM-UQ"="FPKM-UQ\\.txt\\.gz")

# To fix stupid CRAN notes
utils::globalVariables(c("gene_id", "name", "sample_uuid", "V2"))

#' Reads Level 3 RNA sequencing data from GDC.
#'
#' Reads RNASeq Level 3 data downloaded from GDC. The functions does some
#' conservative checking in order to validate that all samples and features are
#' annotated correctly. Also note that duplicate data (files with the same
#' extract name) will be treated as independent duplicates for the same sample.
#'
#' @seealso \code{\link{read_rnaseq_legacy}} to read legacy data.
#' @export
#' @keywords TCGA read RNASeq GDC
#' @param manifest Path to the GDC file manifest.
#' @param folder Folder where the data files reside.
#' @param normalization The normalization method. Must be one of
#'  \describe{
#'      \item{"counts"}{The raw counts.}
#'      \item{"FPKM"}{Fragments Per Kilobase of transcript per Million mapped
#'      reads.}
#'      \item{"FPKM-UQ"}{Fragments Per Kilobase of transcript per Million mapped
#'      reads with Upper Quartile normlization.}
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
#' @importFrom R.utils gunzip
read_rnaseq <- function(manifest, folder, normalization="counts",
    progress=TRUE) {
    man <- fread(manifest)
    afun <- ifelse(progress, pbapply, apply)

    if (normalization %in% names(RNA_FILE_EXT)) {
        files <- man[grep(RNA_FILE_EXT[normalization], filename)]
    } else stop("Not a valid normalization method!")

    files <- files[id %in% tcgar::gdc_files$file_id]
    samples <- tcgar::gdc_files
    setkey(samples, "file_id")
    samples <- samples[files$id]
    setkey(samples, NULL)

    path <- files[1, file.path(folder, id, filename)]
    ext <- sub("\\.gz", "", path)
    if (!file.exists(ext)) gunzip(path)
    feat <- fread(ext, select=1, col.names="ensgene")
    valid_ids <- grep("ENSG", feat$ensgene)
    feat <- data.table(ensgene=sub("\\..+", "", feat$ensgene[valid_ids]))
    feat <- merge(feat, tcgar::genemap, by="ensgene", all.x=TRUE)

    if(progress) {
        cat("Reading RNA-Seq counts:\n")
        pb <- startpb(min=0, max=nrow(files))
        on.exit(closepb(pb))
    }

    counts <- matrix(0, nrow=nrow(feat), ncol=nrow(files))  # pre-allocate
    for(i in 1:nrow(files)) {
        path <- files[i, file.path(folder, id, filename)]
        ext <- sub("\\.gz", "", path)
        if (!file.exists(ext)) gunzip(path)
        if (progress) setpb(pb, i)
        counts[, i] <- fread(ext)[valid_ids, V2]
    }
    dimnames(counts) <- list(feat$ensgene, files$id)

    return(list(counts=counts, features=feat, samples=samples))
}

#' Gene annotations obatined from BioMart.
#'
#' This data set contains additional annotations for genes
#' obtained from BioMart.
#'
#' @format A data frame with 66,203 rows and 4 variables:
#' \describe{
#'   \item{ensgene}{Ensembl Gene ID}
#'   \item{description}{Description of the Gene}
#'   \item{symbol}{The Gene symbol}
#'   \item{entrez}{Entrez Gene ID}
#' }
"genemap"
