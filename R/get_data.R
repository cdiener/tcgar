#  get_data.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

TECHS <- c("RNASeqV2", "HuEx", "clinical")
CLIN <- c("biospecimen_sample_\\w+\\.txt", "clinical_patient_\\w+\\.txt")
RNA <- "unc\\.edu_(\\w+)\\.(\\w+)_RNASeqV2\\.Level_3\\.(\\d+\\.\\d+\\.\\d+)\\.tar\\.gz"
HUEX <- "lbl\\.gov_(\\w+)\\.HuEx-1_0-st-v2\\.Level_3\\.(\\d+\\.\\d+\\.\\d+)\\.tar\\.gz"
MT <- "^.+\\.mage-tab\\.(\\d+\\.\\d+\\.\\d+)/"

download_sdrf <- function(files, base, out, method, quiet) {
    mage_tab <- str_match(files, MT)
    mage_tab <- mage_tab[!is.na(mage_tab[,1]), , drop=F]
    mage_tab <- mage_tab[highest_version(mage_tab[,2]), 1]
    sdrf <- grep("sdrf", dir_list(paste0(base, "/", mage_tab)), value=TRUE)

    download.file(paste0(base, "/", mage_tab, sdrf), file.path(out, sdrf),
        method=method, quiet=quiet)

    return(sdrf)
}

download_clin <- function(base, out, method="auto", quiet=TRUE) {
    files <- dir_list(base)
    files <- sapply(CLIN, grep, x=files, value=TRUE)

    if (any(sapply(files, length) == 0)) return(NULL)

    sapply(files, function(f) download.file(paste0(base, f),
        file.path(out,f), method=method, quiet=quiet))

    invisible(files)
}

download_rna <- function(base, out, method="auto", quiet=TRUE) {
    files <- dir_list(base)
    rna_files <- str_match(files, RNA)
    rna_files <- rna_files[!is.na(rna_files[,1]), ]

    if (length(rna_files) == 0) return(NULL)
    rna_file <- rna_files[highest_version(rna_files[,4]), ]

    if (!file.exists(file.path(out, rna_file[1])))
        download.file(paste0(base, rna_file[1]), file.path(out, rna_file[1]),
            method=method, quiet=quiet)

    sdrf <- download_sdrf(files, base, out, method, quiet)

    invisible(c(rna_file, sdrf))
}

download_huex <- function(base, out, method="auto", quiet=TRUE) {
    files <- dir_list(base)
    huex_files <- str_match(files, HUEX)
    huex_files <- huex_files[!is.na(huex_files[,1]), ]

    if (length(huex_files) == 0) return(NULL)

    sapply(huex_files[,1], function(f) if (!file.exists(file.path(out, f)))
        download.file(paste0(base, f), file.path(out, f), method=method,
            quiet=quiet))

    sdrf <- download_sdrf(files, base, out, method, quiet)

    invisible(c(huex_files, sdrf))
}

#' Downloads data for a given cancer panel and technology.
#'
#' This will download all data for a given cancer panel and technology as
#' a given cancer type and extract it.
#'
#' \emph{Note that you should not extract several downloads of the same technology
#' to the same folder as this will be incompatible with the \code{read_*}
#' functions.}
#'
#' @seealso \code{\link{read_bulk}} to read data from several folders and
#'  technologies.
#' @export
#' @keywords download, TCGA
#' @param panel Name of the cancer panel. See \code{\link{get_panels}} to get
#'  a list of availabe cancer panels.
#' @param tech The technology for which data should be downloaded. Currently only
#'  supports 'HuEx', 'RNASeqV2' and 'clinical'. Can be a single technology or
#'  a character vector of technologies.
#' @param out Optional. The folder to which to extract the data. Default
#'  is to extract into a folder named after the cancer panel within the current
#'  working directory.
#' @param extract Boolean. Whether to extract the data after download.
#' @param method Passed to \code{\link{download.file}}. Specifies the download
#'  method. Defaults to 'auto'.
#' @param quiet Passed to \code{\link{download.file}}. Specifies whether download
#'  status info should be supressed. Defaults to false.
#' @return Nothing.
#' @examples
#' d <- tempdir()
#' get_data("GBM", out=d)
#'
#' @importFrom utils download.file
get_data <- function(panel, tech="clinical", out=toupper(panel), extract=TRUE,
    method="auto", quiet=FALSE) {
    suppressWarnings(dir.create(out, recursive=TRUE))

    panel <- tolower(panel)
    base <- paste0(TCGA_FTP, panel)
    dl_cgcc <- dir_list(paste0(base, "/cgcc/"))

    if ("clinical" %in% tech) {
        download_clin(paste0(base, "/bcr/biotab/clin/"), out, method=method,
            quiet=quiet)
    }
    if ("RNASeqV2" %in% tech && any(grepl("unc\\.edu", dl_cgcc))) {
        rna_folders <- grep("_rnaseqv2", dir_list(paste0(base, "/cgcc/unc.edu/")),
            value=TRUE)

        for (f in rna_folders) {
            download_rna(paste0(base, "/cgcc/unc.edu/", f, "rnaseqv2/"), out,
                method=method, quiet=quiet)
        }
    }
    if ("HuEx" %in% tech && any(grepl("lbl\\.gov", dl_cgcc))) {
        download_huex(paste0(base, "/cgcc/lbl.gov/huex-1_0-st-v2/exon/"), out,
            method=method, quiet=quiet)
    }

    if (extract) {
        tars <- list.files(out, pattern="\\.tar\\.gz", full.names=TRUE)
        sapply(tars, untar, compressed=TRUE, exdir=out)
        dummy <- file.remove(tars)
    }
}
