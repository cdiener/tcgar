#  get_data.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.


#' Downloads data from the Genomic Data Commons Portal.
#'
#' This will download all data specified in a given manifest file. To generate
#' a valid manifest use the portal at \url{https://gdc-portal.nci.nih.gov/}.
#'
#' @seealso \code{\link{read_bulk}} to read data from several folders and
#'  technologies.
#' @export
#' @keywords download, TCGA, GDC
#' @param manifest Path to the GDC file manifest.
#' @param out Optional. The folder to which to download the data. Default
#'  is to extract into a folder named "gdc_data".
#' @param quiet Specifies whether download status info should be supressed (default no).
#' @param np How many parallel download processes should be started.
#' @return Nothing.
#' @examples
#' d <- tempdir()
#' get_data("GBM", out=d)
#'
#' @importFrom utils untar
get_data <- function(manifest, out="gdc_data", quiet=FALSE, np=4) {
    suppressWarnings(dir.create(out, recursive=TRUE))

    man <- fread(manifest)
    bytes <- structure(man[, sum(as.numeric(size))], class="object_size")
    cat(sprintf("Downloading %s...\n", format(bytes, "auto")))
    sout <- ifelse(quiet, FALSE, "")
    command <- ifelse(Sys.info()["sysname"] == "Windows",
        "gdc-client.exe", "gdc-client")
    args <- c("download", "-m", manifest, "-d", out, "-n", as.character(np))
    system2(command, args=args, stdout=sout, stderr=sout)
}
