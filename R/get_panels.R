#  get_panels.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

TCGA_FTP <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"

FLRE <- "<a href=.+>(.+)</a>\\d-\\d-\\d"
order_versions <- function(vers) {
    vers <- strsplit(vers, "\\.")
    n <- max(sapply(vers, length))
    vers <- lapply(vers, function(v){ length(v) <- n; v })
    vers <- data.frame(t(data.frame(vers)))
    
    return(do.call(order, vers, decreasing=T))
}

dir_list <- function(url) {
    if (url.exists(url)) {
        cont <- getURLContent(url)
        ma <- str_match_all(url, FLRE)[[1]]
        return(ma[, 2])
    } else NULL
}

#' Lists all available cancer panels in the public TCGA repository.
#' 
#' \emph{Note that you should not extract several downloads of the same technology 
#' to the same folder as this will be incompatible with the \code{read_*} 
#' functions.}
#' 
#' @seealso \code{\link{get}} to download data for a selected cancer panel.
#' @export
#' @return A character vector containing the available cancer panels.
#' @examples
#' gbm <- system.file("extdata", "GBM", package = "tcgar")
#' rna <- read_rnaseq(gbm)
#'
#' @importFrom stringr str_match_all
#' @importFrom RCurl getURLContent url.exists
get_panels <- function() {
    dirs <- dir_list(TCGA_FTP)
    panels <- dirs[grep("^[a-z_]{2,}$", dirs)]
    return(toupper(panels))
}

