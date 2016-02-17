#  get_panels.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

TCGA_FTP <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"

FLRE <- "<a href=.+>(.+)</a>\\s+\\d+-\\d+-\\d+"
highest_version <- function(vers) {    
    return(order(numeric_version(vers), decreasing=T)[1])
}

dir_list <- function(site) {
    if (url.exists(site)) {
        cont <- getURLContent(site)
        ma <- str_match_all(cont, FLRE)[[1]]
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
    panels <- str_match(dirs, "^([a-z_]{2,})/$")
    panels <- panels[!is.na(panels[,2]),2]
    return(toupper(panels))
}

