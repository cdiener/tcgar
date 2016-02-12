#  get_data.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

TECHS <- c("RNASeqV2", "HuEx", "clinical")

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
#' @keywords panel Name of the cancer panel. See \code{\link{get_list}} to get
#'  a list of availabe cancer panels.
#' @param tech The technology for which data should be downloaded. Currently only
#'  supports 'HuEx', 'RNASeqV2' and 'clinical'.
#' @param out Optional. The folder to which to extract the data. Default (NULL)
#'  is to extract into a folder named after the cancer panel within the current
#'  working directory.
#' @return Nothing.
#' @examples
#' gbm <- system.file("extdata", "GBM", package = "tcgar")
#' rna <- read_rnaseq(gbm)
get_data <- function(panel, tech="clinical", out=NULL) {
    
}

