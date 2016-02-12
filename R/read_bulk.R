#  read_bulk.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

rcomb <- function(x,y) rbindlist(list(x, y), fill=TRUE)
fcomb <- function(x,y) {
    if (is.null(x) || all(dim(x) == dim(y))) y
    else stop("Incompatible feature annotations!") 
}

manager <- list(RNASeqV2 = list(fun="read_rnaseq", 
    summarize=c(features=fcomb, samples=rcomb, counts=cbind)),
    HuEx = list(fun="read_huex", summarize=c(features=fcomb, samples=rcomb, 
    assay=cbind)),
    clinical = list(fun="read_clinical", summarize=c(samples=rcomb, patients=rcomb)))

efun <- function(er) NULL

#' Reads several data types from several TCGA download directories.
#'
#' \code{read_bulk} allows you to combine several data types across several
#' directories containing TCGA data. Data from several directories will be 
#' combined correctly automatically.
#' 
#' @export
#' @keywords TCGA read bulk
#' @param folders A character vector of folders containing downloaded TCGA data.
#' @param what What technologies to read. Must be a character vector consisting 
#'  of any of the following: "HuEx", "RNASeqV2", "clinical".
#' @param args A named list containing additional arguments to the function. Only
#'  the \code{folder} argument will be filled in by the bulk reader.
#' @param stats Logical. Whether to print some basic stats about the bulk read.
#' @return A list of length(what) containing the combined data or NULL if no data
#'  was found.
#' @examples
#' gbm <- system.file("extdata", "GBM", package = "tcgar")
#' rna <- read_bulk(gbm)
#'
#' @importFrom data.table rbindlist
read_bulk <- function(folders, what=c("HuEx", "RNASeqV2", "clinical"), 
    args=list(), stats=TRUE) {
    out <- list()
    for(te in what) out[[te]] <- lapply(manager[[te]]$summarize, function(x) NULL)
    samples_per <- matrix(0, nrow=length(folders), ncol=length(what), 
        dimnames=list(folders, what))
    for (fi in 1:length(folders)) {
        f <- folders[fi]
        cat("                                                        \r")
        cat(sprintf("Reading folder %d/%d...", fi, length(folders)))
        for (tech in what) {
            args[[tech]]$folder <- f 
            dat <- tryCatch(do.call(manager[[tech]]$fun, args[[tech]]), error=efun)
            if (!is.null(dat)) samples_per[f, tech] <- nrow(dat$samples)
            for (field in names(manager[[tech]]$summarize)) {
                sfun <- manager[[tech]]$summarize[[field]]
                out[[c(tech, field)]] <- sfun(out[[c(tech, field)]], dat[[field]])
            }
        }
    }
    cat("\n")
    
    out$summary <- samples_per
    
    if (stats) {
        write("Sample distribution:", file="")
        print(samples_per)
        cat("\n")
        size <- object.size(out)
        write(paste("Size of data set in RAM =", format(size, units="auto")), 
            file="")
    }
    
    if (sum(samples_per) == 0) return(NULL)
    return(out)
}
