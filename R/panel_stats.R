#  panel_stats.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

utils::globalVariables("size")

#' Counts samples per panel across several data sets.
#'
#' @export
#' @keywords TCGA statistics
#' @param ... Several or one data set as returned by the `read_*` functions.
#' @return A list of length(what) containing the combined data or NULL if no data
#'  was found.
#' @examples
#' gbm <- system.file("extdata", "manifest.tsv", package = "tcgar")
#' d <- tempdir()
#' panel_samples(huex=read_huex(gbm, d, progress=FALSE))
#'
#' @importFrom data.table rbindlist
#' @importFrom utils object.size
panel_samples <- function(...) {
    data <- list(...)
    counts <- lapply(data, function(d) {
        if ("list" %in% class(d)) d <- d$samples
        tab <- table(d$panel)
        data.frame(panel=names(tab), counts=as.numeric(tab))
    })
    out <- Reduce(function(df1, df2) merge(df1, df2, by="panel"), counts)
    nam <- names(data)
    if (is.null(nam)) nam <- paste0("exp_", seq_along(data))
    names(out) <- c("panel", nam)
    out <- out[order(out$panel), ]

    return(out)
}
