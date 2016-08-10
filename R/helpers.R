#  helpers.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

ADD_NS <- c("shared_stage", "nte")

highest_version <- function(vers) {
    return(order(numeric_version(vers), decreasing=T)[1])
}

find_dir <- function(path, re) {
    dirs <- list.dirs(path)
    grep(re, dirs, value=TRUE)
}

last_or_na <- function(x) {
    if (length(x) == 0) return(NA)
    val <- x[order(x, decreasing=TRUE)[1]]
    if (nchar(val) == 0) return(NA) else return(val)
}

fix_ns <- function(ns) {
    for(x in ADD_NS) if(is.na(ns[x])) ns[x] = "dummy"
    return(ns)
}
