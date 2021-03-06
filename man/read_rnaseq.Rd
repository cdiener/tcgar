% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_rnaseq.R
\name{read_rnaseq}
\alias{read_rnaseq}
\title{Reads Level 3 RNA sequencing data from GDC.}
\usage{
read_rnaseq(manifest, folder, normalization = "counts", progress = TRUE)
}
\arguments{
\item{manifest}{Path to the GDC file manifest.}

\item{folder}{Folder where the data files reside.}

\item{normalization}{The normalization method. Must be one of
\describe{
    \item{"counts"}{The raw counts.}
    \item{"FPKM"}{Fragments Per Kilobase of transcript per Million mapped
    reads.}
    \item{"FPKM-UQ"}{Fragments Per Kilobase of transcript per Million mapped
    reads with Upper Quartile normlization.}
}}

\item{progress}{Logical. Show progress info?}
}
\value{
A data table containing the features as rows and the samples in the
 columns.
}
\description{
Reads RNASeq Level 3 data downloaded from GDC. The functions does some
conservative checking in order to validate that all samples and features are
annotated correctly. Also note that duplicate data (files with the same
extract name) will be treated as independent duplicates for the same sample.
}
\examples{
# Not run due to large download...
# gbm <- system.file("extdata", "manifest.tsv", package = "tcgar")
# d <- tempdir()
# rna <- read_rnaseq(gbm, d)

}
\seealso{
\code{\link{read_rnaseq_legacy}} to read legacy data.
}
\keyword{GDC}
\keyword{RNASeq}
\keyword{TCGA}
\keyword{read}
