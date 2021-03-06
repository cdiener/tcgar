% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_file_info.R
\name{list_files}
\alias{list_files}
\title{Obtains infromation about the files in GDC.}
\usage{
list_files(query = "default", unique = TRUE, chunk = 10000,
  max_size = Inf)
}
\arguments{
\item{query}{Optional. The GDC query to be used. By default downloads
detailed information for all files in the GDC.}

\item{unique}{optional boolean. Whether only to return file info
for file that can be mapped to a unique patient and sample (if applicable).}

\item{chunk}{optional number. File info is requested in chunks of this size.
Avoids server errors when requesting all files at once.}

\item{max_size}{Maximum number of entries to download.}
}
\value{
A data table mapping files to various ids and the the patient
}
\description{
This function downloads a list of all files in GDC together with additional
metadata. Here the most important metadata is the mapping of individual
experimental data to the corresponding patient.
}
\examples{
 NULL

}
