[![wercker status](https://app.wercker.com/status/7d169bc8d8edd03c8ec8c198038698d4/s
"wercker status")](https://app.wercker.com/project/bykey/7d169bc8d8edd03c8ec8c198038698d4)
[![codecov.io](https://codecov.io/github/cdiener/tcgar/coverage.svg?branch=master)]
(https://codecov.io/github/cdiener/tcgar?branch=master)

# tcgar
R package to read TCGA data and connect it to analysis pipelines (e.g. bioconductor).

## Important!

TCGA has recently migrated all of its data to the [Genomic Data Commons Portal](https://gdc-portal.nci.nih.gov/).
I am currently in the process of implementing the new API so it will work again. This here is the new
GDC branch that can use GDC data.

## Usage

Get and install the [GDC Data Transfer Tool](https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool).
You can download data from within tcgar or by hand. `tcgar` uses a download manifest to download
large selections of data at once. To obtain the manifest go to the [GDC data portal](https://gdc-portal.nci.nih.gov/search/s)
or the [legacy portal](https://gdc-portal.nci.nih.gov/legacy-archive/search/f) (to get the original TCGA data).
This will present tou with a nice UI to select the data set you want. After that click
on "Download Manifest" and save it to your disk. This is all you will need to download
data.

Currently implements download and parsing for:

- Clinical data for samples and patients
- RNASeq data (from level 3 files)
- Exon Expression data (HuEx)
