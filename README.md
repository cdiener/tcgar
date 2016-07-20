[![wercker status](https://app.wercker.com/status/7d169bc8d8edd03c8ec8c198038698d4/s 
"wercker status")](https://app.wercker.com/project/bykey/7d169bc8d8edd03c8ec8c198038698d4)
[![codecov.io](https://codecov.io/github/cdiener/tcgar/coverage.svg?branch=master)]
(https://codecov.io/github/cdiener/tcgar?branch=master)

# tcgar
R package to read TCGA data and connect it to analysis pipelines (e.g. bioconductor).

## Important!

TCGA has recently migrated all of its data to the [Genomic Data Commons Portal](https://gdc-portal.nci.nih.gov/).
Due to this `tcgar` is currently ***not functional anymore***. I am currently in the process of implementing
the new API so it will work again.

Currently implements download and parsing for:

- Clinical data for samples and patients
- RNASeq2 data (version 2 from level 3 files)
- Exon Expression data (HuEx)
