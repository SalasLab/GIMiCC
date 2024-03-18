### =========================================================================
### GIMiCC metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c("Capper_example_betas.rda", "GIMiCC_Library.rda"),
  Description = c(
    paste0(
      "The Capper_example_betas.rda ",
      "contains a matrix of beta values from three ",
      "glioblastoma samples GSM2403088, ",
      "GSM2403089, and GSM2403095. ",
      "These samples are from publicly available, ",
      "data source GEO accession number GSE109381, ",
      "at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109381 "
    ),
    paste0(
      "The GIMiCC_Library.rda ",
      "contains matrices of the the average ",
      "DNA methylation values of the probes included ",
      "in 6 layers of the GIMiCC ",
      "deconvolution for 4 molecular subtypes of glioma. "
    )
  ),
  BiocVersion = c("3.19"), # check
  Genome = rep("hg19", 1),
  SourceType = rep("tar.gz", 1),
  SourceUrl = c("https://bit.ly/4al2O27", paste0(
    "https://bit.ly/4al2O27, ",
    "https://bit.ly/42HtPK9, ",
    "https://bit.ly/4aBpoUH, ",
    "https://bit.ly/3ThjgK2, ",
    "https://bit.ly/3JygxbC, ",
    "https://bit.ly/3ZYMtuX, ",
    "https://bit.ly/3mMFB5P, ",
    "https://bit.ly/429rSpB, ",
    "https://bit.ly/3vhZFS0"
  )), # check
  SourceVersion = "Mar 15 2024",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = paste0(
    "Synapse, ", "Bioconductor, ", "GEO, ",
    "ArrayExpress, ", "TCGA"
  ),
  Maintainer = "Steven C. Pike <steven.c.pike.gr@dartmouth.edu>", # check
  RDataClass = c("Beta Matrix"),
  DispatchClass = c(rep("Rda", 1)),
  RDataPath = c(paste0(
    "GIMiCC/",
    "Capper_example_betas.rda"
  ), paste0(
    "GIMiCC/",
    "GIMiCC_Library.rda"
  )),
  Tags = "",
  Notes = paste0(
    "Capper D et al 2018, ", "Zhang Z et al 2023, ", "Zhang Z et al 2022, ",
    "Weightman Potter PG et al 2021, ", "de Witte et al 2022, ",
    "Mendizabal et al 2015, ", "Kozlenkov et al 2018, ", "Salas et al 2022"
  ) # check
)

write.csv(meta, file = "inst/extdata/metadata.csv", row.names = FALSE)
