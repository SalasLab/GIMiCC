### =========================================================================
### GIMiCC metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c("Capper_example_betas.rda","GIMiCC_Library.rda"),
  Description = c(paste0("The Capper_example_betas.rda ",
                         "contains a matrix of beta values from three ",
                         "glioblastoma samples GSM2403088, ",
                         "GSM2403089, and GSM2403095. ",
                         "These samples are from publicly available, ",
                         "data source GEO accession number GSE109381, ",
                         "at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109381 "),
                  paste0("The GIMiCC_Library.rda ",
                         "contains matrices of the the average ",
                         "DNA methylation values of the probes included ",
                         "in 6 layers of the GIMiCC ",
                         "deconvolution for 4 molecular subtypes of glioma. "
                         )),
  BiocVersion = c("3.17"), # check
  Genome = rep("hg19", 1),
  SourceType = rep("tar.gz", 1),
  SourceUrl = c("https://bit.ly/42HtPK9",paste0("https://bit.ly/42HtPK9, ", "https://bit.ly/3TZpTkl, ",
                     "https://bit.ly/3lBBHMR, ", "https://bit.ly/3z9P8GI, ",
                     "https://bit.ly/40hgaYO"
                     )), # check
  SourceVersion = "Mar 15 2024",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = paste0("GEO"),
  Maintainer = "", # check
  RDataClass = c("Beta Matrix") ,
  DispatchClass = c(rep("Rda",1)),
  RDataPath = c(paste0("GIMiCC/",
                       "Capper_example_betas.rda"),paste0("GIMiCC/",
                                                  "GIMiCC_Library.rda")),
  Tags = "",
  Notes = paste0("Zhang Z et al 2023, ", "Zhang Z et al 2022, ",
                 "Farkas SA et al 2013, ", "Zhang W et al 2020, ",
                 "Timp W et al 2014, ", "Lennard K et al 2016") #check
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)

