## code to prepare `Capper_test_data` dataset goes here

library(GEOquery)
gse <- getGEO("GSE109379", getGPL = FALSE)

pheno_df <- pData(gse[[1]])

pheno_df <- pheno_df[pheno_df$`methylation class:ch1`,]


usethis::use_data(Capper_test_data, overwrite = TRUE)
