#' Test GBM DNAm data used for testing
#'
#' A random subset of the Capper et al. dataset (GSE90496) including 10 
#' glioblastoma (GBM) samples. DNAm data is represented in a data matrix of
#' beta values, or the fraction of methylated alleles at the given CpG.
#'
#' @docType data
#'
#' @format A matrix with 485512 rows (CpGs) and 10 columns (samples).
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90496}
"Capper_test_betas"