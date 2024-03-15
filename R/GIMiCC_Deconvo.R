#' @title
#' Function to perform GIMiCC hierarchical deconvolution
#'
#' @description
#' Given DNA methylation microarray data as an input, this function will
#' estimate the relative composition of the sample. At the highest resolution,
#' this will include 18 cell-types. Users may also specify lower resolution
#' estimates as well to observe trends in broader cell populations.
#'
#' @param
#' Betas A matrix (CpG as rows, samples as columns) of DNA methylation data
#' presented as beta values. 450k, EPIC and EPICv2 platforms are accepted.
#'
#' @param
#' h The hierarchical level of resolution desired. Parameter must equal one of
#' the following {0,1,2,3,4,5}. The larger the number, the higher the
#' resolution where level 0 provides estimates for only "Tumor" and "Non-Tumor"
#' fractions whereas level 5 provides estimates for all 18 cell types.
#'
#' @param
#' tumor.type Must be one of the following {"AST-HG","AST-LG","GBM","OLG"} as
#' defined by the Capper et al molecular classifier. AST-HG is for high-grade 
#' astrocytoma, AST-LG is for low-grade astrocytoma, GBM is for glioblastoma,
#' and OLG is for oligodendroglioma.
#' 
#' @return
#' A dataframe of cellular composition estimates per sample, scaled to sum of
#' 100%.
#'
#' @import InfiniumPurify
#' @import dplyr
#' @import minfi
#' @importFrom FlowSorted.Blood.EPIC projectCellType_CP
#'
#' @examples
#' ## Grabbing data
#' library(ExperimentHub)
#' Example_Beta<-query(ExperimentHub(), "GIMiCC")[["EH????????????"]]
#'
#' ## Using GIMiCC
#' GIMiCC_results <- GIMiCC_Deconvo(Example_Beta, tumor.type = "GBM")
#' head(GIMiCC_results)
#'
#' @export

GIMiCC_Deconvo <- function(Betas, h = 5, tumor.type = "") {
  ## testing if input is correct
  if (sum(Betas < 0, na.rm = TRUE) > 0) {
    stop("Negative values detected - please use Betas as input")
  }
  if (sum(Betas > 1, na.rm = TRUE) > 0) {
    stop("Values greater than one detected - please use Betas as input")
  }

  ## Loading the correct libraries
  if (tumor.type == "AST-LG") {
    message("\nLoading libraries for Low-Grade Astrocytoma Tumor Deconvolution")
    GIMiCC_libraries <- data_env[["AST-LG_libraries"]]
  } else if (tumor.type == "AST-HG") {
    message("\nLoading libraries for High-Grade Astrocytoma Tumor Deconvolution")
    GIMiCC_libraries <- data_env[["AST-HG_libraries"]]
  } else if (tumor.type == "GBM") {
    message("\nLoading libraries for Glioblastoma Tumor Deconvolution")
    GIMiCC_libraries <- data_env[["GBM_libraries"]]
  } else if (tumor.type == "OLG") {
    message("\nLoading libraries for Oligodendroglioma Tumor Deconvolution")
    GIMiCC_libraries <- data_env[["OLG_libraries"]]
  } else {
    stop("Invalid tumor type")
  }

  Library_Layer0 <- GIMiCC_libraries[["Library_Layer0"]]
  Library_Layer1 <- GIMiCC_libraries[["Library_Layer1"]]
  Library_Layer2A <- GIMiCC_libraries[["Library_Layer2A"]]
  Library_Layer2B <- GIMiCC_libraries[["Library_Layer2B"]]
  Library_Layer2C <- GIMiCC_libraries[["Library_Layer2C"]]
  Library_Layer2D <- GIMiCC_libraries[["Library_Layer2D"]]
  Library_Layer3A <- GIMiCC_libraries[["Library_Layer3A"]]
  Library_Layer3B <- GIMiCC_libraries[["Library_Layer3B"]]
  Library_Layer4 <- GIMiCC_libraries[["Library_Layer4"]]
  Library_Layer5A <- GIMiCC_libraries[["Library_Layer5A"]]
  Library_Layer5B <- GIMiCC_libraries[["Library_Layer5B"]]
  Library_Layer5C <- GIMiCC_libraries[["Library_Layer5C"]]
  rm(GIMiCC_libraries)

  ## calculating tumor purity
  tumor_sample <- colnames(Betas)
  tumor_beta_iDMC <- Betas[rownames(Betas) %in% rownames(Library_Layer0), ]
  idmc_dat <- Library_Layer0[rownames(tumor_beta_iDMC), ]
  purity <- c()

  for (t in tumor_sample) {
    beta_adj <- c(
      tumor_beta_iDMC[idmc_dat$hyper == TRUE, t],
      1 - tumor_beta_iDMC[idmc_dat$hyper == FALSE, t]
    )
    pu <- InfiniumPurify:::.get_peak(beta_adj)
    purity[t] <- pu
  }
  purity_iDMC <- as.data.frame(purity)
  proj <- purity_iDMC
  proj$NonTumor <- 1 - proj$purity
  colnames(proj)[1] <- "Tumor"
  h0_proj <- proj

  ## projecting each hierarchial layer
  ## Layer 1
  proj1 <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer1[which(rownames(Library_Layer1) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer1[which(rownames(Library_Layer1) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj1[proj1 < 1e-05] <- 0

  message(
    nrow(Library_Layer1) -
      nrow(Library_Layer1[which(rownames(Library_Layer1) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer1),
    "probes were missing in the Betas matrix for L1"
  )

  ## Layer 2
  proj2A <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer2A[which(rownames(Library_Layer2A) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer2A[which(rownames(Library_Layer2A) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj2A[proj2A < 1e-05] <- 0
  for (i in 1:nrow(proj2A)) {
    z <- 1 / sum(proj2A[i, ])
    proj2A[i, ] <- z * proj2A[i, ]
  }
  message(
    nrow(Library_Layer2A) -
      nrow(Library_Layer2A[which(rownames(Library_Layer2A) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer2A),
    "probes were missing in the Betas matrix for L2A"
  )

  proj2B <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer2B[which(rownames(Library_Layer2B) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer2B[which(rownames(Library_Layer2B) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj2B[proj2B < 1e-05] <- 0
  for (i in 1:nrow(proj2B)) {
    z <- 1 / sum(proj2B[i, ])
    proj2B[i, ] <- z * proj2B[i, ]
  }
  message(
    nrow(Library_Layer2B) -
      nrow(Library_Layer2B[which(rownames(Library_Layer2B) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer2B),
    "probes were missing in the Betas matrix for L2B"
  )

  proj2C <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer2C[which(rownames(Library_Layer2C) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer2C[which(rownames(Library_Layer2C) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj2C[proj2C < 1e-05] <- 0
  for (i in 1:nrow(proj2C)) {
    z <- 1 / sum(proj2C[i, ])
    proj2C[i, ] <- z * proj2C[i, ]
  }
  message(
    nrow(Library_Layer2C) -
      nrow(Library_Layer2C[which(rownames(Library_Layer2C) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer2C),
    "probes were missing in the Betas matrix for L2C"
  )

  proj2D <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer2D[which(rownames(Library_Layer2D) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer2D[which(rownames(Library_Layer2D) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj2D[proj2D < 1e-05] <- 0
  for (i in 1:nrow(proj2D)) {
    z <- 1 / sum(proj2D[i, ])
    proj2D[i, ] <- z * proj2D[i, ]
  }
  message(
    nrow(Library_Layer2D) -
      nrow(Library_Layer2D[which(rownames(Library_Layer2D) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer2D),
    "probes were missing in the Betas matrix for L2D"
  )

  ## Layer 3
  proj3A <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer3A[which(rownames(Library_Layer3A) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer3A[which(rownames(Library_Layer3A) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj3A[proj3A < 1e-05] <- 0
  for (i in 1:nrow(proj3A)) {
    z <- 1 / sum(proj3A[i, ])
    proj3A[i, ] <- z * proj3A[i, ]
  }
  message(
    nrow(Library_Layer3A) -
      nrow(Library_Layer3A[which(rownames(Library_Layer3A) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer3A),
    "probes were missing in the Betas matrix for L3A"
  )

  proj3B <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer3B[which(rownames(Library_Layer3B) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer3B[which(rownames(Library_Layer3B) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj3B[proj3B < 1e-05] <- 0
  for (i in 1:nrow(proj3B)) {
    z <- 1 / sum(proj3B[i, ])
    proj3B[i, ] <- z * proj3B[i, ]
  }
  message(
    nrow(Library_Layer3B) -
      nrow(Library_Layer3B[which(rownames(Library_Layer3B) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer3B),
    "probes were missing in the Betas matrix for L3B"
  )

  ## Layer 4
  proj4 <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer4[which(rownames(Library_Layer4) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer4[which(rownames(Library_Layer4) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj4[proj4 < 1e-05] <- 0
  for (i in 1:nrow(proj4)) {
    z <- 1 / sum(proj4[i, ])
    proj4[i, ] <- z * proj4[i, ]
  }
  message(
    nrow(Library_Layer4) -
      nrow(Library_Layer4[which(rownames(Library_Layer4) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer4),
    "probes were missing in the Betas matrix for L4"
  )

  ## Layer 5
  proj5A <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer5A[which(rownames(Library_Layer5A) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer5A[which(rownames(Library_Layer5A) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj5A[proj5A < 1e-05] <- 0
  for (i in 1:nrow(proj5A)) {
    z <- 1 / sum(proj5A[i, ])
    proj5A[i, ] <- z * proj5A[i, ]
  }
  message(
    nrow(Library_Layer5A) -
      nrow(Library_Layer5A[which(rownames(Library_Layer5A) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer5A),
    "probes were missing in the Betas matrix for L5A"
  )

  proj5B <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer5B[which(rownames(Library_Layer5B) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer5B[which(rownames(Library_Layer5B) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj5B[proj5B < 1e-05] <- 0
  for (i in 1:nrow(proj5B)) {
    z <- 1 / sum(proj5B[i, ])
    proj5B[i, ] <- z * proj5B[i, ]
  }
  message(
    nrow(Library_Layer5B) -
      nrow(Library_Layer5B[which(rownames(Library_Layer5B) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer5B),
    "probes were missing in the Betas matrix for L5B"
  )

  proj5C <- as.data.frame(
    projectCellType_CP(
      Betas[rownames(Library_Layer5C[which(rownames(Library_Layer5C) %in%
        rownames(Betas)), ]), ],
      as.matrix(Library_Layer5C[which(rownames(Library_Layer5C) %in%
        rownames(Betas)), ]),
      lessThanOne = TRUE
    )
  )
  proj5C[proj5C < 1e-05] <- 0
  for (i in 1:nrow(proj5C)) {
    z <- 1 / sum(proj5C[i, ])
    proj5C[i, ] <- z * proj5C[i, ]
  }
  message(
    nrow(Library_Layer5C) -
      nrow(Library_Layer5C[which(rownames(Library_Layer5C) %in%
        rownames(Betas)), ]),
    " out of ",
    nrow(Library_Layer5C),
    "probes were missing in the Betas matrix for L5C"
  )

  ## Scaling each layer by the proportions of the previous layer
  proj1 <- proj1 / rowSums(proj1)
  proj <- cbind(
    proj[, c("Tumor")],
    proj[, c("NonTumor")] * proj1
  )
  colnames(proj)[1] <- "Tumor"
  h1_proj <- proj

  proj2A <- proj2A[, -which(colnames(proj2A) %in%
    c("Glial", "Neuronal", "Immune")),
  drop = FALSE
  ]
  proj2A <- proj2A / rowSums(proj2A)
  proj <- cbind(
    proj[, -which(colnames(proj) %in%
      c("Endothelial and Stromal"))],
    proj[, c("Endothelial and Stromal")] * proj2A
  )
  h2A_proj <- proj

  proj2B <- proj2B[, -which(colnames(proj2B) %in%
    c(
      "Endothelial and Stromal",
      "Neuronal", "Immune"
    )),
  drop = FALSE
  ]
  proj2B <- proj2B / rowSums(proj2B)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("Glial"))],
    proj[, c("Glial")] * proj2B
  )
  h2B_proj <- proj

  proj2C <- proj2C[, -which(colnames(proj2C) %in%
    c("Endothelial and Stromal", "Glial", "Immune")),
  drop = FALSE
  ]
  proj2C <- proj2C / rowSums(proj2C)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("Neuronal"))],
    proj[, c("Neuronal")] * proj2C
  )
  h2C_proj <- proj

  proj2D <- proj2D[, -which(colnames(proj2D) %in%
    c(
      "Endothelial and Stromal",
      "Glial", "Neuronal"
    )),
  drop = FALSE
  ]
  proj2D <- proj2D / rowSums(proj2D)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("Immune"))],
    proj[, c("Immune")] * proj2D
  )
  h2D_proj <- proj

  proj3A <- proj3A[, which(colnames(proj3A) %in% c("Tcell", "NK", "Bcell")),
    drop = FALSE
  ]
  proj3A <- proj3A / rowSums(proj3A)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("Lymphoid"))],
    proj[, c("Lymphoid")] * proj3A
  )
  h3A_proj <- proj

  proj3B <- proj3B[, which(colnames(proj3B) %in% c("Mono", "Neu")),
    drop = FALSE
  ]
  proj3B <- proj3B / rowSums(proj3B)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("Myeloid"))],
    proj[, c("Myeloid")] * proj3B
  )
  h3B_proj <- proj

  proj4 <- proj4[, which(colnames(proj4) %in% c("CD4Tcell", "CD8Tcell")),
    drop = FALSE
  ]
  proj4 <- proj4 / rowSums(proj4)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("Tcell"))],
    proj[, c("Tcell")] * proj4
  )
  h4_proj <- proj

  proj5A <- proj5A[, which(colnames(proj5A) %in% c("CD4mem", "CD4nv", "Treg")),
    drop = FALSE
  ]
  proj5A <- proj5A / rowSums(proj5A)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("CD4Tcell"))],
    proj[, c("CD4Tcell")] * proj5A
  )
  h5A_proj <- proj

  proj5B <- proj5B[, which(colnames(proj5B) %in% c("CD8mem", "CD8nv")),
    drop = FALSE
  ]
  proj5B <- proj5B / rowSums(proj5B)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("CD8Tcell"))],
    proj[, c("CD8Tcell")] * proj5B
  )
  h5B_proj <- proj

  proj5C <- proj5C[, which(colnames(proj5C) %in% c("Bmem", "Bnv")),
    drop = FALSE
  ]
  proj5C <- proj5C / rowSums(proj5C)
  proj <- cbind(
    proj[, -which(colnames(proj) %in% c("Bcell"))],
    proj[, c("Bcell")] * proj5C
  )
  h5C_proj <- proj

  ## cleaning the data
  is.nan.data.frame <- function(x) {
    do.call(cbind, lapply(x, is.nan))
  }

  proj[is.nan.data.frame(proj)] <- 0
  proj$Sum <- round(rowSums(proj), 2)
  proj_low <- proj %>% filter(Sum < 1)
  ID_low <- rownames(proj_low)
  proj <- proj[, !colnames(proj) == "Sum"]
  proj[ID_low, ] <- h2C_proj[ID_low, ]

  proj$NonTumor <- h0_proj$NonTumor
  proj$Neuronal <- h1_proj$Neuronal
  proj$Glial <- h1_proj$Glial
  proj$`Endothelial and Stromal` <- h1_proj$`Endothelial and Stromal`
  proj$Immune <- h1_proj$Immune
  proj$Lymphoid <- h2D_proj$Lymphoid
  proj$Myeloid <- h2D_proj$Myeloid
  proj$Tcell <- h3B_proj$Tcell
  proj$Bcell <- h3B_proj$Bcell
  proj$CD8Tcell <- h4_proj$CD8Tcell
  proj$CD4Tcell <- h4_proj$CD4Tcell

  proja <- proj[!rownames(proj) %in% ID_low, ]
  proja[is.nan.data.frame(proja)] <- 0
  proj[rownames(proja), ] <- proja

  # selected the info for the user input layer
  if (h == "0") {
    output <- proj[, c("Tumor", "NonTumor")]
  } else {
    if (h == "1") {
      output <- proj[, c(
        "Tumor", "Immune", "Endothelial and Stromal",
        "Glial", "Neuronal"
      )]
    } else {
      if (h == "2") {
        output <- proj[, c(
          "Tumor", "Lymphoid", "Myeloid", "Endothelial",
          "Stromal", "Astrocyte", "Microglia",
          "Oligodendrocyte", "GABA", "GLU"
        )]
      } else {
        if (h == "3") {
          output <- proj[, c(
            "Tumor", "Tcell", "Bcell", "NK", "Mono",
            "Neu", "Endothelial", "Stromal", "Astrocyte",
            "Microglia", "Oligodendrocyte", "GABA", "GLU"
          )]
        } else {
          if (h == "4") {
            output <- proj[, c(
              "Tumor", "CD4Tcell", "CD8Tcell", "Bcell",
              "NK", "Mono", "Neu", "Endothelial", "Stromal",
              "Astrocyte", "Microglia", "Oligodendrocyte",
              "GABA", "GLU"
            )]
          } else {
            if (h == "5") {
              output <- proj[, c(
                "Tumor", "CD4nv", "CD4mem", "Treg", "CD8nv",
                "CD8mem", "Bnv", "Bmem", "NK", "Mono", "Neu",
                "Endothelial", "Stromal", "Astrocyte",
                "Microglia", "Oligodendrocyte", "GABA",
                "GLU"
              )]
            }
          }
        }
      }
    }
  }



  output$Sum <- rowSums(output)

  output_low <- output %>% filter(Sum == "NaN")
  ID_low <- rownames(output_low)
  if (length(ID_low) != 0) {
    print("NaN indicates noisy deconvolution signal, thus removed")
  }
  return(output[, !colnames(output) == "Sum"] * 100)
}
