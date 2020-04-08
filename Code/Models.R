##################################################
# Model the relationship between CNA, methylation and gene expression
##################################################

## Set the working directory
setwd("/open/tmp/Christian/GeneExpressionModulation/")

##################################################
## Load/install libraries
if (!requireNamespace("doMC", quietly = TRUE))
  install.packages("doMC")
library(doMC)

if (!requireNamespace("mgcv", quietly = TRUE))
  install.packages("mgcv")
library(mgcv)

##################################################
## Parameters

# Number of principal components to use
numPCs <- 5

# Include genes with less than numPCs principal components?
includeFewPCs <- FALSE

# Max and min value for CNA (to add jitter)
cna_cap <- 3.657
cna_cap_low <- -1.293

# Number of base pairs around gene (window) for which a CpG should be classified as "belonging" to a gene
windowBP <- 50000
w_BP_short <- windowBP/1000

# Cores for running in parallel
crs <- 23
registerDoMC(crs)
print(paste0("Cores: ", crs))

# Sample size for models (and cutoff for number of tumors to use)
model_sampleSize <- 100

# Number of repeated runs
runs <- 100

# K value for GAMs (number of basis functions)
gam_k <- 4

# Constant to scale the log-transformed methylation values by
logConst <- 10

##################################################
## Functions

# Normalize between 0 and 1
normalize <- function(x){
  return((x - min(x))/(max(x) - min(x)))
} 

##################################################
## Load data

# Read the sample overview
patientData <- read.table(file = "./Output/SampleOverview.txt",
                          header = TRUE,
                          stringsAsFactors = FALSE)

# Read the gene expression data
exprs <- read.table(file = "./Output/GeneExpressionFormatted_sdFiltered.txt",
                    header = TRUE,
                    stringsAsFactors = FALSE)


# Read the copy number data
copyNumber <- read.table(file = "./Data/ISAR_GISTIC.all_data_by_genes.txt",
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         sep = "\t")

##################################################
## Format data

## Gene expression
# Change periods to dashes in sample names
colnames(exprs) <- gsub("\\.", "-", colnames(exprs))

## Copy number
# Remove first 3 columns, set row names to Gene.Symbol
copyNumber <- data.frame(copyNumber[, 4:ncol(copyNumber)], row.names = copyNumber$Gene.Symbol)

# Change periods to dashes in sample names
colnames(copyNumber) <- gsub("\\.", "-", colnames(copyNumber))

## Sample data
# Remove tumor types with less than n = model_sampleSize tumors
# Get tumor types
tumorTypes <- unique(patientData$TumorType)[order(unique(patientData$TumorType))]

# Get sample numbers for each tumor type
sn <- c()
for(i in tumorTypes){
  sn_temp <- nrow(patientData[patientData$inAll == 1 & patientData$TumorType == i, ])
  sn <- c(sn, sn_temp)
}

# Assign names to the vector 
names(sn) <- tumorTypes

# Get the tumor types with sufficient samples
tumorTypes <- names(sn)[sn > model_sampleSize]

##################################################
## Models

# Find existing files from previous runs
existingFiles <- list.files("./Output/TempFiles/", pattern = "AllModels_")

# For tt = each tumor type 
foreach(tt = tumorTypes) %dopar% {
  
  print(paste("Start", tt))
  
  # Read the methylation data for the given tumor type
  currentMethylation <- readRDS(file = paste0("./Output/MethylationAllPCValues_", w_BP_short, "kWindow_", tt, ".rds"))
  
  # Read the cumulative variances for the given tumor type
  currentVariance <- readRDS(file = paste0("./Output/MethylationAllCumulativeVariances_", w_BP_short, "kWindow_", tt, ".rds"))
  
  ## Find genes for which there are no probes within the chosen window
  # Create an empty list to place NA values in
  nas <- c()
  # For j = each gene in the list of methylation data (j = index in list, not gene name)
  for(j in 1:length(currentMethylation)){
    # If there are any NA values for the given gene, place the gene/index in the list of NA genes
    if(anyNA(currentMethylation[j])){nas <- c(nas, j)}
  }
  
  # Remove the NA genes in the methylation data 
  currentMethylation <- currentMethylation[-nas]
  
  # Remove the NA genes in the PCA variance data 
  currentVariance <- currentVariance[-nas] 
  
  # Find the genes which are found in copy number, gene expression and metyhlation data
  overlappingGenes <- rownames(exprs)[rownames(exprs) %in% names(currentMethylation)]
  overlappingGenes <- overlappingGenes[overlappingGenes %in% rownames(copyNumber)]
  
  for(run in 1:runs){
    
    # Go to next run if previously carried out
    if(paste0("AllModels_", w_BP_short, "kWindow_", tt, "_run", run, ".txt") %in% existingFiles){next}
    
    # Skip problematic runs due to error: "Supplied matrix not symmetric" (https://github.com/cran/mgcv/blob/master/R/mgcv.r)
    if(tt == "CESC" & run == 71){run <- run + runs}
    if(tt == "KIRP" & run == 77){run <- run + runs}
    
    # Start time
    a <- Sys.time()
    
    # Set seed
    set.seed(run)
    
    # Choose samples
    currentSamplesNums <- sample(x = 1:ncol(currentMethylation[[1]]),
                                 size = model_sampleSize)
    currentSamplesNums <- currentSamplesNums[order(currentSamplesNums)]
    
    # Create an empty data frame in which all necessary data can be placed
    allDF_cols <- c("gene",
                    #Variance
                    "expression_variance",
                    "cna_variance",
                    "meth_mean_variance",
                    "meth_variance_captured",
                    "meth_pcs_used",
                    # CNA linear models
                    "lm_logE_C_r2",
                    "lm_logE_C_r2_adj",
                    "lm_logE_C_aic",
                    "lm_logE_logC_r2",
                    "lm_logE_logC_r2_adj",
                    "lm_logE_logC_aic",
                    # CNA GAMs
                    "gam_logE_C_r2",
                    "gam_logE_C_r2_adj",
                    "gam_logE_C_aic",
                    "gam_logE_logC_r2",
                    "gam_logE_logC_r2_adj",
                    "gam_logE_logC_aic",
                    # Methylation linear models
                    "lm_logE_M_r2",
                    "lm_logE_M_r2_adj",
                    "lm_logE_M_aic",
                    "lm_logE_logM_r2",
                    "lm_logE_logM_r2_adj",
                    "lm_logE_logM_aic",
                    # Methylation GAMs
                    "gam_logE_M_r2",
                    "gam_logE_M_r2_adj",
                    "gam_logE_M_aic",
                    "gam_logE_logM_r2",
                    "gam_logE_logM_r2_adj",
                    "gam_logE_logM_aic",
                    # Combined linear models
                    "lm_logE_C_M_r2",
                    "lm_logE_C_M_r2_adj",
                    "lm_logE_C_M_aic",
                    "lm_logE_logC_M_r2",
                    "lm_logE_logC_M_r2_adj",
                    "lm_logE_logC_M_aic",
                    "lm_logE_C_logM_r2",
                    "lm_logE_C_logM_r2_adj",
                    "lm_logE_C_logM_aic",
                    "lm_logE_logC_logM_r2",
                    "lm_logE_logC_logM_r2_adj",
                    "lm_logE_logC_logM_aic",
                    # Both GAM/spline
                    "gam_logE_C_M_r2",
                    "gam_logE_C_M_r2_adj",
                    "gam_logE_C_M_aic",
                    "gam_logE_logC_M_r2",
                    "gam_logE_logC_M_r2_adj",
                    "gam_logE_logC_M_aic",
                    "gam_logE_C_logM_r2",
                    "gam_logE_C_logM_r2_adj",
                    "gam_logE_C_logM_aic",
                    "gam_logE_logC_logM_r2",
                    "gam_logE_logC_logM_r2_adj",
                    "gam_logE_logC_logM_aic",
                    # Methylation spline, CNA linear
                    "methGAM_cnaLM_logE_C_M_r2",
                    "methGAM_cnaLM_logE_C_M_r2_adj",
                    "methGAM_cnaLM_logE_C_M_aic",
                    "methGAM_cnaLM_logE_logC_M_r2",
                    "methGAM_cnaLM_logE_logC_M_r2_adj",
                    "methGAM_cnaLM_logE_logC_M_aic",
                    "methGAM_cnaLM_logE_C_logM_r2",
                    "methGAM_cnaLM_logE_C_logM_r2_adj",
                    "methGAM_cnaLM_logE_C_logM_aic",
                    "methGAM_cnaLM_logE_logC_logM_r2",
                    "methGAM_cnaLM_logE_logC_logM_r2_adj",
                    "methGAM_cnaLM_logE_logC_logM_aic",
                    # CNA spline, methylation linear
                    "cnaGAM_methLM_logE_C_M_r2",
                    "cnaGAM_methLM_logE_C_M_r2_adj",
                    "cnaGAM_methLM_logE_C_M_aic",
                    "cnaGAM_methLM_logE_logC_M_r2",
                    "cnaGAM_methLM_logE_logC_M_r2_adj",
                    "cnaGAM_methLM_logE_logC_M_aic",
                    "cnaGAM_methLM_logE_C_logM_r2",
                    "cnaGAM_methLM_logE_C_logM_r2_adj",
                    "cnaGAM_methLM_logE_C_logM_aic",
                    "cnaGAM_methLM_logE_logC_logM_r2",
                    "cnaGAM_methLM_logE_logC_logM_r2_adj",
                    "cnaGAM_methLM_logE_logC_logM_aic")

    allDF <- data.frame(matrix(data = NA,
                               ncol = length(allDF_cols),
                               nrow = length(overlappingGenes),
                               dimnames = list(overlappingGenes, allDF_cols)))
    
    # Enter gene name
    allDF$gene <- overlappingGenes
    
    # For j = each gene
    for(j in 1:length(overlappingGenes)){
      
      # Monitor how quickly models are running
      print(paste0(tt, ", run: ", run, ", gene: ", j))
      
      # Get the name of the current gene
      currentGene <- overlappingGenes[j]
      
      ## Methylation
      # Subset the methylation data in a temporary data frame
      subDF_meth <- data.frame(t(currentMethylation[[currentGene]]))
      
      # Subset sample
      coln <- colnames(subDF_meth)
      subDF_meth <- data.frame(subDF_meth[currentSamplesNums, ], row.names = rownames(subDF_meth)[currentSamplesNums])
      colnames(subDF_meth) <- coln
      
      # Set the row names for the methylation data to the sample names (+ change periods to dashes in the names)
      rownames(subDF_meth) <- substr(gsub("\\.", "-", rownames(subDF_meth)), start = 1, stop = 16)
      
      ## Expression
      # Get the gene expression sample IDs for the relevant data
      exprIDs <- patientData$ExprID[patientData$SampleID %in% rownames(subDF_meth)]
      
      # Create a data frame with the gene expression data, set row names/sample IDs to match methylation data
      subDF_expr <- matrix((t(exprs[currentGene, exprIDs])),
                           ncol = 1,
                           dimnames = list(substr(exprIDs, start = 1, stop = 16), "GeneExpression"))
      
      ## CNA
      # Get copy number IDs for the relevant data
      cna_ID <- patientData$CNA_ID[patientData$SampleID %in% rownames(subDF_meth)]
      
      # Create a data frame with the copy number data, set row names/sample IDs to match methylation data
      subDF_cna <- matrix((t(copyNumber[currentGene, cna_ID])),
                           ncol = 1,
                           dimnames = list(substr(cna_ID, start = 1, stop = 16), "CNA"))
      
      # Add jitter to max CNA value
      if(max(subDF_cna == cna_cap)){
      subDF_cna[subDF_cna == cna_cap] <- jitter(x = rep(cna_cap, times = length(subDF_cna[subDF_cna == cna_cap])))
      }
      
      # Add jitter to min CNA value
      if(max(subDF_cna == cna_cap_low)){
        subDF_cna[subDF_cna == cna_cap_low] <- jitter(x = rep(cna_cap_low, times = length(subDF_cna[subDF_cna == cna_cap_low])))
      }
      
      # Some cancer types are missing expression data for certain genes; skip gene if only NA values
      if(all(is.na(subDF_expr[, "GeneExpression"]))){next}
      
      # Some genes have identical expression values for all (or almost all) genes; skip if more than 80% of all gene expression values are duplicated
      if(length(unique(subDF_expr[, "GeneExpression"][complete.cases(subDF_expr[, "GeneExpression"])])) < (0.2 * length(subDF_expr[, "GeneExpression"][complete.cases(subDF_expr[, "GeneExpression"])]))){next}
      
      # If there are more principal components than "required" (defined by numPCs)
      if(ncol(subDF_meth) >= numPCs){
        # Use principal components 1 to numPCs (defined earlier) in the analyses
        methCols <- 1:numPCs
      }
      
      # If there are fewer principal components than "required" (defined by numPCs)
      if(ncol(subDF_meth) < numPCs){
        ## Should genes with fewer than numPCs principal components be included?
        # FALSE = No
        if(includeFewPCs == FALSE){
          # Leave blank/NA in allDF, go to next
          next
        }
        # TRUE = Yes
        if(includeFewPCs == TRUE){
          # Use the principal components that are available
          methCols <- 1:ncol(subDF_meth)
        }
      }
      
      # Enter mean methylation variance into the output data frame (prior to normalizing)
      methVars <- c()
      for(k in methCols){
          methVars <- c(methVars, var(subDF_meth[, k]))
        }
      allDF[currentGene, "meth_mean_variance"] <- mean(methVars)
      
      # Normalize methylation data from 0 to 1
      subDF_meth <- apply(subDF_meth[, methCols], 2, normalize)
      
      # Bind the gene methylation and gene expression data frames together
      subDF <- data.frame(cbind(subDF_meth[, methCols], subDF_cna, subDF_expr))
      
      # Remove any rows with missing expression data
      subDF <- subDF[!is.na(subDF$GeneExpression), ]
      
      # Skip if more than 20% of expression data is missing
      if(nrow(subDF) < 0.8 * model_sampleSize){next}
      
      # Flip the methylation data if necessary (so that higher methylation -> lower expression)
      for(k in methCols){
        # Get the median gene expression for tumors which have an above median methylation signature
        exprs_methHigh <- median(subDF$GeneExpression[subDF[, k] > median(subDF[, k])])
        # Get the median gene expression for tumors which have a below median methlyation signature
        exprs_methLow <- median(subDF$GeneExpression[subDF[, k] < median(subDF[, k])])
        # If the samples with a high methylation signature have a higher median gene expression level, flip the direction of the principal component
        if(exprs_methHigh > exprs_methLow){
          subDF[, k] <- 1 - subDF[, k]
        }
      }
      
      # If there is only one principal component available, the name ("PC1") is removed earlier; fix this
      if(ncol(subDF) == 3){colnames(subDF)[1] <- "PC1"}
      
      # Add log2(methylation * constant + 1) to the data frame
      for(k in methCols){
        subDF[, paste0("PC", k, "_log")] <- log2((subDF[, k] * logConst) + 1)
      }
      
      # Add non-log CN to the data frame (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/)
      subDF$CNA_noLog <- 2^subDF$CNA * 2
      
      ##################################################
      ## Models
      
      ## CNA linear regressions
      lm_logE_C <- lm(GeneExpression ~ CNA_noLog,
                      data = subDF)
      
      lm_logE_logC <- lm(GeneExpression ~ CNA,
                      data = subDF)
      
      ## CNA GAMs
      gam_logE_C <- gam(GeneExpression ~ s(CNA_noLog, k = gam_k),
                      data = subDF,
                      method = "REML")
      
      gam_logE_logC <- gam(GeneExpression ~ s(CNA, k = gam_k),
                         data = subDF,
                         method = "REML")
      
      
      ## Methylation linear regressions
      lm_logE_M <- lm(reformulate(termlabels = colnames(subDF)[methCols],
                                  response = "GeneExpression"),
                      data = subDF)
      
      lm_logE_logM <- lm(reformulate(termlabels = paste0(colnames(subDF)[methCols], "_log"),
                                  response = "GeneExpression"),
                      data = subDF)
      
      ## Methylation GAMs
      gam_logE_M <- gam(reformulate(termlabels = paste0("s(", colnames(subDF)[methCols], ", k = ", gam_k, ")"),
                                  response = "GeneExpression"),
                      data = subDF,
                      method = "REML")
      
      gam_logE_logM <- gam(reformulate(termlabels = paste0("s(", colnames(subDF)[methCols], "_log", ", k = ", gam_k, ")"),
                                     response = "GeneExpression"),
                         data = subDF,
                         method = "REML")
      
      ## Combined: only linear regressions
      lm_logE_C_M <- lm(reformulate(termlabels = c(colnames(subDF)[methCols], "CNA_noLog"),
                                 response = "GeneExpression"),
                     data = subDF)
      
      lm_logE_logC_M <- lm(reformulate(termlabels = c(colnames(subDF)[methCols], "CNA"),
                                    response = "GeneExpression"),
                        data = subDF)
      
      lm_logE_C_logM <- lm(reformulate(termlabels = c(paste0(colnames(subDF)[methCols], "_log"), "CNA_noLog"),
                                    response = "GeneExpression"),
                        data = subDF)
      
      lm_logE_logC_logM <- lm(reformulate(termlabels = c(paste0(colnames(subDF)[methCols], "_log"), "CNA"),
                                       response = "GeneExpression"),
                           data = subDF)
      
      
      ## Combined: only GAMs
      gam_logE_C_M <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], ", k = ", gam_k, ")"), paste0("s(CNA_noLog, k =", gam_k, ")")),
                                      response = "GeneExpression"),
                          data = subDF,
                          method = "REML")
      
      gam_logE_logC_M <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], ", k = ", gam_k, ")"), paste0("s(CNA, k =", gam_k, ")")),
                                         response = "GeneExpression"),
                             data = subDF,
                             method = "REML")
      
      
      gam_logE_C_logM <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], "_log, k = ", gam_k, ")"), paste0("s(CNA_noLog, k =" ,gam_k, ")")),
                                         response = "GeneExpression"),
                             data = subDF,
                             method = "REML")
      
      gam_logE_logC_logM <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], "_log, k = ", gam_k, ")"), paste0("s(CNA, k =" ,gam_k, ")")),
                                            response = "GeneExpression"),
                                data = subDF,
                                method = "REML")
      
      ## Combined: methylation GAM, linear CNA
      methGAM_cnaLM_logE_C_M <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], ", k = ", gam_k, ")"), "CNA_noLog"),
                                      response = "GeneExpression"),
                          data = subDF,
                          method = "REML")
      
      methGAM_cnaLM_logE_logC_M <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], ", k = ", gam_k, ")"), "CNA"),
                                         response = "GeneExpression"),
                             data = subDF,
                             method = "REML")
      
      methGAM_cnaLM_logE_C_logM <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], "_log, k = ", gam_k, ")"), "CNA_noLog"),
                                         response = "GeneExpression"),
                             data = subDF,
                             method = "REML")
      
      methGAM_cnaLM_logE_logC_logM <- gam(reformulate(termlabels = c(paste0("s(", colnames(subDF)[methCols], "_log, k = ", gam_k, ")"), "CNA"),
                                            response = "GeneExpression"),
                                data = subDF,
                                method = "REML")
      
      ## Combined: LM for methylation, spline for CNA
      cnaGAM_methLM_logE_C_M <- gam(reformulate(termlabels = c(colnames(subDF)[methCols], paste0("s(CNA_noLog, k =", gam_k, ")")),
                                      response = "GeneExpression"),
                          data = subDF,
                          method = "REML")
      
      cnaGAM_methLM_logE_logC_M <- gam(reformulate(termlabels = c(colnames(subDF)[methCols], paste0("s(CNA, k =", gam_k, ")")),
                                         response = "GeneExpression"),
                             data = subDF,
                             method = "REML")
      
      
      cnaGAM_methLM_logE_C_logM <- gam(reformulate(termlabels = c(paste0(colnames(subDF)[methCols], "_log"), paste0("s(CNA_noLog, k =" ,gam_k, ")")),
                                         response = "GeneExpression"),
                             data = subDF,
                             method = "REML")
      
      cnaGAM_methLM_logE_logC_logM <- gam(reformulate(termlabels = c(paste0(colnames(subDF)[methCols], "_log"), paste0("s(CNA, k =" ,gam_k, ")")),
                                            response = "GeneExpression"),
                                data = subDF,
                                method = "REML")
      
      ##################################
      # Enter data into the data frame
      
      # Expression variance
      allDF[currentGene, "expression_variance"] <- var(subDF$GeneExpression)
      
      # CNA variance
      allDF[currentGene, "cna_variance"] <- var(subDF$CNA)
      
      # Methylation: variance captured by PCs
      allDF[currentGene, "meth_variance_captured"] <- currentVariance[[currentGene]][max(methCols)]
      
      # Methylation: number of PCs used
      allDF[currentGene, "meth_pcs_used"] <- max(methCols)
      
      
      # 1
      # Rsq
      allDF[currentGene, "lm_logE_C_r2"] <- summary(lm_logE_C)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_C_r2_adj"] <- summary(lm_logE_C)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_C_aic"] <- AIC(lm_logE_C)
      
      
      # 2
      # Rsq
      allDF[currentGene, "lm_logE_logC_r2"] <- summary(lm_logE_logC)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_logC_r2_adj"] <- summary(lm_logE_logC)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_logC_aic"] <- AIC(lm_logE_logC)
      
      
      # 3
      # Rsq
      allDF[currentGene, "gam_logE_C_r2"] <- summary(gam_logE_C)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_C_r2_adj"] <- summary(gam_logE_C)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_C_aic"] <- AIC(gam_logE_C)
      
      
      # 4
      # Rsq
      allDF[currentGene, "gam_logE_logC_r2"] <- summary(gam_logE_logC)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_logC_r2_adj"] <- summary(gam_logE_logC)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_logC_aic"] <- AIC(gam_logE_logC)
      
      
      # 5
      # Rsq
      allDF[currentGene, "lm_logE_M_r2"] <- summary(lm_logE_M)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_M_r2_adj"] <- summary(lm_logE_M)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_M_aic"] <- AIC(lm_logE_M)
      
      
      # 6
      # Rsq
      allDF[currentGene, "lm_logE_logM_r2"] <- summary(lm_logE_logM)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_logM_r2_adj"] <- summary(lm_logE_logM)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_logM_aic"] <- AIC(lm_logE_logM)
      
      
      # 7
      # Rsq
      allDF[currentGene, "gam_logE_M_r2"] <- summary(gam_logE_M)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_M_r2_adj"] <- summary(gam_logE_M)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_M_aic"] <- AIC(gam_logE_M)
      
      
      # 8
      # Rsq
      allDF[currentGene, "gam_logE_logM_r2"] <- summary(gam_logE_logM)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_logM_r2_adj"] <- summary(gam_logE_logM)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_logM_aic"] <- AIC(gam_logE_logM)
      
      
      # 9
      # Rsq
      allDF[currentGene, "lm_logE_C_M_r2"] <- summary(lm_logE_C_M)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_C_M_r2_adj"] <- summary(lm_logE_C_M)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_C_M_aic"] <- AIC(lm_logE_C_M)
      
      
      # 10
      # Rsq
      allDF[currentGene, "lm_logE_logC_M_r2"] <- summary(lm_logE_logC_M)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_logC_M_r2_adj"] <- summary(lm_logE_logC_M)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_logC_M_aic"] <- AIC(lm_logE_logC_M)
      
      
      # 11
      # Rsq
      allDF[currentGene, "lm_logE_C_logM_r2"] <- summary(lm_logE_C_logM)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_C_logM_r2_adj"] <- summary(lm_logE_C_logM)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_C_logM_aic"] <- AIC(lm_logE_C_logM)
      
      
      # 12
      # Rsq
      allDF[currentGene, "lm_logE_logC_logM_r2"] <- summary(lm_logE_logC_logM)$r.squared
      
      # Rsq_adj
      allDF[currentGene, "lm_logE_logC_logM_r2_adj"] <- summary(lm_logE_logC_logM)$adj.r.squared
      
      # AIC
      allDF[currentGene, "lm_logE_logC_logM_aic"] <- AIC(lm_logE_logC_logM)
      
      
      # 13
      # Rsq
      allDF[currentGene, "gam_logE_C_M_r2"] <- summary(gam_logE_C_M)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_C_M_r2_adj"] <- summary(gam_logE_C_M)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_C_M_aic"] <- AIC(gam_logE_C_M)
      
      
      # 14
      # Rsq
      allDF[currentGene, "gam_logE_logC_M_r2"] <- summary(gam_logE_logC_M)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_logC_M_r2_adj"] <- summary(gam_logE_logC_M)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_logC_M_aic"] <- AIC(gam_logE_logC_M)
      
      
      # 15
      # Rsq
      allDF[currentGene, "gam_logE_C_logM_r2"] <- summary(gam_logE_C_logM)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_C_logM_r2_adj"] <- summary(gam_logE_C_logM)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_C_logM_aic"] <- AIC(gam_logE_C_logM)
      
      
      # 16
      # Rsq
      allDF[currentGene, "gam_logE_logC_logM_r2"] <- summary(gam_logE_logC_logM)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "gam_logE_logC_logM_r2_adj"] <- summary(gam_logE_logC_logM)$r.sq
      
      # AIC
      allDF[currentGene, "gam_logE_logC_logM_aic"] <- AIC(gam_logE_logC_logM)
      
      
      # 17
      # Rsq
      allDF[currentGene, "methGAM_cnaLM_logE_C_M_r2"] <- summary(methGAM_cnaLM_logE_C_M)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "methGAM_cnaLM_logE_C_M_r2_adj"] <- summary(methGAM_cnaLM_logE_C_M)$r.sq
      
      # AIC
      allDF[currentGene, "methGAM_cnaLM_logE_C_M_aic"] <- AIC(methGAM_cnaLM_logE_C_M)
      
      
      # 18
      # Rsq
      allDF[currentGene, "methGAM_cnaLM_logE_logC_M_r2"] <- summary(methGAM_cnaLM_logE_logC_M)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "methGAM_cnaLM_logE_logC_M_r2_adj"] <- summary(methGAM_cnaLM_logE_logC_M)$r.sq
      
      # AIC
      allDF[currentGene, "methGAM_cnaLM_logE_logC_M_aic"] <- AIC(methGAM_cnaLM_logE_logC_M)
      
      
      # 19
      # Rsq
      allDF[currentGene, "methGAM_cnaLM_logE_C_logM_r2"] <- summary(methGAM_cnaLM_logE_C_logM)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "methGAM_cnaLM_logE_C_logM_r2_adj"] <- summary(methGAM_cnaLM_logE_C_logM)$r.sq
      
      # AIC
      allDF[currentGene, "methGAM_cnaLM_logE_C_logM_aic"] <- AIC(methGAM_cnaLM_logE_C_logM)
      
      
      # 20
      # Rsq
      allDF[currentGene, "methGAM_cnaLM_logE_logC_logM_r2"] <- summary(methGAM_cnaLM_logE_logC_logM)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "methGAM_cnaLM_logE_logC_logM_r2_adj"] <- summary(methGAM_cnaLM_logE_logC_logM)$r.sq
      
      # AIC
      allDF[currentGene, "methGAM_cnaLM_logE_logC_logM_aic"] <- AIC(methGAM_cnaLM_logE_logC_logM)
      
      
      # 21
      # Rsq
      allDF[currentGene, "cnaGAM_methLM_logE_C_M_r2"] <- summary(cnaGAM_methLM_logE_C_M)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "cnaGAM_methLM_logE_C_M_r2_adj"] <- summary(cnaGAM_methLM_logE_C_M)$r.sq
      
      # AIC
      allDF[currentGene, "cnaGAM_methLM_logE_C_M_aic"] <- AIC(cnaGAM_methLM_logE_C_M)
      
      
      # 22
      # Rsq
      allDF[currentGene, "cnaGAM_methLM_logE_logC_M_r2"] <- summary(cnaGAM_methLM_logE_logC_M)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "cnaGAM_methLM_logE_logC_M_r2_adj"] <- summary(cnaGAM_methLM_logE_logC_M)$r.sq
      
      # AIC
      allDF[currentGene, "cnaGAM_methLM_logE_logC_M_aic"] <- AIC(cnaGAM_methLM_logE_logC_M)
      
      
      # 23
      # Rsq
      allDF[currentGene, "cnaGAM_methLM_logE_C_logM_r2"] <- summary(cnaGAM_methLM_logE_C_logM)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "cnaGAM_methLM_logE_C_logM_r2_adj"] <- summary(cnaGAM_methLM_logE_C_logM)$r.sq
      
      # AIC
      allDF[currentGene, "cnaGAM_methLM_logE_C_logM_aic"] <- AIC(cnaGAM_methLM_logE_C_logM)
      
      
      # 24
      # Rsq
      allDF[currentGene, "cnaGAM_methLM_logE_logC_logM_r2"] <- summary(cnaGAM_methLM_logE_logC_logM)$dev.expl
      
      # Rsq_adj
      allDF[currentGene, "cnaGAM_methLM_logE_logC_logM_r2_adj"] <- summary(cnaGAM_methLM_logE_logC_logM)$r.sq
      
      # AIC
      allDF[currentGene, "cnaGAM_methLM_logE_logC_logM_aic"] <- AIC(cnaGAM_methLM_logE_logC_logM)
      
    }
    
    # Remove genes which have been skipped
    allDF <- allDF[complete.cases(allDF$expression_variance), ]
    
    # Write output to file
      write.table(allDF,
                  file = paste0("./Output/TempFiles/AllModels_", w_BP_short, "kWindow_", tt, "_run", run, ".txt"),
                  sep = "\t",
                  quote = FALSE)
      
      # Print run statistics
      print(paste0(tt, ": run ", run, ", time taken: ", Sys.time() - a))
  }
}
