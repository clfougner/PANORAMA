##################################################
# Gather the outputs from the models into .rds objects (one for each tumor type)
##################################################

## Set working directory
setwd("/path/to/PANORAMA/")

##################################################
## Load/install libraries

if (!requireNamespace("doMC", quietly = TRUE))
  install.packages("doMC")
library(doMC)

if (!requireNamespace("gtools", quietly = TRUE))
  install.packages("gtools")
library(gtools)


##################################################
## Parameters

# Number of principal components to use
numPCs <- 5

# Number of base pairs around gene (window) for which a CpG should be classified as "belonging" to a gene
windowBP <- 50000
w_BP_short <- windowBP/1000

# Register cores
crs <- 1
registerDoMC(crs)
print(paste0("Cores: ", crs))

##################################################
## Functions

# Median with na.rm set to TRUE (for apply)
median_naRM <-function(inpt){
  return(median(inpt, na.rm = TRUE))
}

##################################################
## Gather data from all models

# Get list of all model files
modelFiles <- list.files(path = "./Output/TempFiles/", pattern = paste0("AllModels_", w_BP_short, "kWindow_"))

# Get the available tumor types
tts <- gsub(modelFiles, pattern =  paste0("AllModels_", w_BP_short, "kWindow_"), replacement = "")
tts <- gsub(tts, pattern =  "\\d", replacement = "")
tts <- unique(gsub(tts, pattern =  "_run.txt", replacement = ""))

# For each tumor type
foreach(i = 1:length(tts)) %dopar% {
  print(tts[i])
  
  # Get the files for the given tumor type
  currentFiles <- mixedsort(modelFiles[grep(x = modelFiles, pattern = tts[i])])
  
  # Read the first data frame
  allDat <- read.table(file = paste0("./Output/TempFiles/", currentFiles[1]),
                                      header = TRUE,
                                      stringsAsFactors = FALSE)
  
  # Get the number of columns in allDat
  colLengths <- ncol(allDat)
  
  # Set the run number in allDat
  allDat$run <- 1
  
  # Remove row names (else rbind will give issues with repeated row names)
  rownames(allDat) <- NULL
  
  # For each model file after the first one
  for(j in 2:length(currentFiles)){
    
    # Read the model file
    currentDat <- read.table(file = paste0("./Output/TempFiles/", currentFiles[j]),
                              header = TRUE,
                              stringsAsFactors = FALSE)
    
    # Set the run number
    currentDat$run <- j
    
    # Remove row names
    rownames(currentDat) <- NULL
    
    # Bind the outputs together
    allDat <- rbind(allDat, currentDat)
  }
  
  # Create an empty data frame for the median values
  medianDat <- data.frame(matrix(data = NA,
                                 ncol = ncol(allDat) - 1,
                                 nrow = length(unique(allDat$gene)),
                                 dimnames = list(unique(allDat$gene), colnames(allDat[1:ncol(allDat) - 1]))))
  
  # Enter data for the "gene" column into the data frame 
  medianDat$gene <- rownames(medianDat)
  
  ## Convert AIC to deltaAIC
  aicCols <- grep(colnames(medianDat), pattern = "aic", value = TRUE)
  
  # CNA
  aicCols_cna <- grep(aicCols, pattern = "C", value = TRUE)
  aicCols_cna <- aicCols_cna[-grep(aicCols_cna, pattern = "M")]
  
  cnaMins_logE <- apply(allDat[, aicCols_cna], 1, min)
  allDat[, aicCols_cna] <- allDat[, aicCols_cna] - cnaMins_logE
  
  # Methylation
  aicCols_meth <- grep(aicCols, pattern = "M", value = TRUE)
  aicCols_meth <- aicCols_meth[-grep(aicCols_meth, pattern = "C")]
  
  methMins_logE <- apply(allDat[, aicCols_meth], 1, min)
  allDat[, aicCols_meth] <- allDat[, aicCols_meth] - methMins_logE
  
  # Combined
  aicCols_combined <- grep(aicCols, pattern = "M", value = TRUE)
  aicCols_combined <- aicCols_combined[grep(aicCols_combined, pattern = "C")]
  
  combinedMins_logE <- apply(allDat[, aicCols_combined], 1, min)
  allDat[, aicCols_combined] <- allDat[, aicCols_combined] - combinedMins_logE
  
  # Get medians for all
  # For each gene
  for(j in medianDat$gene){
    
    # Get the median for each relevant column, for all runs for each gene
    currentGeneRowsMedian <- apply(allDat[allDat$gene == j, 2:colLengths], 2, FUN = median_naRM)
    
    # Enter the medians into the data frame
    medianDat[j, 2:colLengths] <- currentGeneRowsMedian
  }
  
  # Gather the medians and allDat into a single list
  together <- list(medians = medianDat, all = allDat)
  
  # Save to file as an R object
  saveRDS(together, file = paste0("./Output/AllModels_", w_BP_short, "kWindow_", tts[i], ".rds"))
}
