##################################################
# Visualize the results from the models
##################################################

## Set the working directory
setwd("/open/tmp/Christian/GeneExpressionModulation/")

##################################################
## Load libraries

if (!requireNamespace("vioplot", quietly = TRUE))
  install.packages("vioplot")
library(vioplot)

if (!requireNamespace("gridExtra", quietly = TRUE))
  install.packages("gridExtra")
library(gridExtra)

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  install.packages("ComplexHeatmap")
library(ComplexHeatmap)

if (!requireNamespace("circlize", quietly = TRUE))
  install.packages("circlize")
library(circlize)

##################################################
## Parameters
# Number of principal components to use
numPCs <- 5

# Sample size for models (and cutoff for number of tumors to use)
model_sampleSize <- 100

# Number of base pairs around gene (window) for which a CpG should be classified as "belonging" to a gene
windowBP <- 50000
w_BP_short <- windowBP/1000

# "Significant" difference for AIC
AIC_diffLimit <- 3

# Cut-off for percentile considered highly associated
highAssociationCutOff <- 0.8

# Colors for scatter plots
cols <- c("#FF8E1E", # orange
          "#2F9830", # green
          "#2178AE", # blue
          "#673C8E", # purple
          "#E11C1A") # red

# Point types for scatter plots
pointTypes <- rep(c(21:25), each = 5)

# CNA color
cnaColor <- "#046C9A"

# Methylation color
methylationColor <- "#4C856E"

# Combined model color
combinedColor <- "#C93312"


##################################################
## Functions

# Get p-value from a linear model
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

##################################################
## Sample data

# Read the sample overview
patientData <- read.table(file = "./Output/SampleOverview.txt", header = TRUE, stringsAsFactors = FALSE)

## Get number of patients/samples from each cancer type
# Create empty vector for sample numbers
allSampleNumbers <- c()

# For each tumor type in the data set
for(i in unique(patientData$TumorType)){
  
  # Get the number of tumors with data available for all required levels, for the given tumor type
  currentNum <- nrow(patientData[patientData$TumorType == i & patientData$inAll == 1, ])
  
  # Add to the vector
  allSampleNumbers <- c(allSampleNumbers, currentNum)
}

# Add names to the vector of sample numbers
names(allSampleNumbers) <- unique(patientData$TumorType)

# Order according to sample numbers
allSampleNumbers <- allSampleNumbers[order(names(allSampleNumbers))]

# Keep a vector of all sample numbers 
allSampleNumbersIncl <- allSampleNumbers

# Remove tumor types with too few samples
allSampleNumbers <- allSampleNumbers[allSampleNumbers > model_sampleSize]

# Names of tumor types
tts <- names(allSampleNumbers)

##################################################
## Models

# Get names of model files
modelFiles <- paste0("./Output/AllModels_", w_BP_short, "kWindow_", tts, ".rds")

# Read first file
modelsFirst <- readRDS(file = modelFiles[1])

# Get medians only (not individual runs)
models <- modelsFirst[[1]]

# Set row names to null
rownames(models) <- NULL

# Add tumor type to the row
models$TumorType <- tts[1]

# For each remaining file
for(i in 2:length(modelFiles)){

  # Read the file
  currentModels <- readRDS(file = modelFiles[i])

  # Get medians only (not individual runs)
  currentModels <- currentModels[[1]]

  # Set row names to null
  rownames(currentModels) <- NULL

  # Add tumor type to the row
  currentModels$TumorType <- tts[i]

  # Bind to all other runs
  models <- rbind(models, currentModels)

}

# Find Rsq_adj columns
adjCols <- grep(x = colnames(models), pattern = "adj", value = TRUE)

# Set negative Rsq_adj values to zero
for(i in adjCols){
  models[, i][models[, i] < 0] <- 0
}


##################################################
## Should GAMs or linear models be used?

## Get column names

# Get names of AIC columns
aicCols <- grep(colnames(models), pattern = "aic", value = TRUE)

# CNA
aicCols_cna <- grep(aicCols, pattern = "C", value = TRUE)
aicCols_cna <- aicCols_cna[-grep(aicCols_cna, pattern = "M")]

# Methylation
aicCols_meth <- grep(aicCols, pattern = "M", value = TRUE)
aicCols_meth <- aicCols_meth[-grep(aicCols_meth, pattern = "C")]

# Combined
aicCols_combined <- grep(aicCols, pattern = "M", value = TRUE)
aicCols_combined <- aicCols_combined[grep(aicCols_combined, pattern = "C")]

## Create empty data frames

# CNA, AIC
cna_agg_df_aic <- data.frame(matrix(data = NA,
                                nrow = length(aicCols_cna),
                                ncol = length(tts),
                                dimnames = list(aicCols_cna, tts)))

# CNA, Rsq_adj
cna_agg_df_rsq <- data.frame(matrix(data = NA,
                                nrow = length(aicCols_cna),
                                ncol = length(tts),
                                dimnames = list(gsub(aicCols_cna, pattern = "_aic", replacement = "_r2_adj"), tts)))

# Methylation, AIC
meth_agg_df_aic <- data.frame(matrix(data = NA,
                                    nrow = length(aicCols_meth),
                                    ncol = length(tts),
                                    dimnames = list(aicCols_meth, tts)))

# Methylation, Rsq_adj
meth_agg_df_rsq <- data.frame(matrix(data = NA,
                                    nrow = length(aicCols_meth),
                                    ncol = length(tts),
                                    dimnames = list(gsub(aicCols_meth, pattern = "_aic", replacement = "_r2_adj"), tts)))

# Combined, AIC
cmb_agg_df_aic <- data.frame(matrix(data = NA,
                                     nrow = length(aicCols_combined),
                                     ncol = length(tts),
                                     dimnames = list(aicCols_combined, tts)))

# Combined, Rsq_adj
cmb_agg_df_rsq <- data.frame(matrix(data = NA,
                                     nrow = length(aicCols_combined),
                                     ncol = length(tts),
                                     dimnames = list(gsub(aicCols_combined, pattern = "_aic", replacement = "_r2_adj"), tts)))

# For each tumor type
for(i in tts){

  # Subset data
  currentData <- models[models$TumorType == i, ]
    
  ## CNA
  # Make an empty data frame AIC and Rsq values
  cnaDF <- data.frame(matrix(data = NA,
                             nrow = length(aicCols_cna),
                             ncol = 4,
                             dimnames = list(aicCols_cna, c("Formula", "nBest", "r2_avg", "r2_adj_avg"))))
  
  # Enter the formula i.e. model type
  cnaDF$Formula <- gsub(rownames(cnaDF), pattern = "_aic", replacement = "")
  
  # For each model type
  for(j in rownames(cnaDF)){
    
    # Subset AIC values for the model
    currentVals <- currentData[, j]
    
    # Find the number of genes with delta AIC less than AIC_diffLimit and enter into the data frame
    cnaDF[j, "nBest"] <- length(currentVals[currentVals < AIC_diffLimit])
    
    # Get the correct column name for Rsq_adj columns
    currentRsqAdjName <- paste0(cnaDF[j, "Formula"], "_r2_adj")
    
    # Find the mean Rsq_adj for the model and enter into the data frame
    cnaDF[j, "r2_adj_avg"] <- mean(currentData[, currentRsqAdjName])
  }
  
  # Enter data into the aggregate (pan-cancer) data frames for AIC
  cna_agg_df_aic[, i] <- cnaDF[, "nBest"]
  
  # Enter data into the aggregate (pan-cancer) data frames for Rsq_adj
  cna_agg_df_rsq[, i] <- cnaDF[, "r2_adj_avg"]
  
  # Create vector for plotting number of genes with delta AIC less than AIC_diffLimit
  intermediate <- cnaDF$nBest
  
  # Set names for the above vector as the formula (model type)
  names(intermediate) <- cnaDF$Formula
  
  # Get vector of mean Rsq_adj (for pplotting), and order 
  rsqs_adj <- cnaDF$r2_adj_avg[rev(c(4, 3, 2, 1))]
  
  # Reorder the AIC vector
  intermediate <- intermediate[rev(c(4, 3, 2, 1))]
  
  # Start deltaAIC figure
  pdf(paste0("./Output/Figures/BestEquationCNA_",
             w_BP_short,
             "kWindow_", i, ".pdf"),
      pointsize = 8,
      width = 3,
      height = 2)
      
    # Set plotting parameters
    par(mar = c(bottom = 4,
                left = 7.5,
                top = 3,
                right = 2.5),
        xpd = TRUE)
      
    # Make barplot for number of genes with delta AIC less than AIC_diffLimit
    bp <- barplot(intermediate,
              horiz = TRUE,
              las = 2,
              space = c(0.6, 0.2, 0.6, 0.2),
              xlim = c(0, 15000),
              col = cnaColor,
              main = paste0("CNA: ", i))
    
    # Add Rsq_adj annotation
    text(y = bp, x = intermediate + 1500, labels = paste0(round(rsqs_adj, digits = 3)))
    
  # End figure
  dev.off()

  
  ## Methylation
  # Make an empty data frame AIC and Rsq values
  methDF <- data.frame(matrix(data = NA,
                             nrow = length(aicCols_meth),
                             ncol = 4,
                             dimnames = list(aicCols_meth, c("Formula", "nBest", "r2_avg", "r2_adj_avg"))))
  
  # Enter the formula i.e. model type
  methDF$Formula <- gsub(rownames(methDF), pattern = "_aic", replacement = "")
  
  # For each model type
  for(j in rownames(methDF)){
    
    # Subset AIC values for the model
    currentVals <- currentData[, j]
    
    # Find the number of genes with delta AIC less than AIC_diffLimit and enter into the data frame
    methDF[j, "nBest"] <- length(currentVals[currentVals < AIC_diffLimit])
    
    # Get the correct column name for Rsq_adj columns
    currentRsqAdjName <- paste0(methDF[j, "Formula"], "_r2_adj")
    
    # Find the mean Rsq_adj for the model and enter into the data frame
    methDF[j, "r2_adj_avg"] <- mean(currentData[, currentRsqAdjName])
  }
  
  # Enter data into the aggregate (pan-cancer) data frames for AIC
  meth_agg_df_aic[, i] <- methDF[, "nBest"]
  
  # Enter data into the aggregate (pan-cancer) data frames for Rsq_adj
  meth_agg_df_rsq[, i] <- methDF[, "r2_adj_avg"]
  
  # Create vector for plotting number of genes with delta AIC less than AIC_diffLimit
  intermediate <- methDF$nBest
  
  # Set names for the above vector as the formula (model type)
  names(intermediate) <- methDF$Formula
  
  # Get vector of mean Rsq_adj (for pplotting), and order 
  rsqs_adj <- methDF$r2_adj_avg[rev(c(3, 4, 1, 2))]
  
  # Reorder the AIC vector
  intermediate <- intermediate[rev(c(3, 4, 1, 2))]
  
  # Start deltaAIC figure
  pdf(paste0("./Output/Figures/BestEquationMeth_",
             w_BP_short,
             "kWindow_", i, ".pdf"),
      pointsize = 8,
      width = 3,
      height = 2)
  
    # Set plotting parameters
    par(mar = c(bottom = 4,
                left = 7,
                top = 3,
                right = 2.5),
        xpd = TRUE)
  
    # Make barplot for number of genes with delta AIC less than AIC_diffLimit
    bp <- barplot(intermediate,
                  horiz = TRUE,
                  las = 2,
                  space = c(0.6, 0.2, 0.6, 0.2),
                  col = methylationColor,
                  xlim = c(0, 15000),
                  main = paste0("Methylation: ", i))
    
    # Add Rsq_adj annotation
    text(y = bp, x = intermediate + 1500, labels = paste0(round(rsqs_adj, digits = 3)))
  
  # End figure
  dev.off()
  
  
  ## Combined model
  # Make an empty data frame AIC and Rsq values
  bothDF <- data.frame(matrix(data = NA,
                              nrow = length(aicCols_combined),
                              ncol = 4,
                              dimnames = list(aicCols_combined, c("Formula", "nBest", "r2_avg", "r2_adj_avg"))))
  
  # Enter the formula i.e. model type
  bothDF$Formula <- gsub(rownames(bothDF), pattern = "_aic", replacement = "")
  
  # For each model type
  for(j in rownames(bothDF)){
    
    # Subset AIC values for the model
    currentVals <- currentData[, j]
    
    # Find the number of genes with delta AIC less than AIC_diffLimit and enter into the data frame
    bothDF[j, "nBest"] <- length(currentVals[currentVals < AIC_diffLimit])
    
    # Get the correct column name for Rsq_adj columns
    currentRsqAdjName <- paste0(bothDF[j, "Formula"], "_r2_adj")
    
    # Find the mean Rsq_adj for the model and enter into the data frame
    bothDF[j, "r2_adj_avg"] <- mean(currentData[, currentRsqAdjName])
  }
  
  # Enter data into the aggregate (pan-cancer) data frames for AIC
  cmb_agg_df_aic[, i] <- bothDF[, "nBest"]
  
  # Enter data into the aggregate (pan-cancer) data frames for Rsq_adj
  cmb_agg_df_rsq[, i] <- bothDF[, "r2_adj_avg"]
  
  # Create vector for plotting number of genes with delta AIC less than AIC_diffLimit
  intermediate <- bothDF$nBest
  
  # Set names for the above vector as the formula (model type)
  names(intermediate) <- bothDF$Formula
  
  # Get vector of mean Rsq_adj (for pplotting), and order 
  rsqs_adj <- bothDF$r2_adj_avg[rev(c(6, 8, 5, 7,
                                  10, 12, 9, 11,
                                  14, 16, 13, 15,
                                  2, 4, 1, 3))]
  
  # Reorder the AIC vector
  intermediate <- intermediate[rev(c(6, 8, 5, 7,
                                 10, 12, 9, 11,
                                 14, 16, 13, 15,
                                 2, 4, 1, 3))]
  
  # Start deltaAIC figure
  pdf(paste0("./Output/Figures/BestEquationBoth_",
             w_BP_short,
             "kWindow_", i, ".pdf"),
      pointsize = 8,
      width = 4,
      height = 4)
  
    # Set plotting parameters
    par(mar = c(bottom = 4.5,
                left = 15,
                top = 3,
                right = 2.5),
        xpd = TRUE)
    
    # Make barplot for number of genes with delta AIC less than AIC_diffLimit
    bp <- barplot(intermediate,
                  horiz = TRUE,
                  las = 2,
                  space = c(0.6, 0.2, 0.2, 0.2),
                  xlim = c(0, 15000),
                  col = combinedColor,
                  main = paste0("Combined: ", i))
    
    # Add Rsq_adj annotation
    text(y = bp, x = intermediate + 1800, labels = paste0(round(rsqs_adj, digits = 3)))
  
  # End figure
  dev.off()
  
}


## CNA, pan-cancer deltaAIC and Rsq_adj
# Re-order pan-cancer AIC data frame
cna_agg_df_aic <- cna_agg_df_aic[aicCols_cna[c(4, 3, 2, 1)], ]

# Re-order pan-cancer Rsq_adj data frame
cna_agg_df_rsq <- cna_agg_df_rsq[rownames(cna_agg_df_rsq)[c(4, 3, 2, 1)], ]

# Start figure
pdf(paste0("./Output/Figures/CNA_AllAICs",
           w_BP_short,
           "kWindow.pdf"),
    width = 4.2,
    height = 1.8,
    pointsize = 8)

  # Heatmap
  Heatmap(as.matrix(cna_agg_df_aic),
          col = colorRamp2(c(0, max(cna_agg_df_aic)), c("white", cnaColor)),
          cluster_rows = FALSE,
          split = rep(c("A", "B"), each = 2),
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 7),
          cluster_columns = FALSE)

# End figure
dev.off()


## Methylation, pan-cancer deltaAIC and Rsq_adj
# Re-order pan-cancer AIC data frame
meth_agg_df_aic <- meth_agg_df_aic[aicCols_meth[c(3, 4, 1, 2)], ]

# Re-order pan-cancer Rsq_adj data frame
meth_agg_df_rsq <- meth_agg_df_rsq[rownames(meth_agg_df_rsq)[c(3, 4, 1, 2)], ]

# Start figure
pdf(paste0("./Output/Figures/Meth_AllAICs",
           w_BP_short,
           "kWindow.pdf"),
    width = 4.2,
    height = 1.8,
    pointsize = 8)
  
  # Heatmap
  Heatmap(as.matrix(meth_agg_df_aic),
          col = colorRamp2(c(0, max(meth_agg_df_aic)), c("white", methylationColor)),
          cluster_rows = FALSE,
          split = rep(c("A", "B"), each = 2),
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 7),
          cluster_columns = FALSE)
  
# End figure
dev.off()


## Combined, pan-cancer deltaAIC and Rsq_adj
# Re-order pan-cancer AIC data frame
cmb_agg_df_aic <- cmb_agg_df_aic[aicCols_combined[c(6, 8, 5, 7,
                                                    10, 12, 9, 11,
                                                    14, 16, 13, 15,
                                                    2, 4, 1, 3)], ]

# Re-order pan-cancer Rsq_adj data frame
cmb_agg_df_rsq <- cmb_agg_df_rsq[rownames(cmb_agg_df_rsq)[c(6, 8, 5, 7,
                                                            10, 12, 9, 11,
                                                            14, 16, 13, 15,
                                                            2, 4, 1, 3)], ]

# Start figure
pdf(paste0("./Output/Figures/Combined_AllAICs",
           w_BP_short,
           "kWindow.pdf"),
    width = 5,
    height = 4,
    pointsize = 8)

  # Heatmap
  Heatmap(as.matrix(cmb_agg_df_aic),
          col = colorRamp2(c(0, max(cmb_agg_df_aic)), c("white", combinedColor)),
          cluster_rows = FALSE,
          split = rep(c("A", "B", "C", "D"), each = 4),
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 7),
          cluster_columns = FALSE)
  
# End figure
dev.off()

## Write pan-cancer deltaAIC and Rsq_adj to file
# CNA, deltaAIC
write.table(x = cna_agg_df_aic, file = "./Output/cna_agg_df_aic.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# CNA, Rsq_dj
write.table(x = cna_agg_df_rsq, file = "./Output/cna_agg_df_rsq.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Meth, deltaAIC
write.table(x = meth_agg_df_aic, file = "./Output/meth_agg_df_aic.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Meth, Rsq_adj
write.table(x = meth_agg_df_rsq, file = "./Output/meth_agg_df_rsq.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Combined, deltaAIC
write.table(x = cmb_agg_df_aic, file = "./Output/cmb_agg_df_aic.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Combined, Rsq_adj
write.table(x = cmb_agg_df_rsq, file = "./Output/cmb_agg_df_rsq.txt", sep = "\t", quote = FALSE, row.names = TRUE)


##################################################
##################################################
## GAMs show a meaningful improvement over LMs for methylation data, but not CNA.
## Continue with GAM for methylation and LM for CNA, and mixed for the combined model
##################################################
##################################################

##################################################
## Model statistics

## Find mean methylation Rsq_adj for each tumor type
# Make empty vector
meth_rsq_mean <- c()

# For each tumor type
for(i in tts){
  
  # Subset model data for the tumor type, find mean Rsq_adj, add to vector
  meth_rsq_mean <- c(meth_rsq_mean, mean(models$gam_logE_M_r2_adj[models$TumorType == i]))
  
}

# Set vector names
names(meth_rsq_mean) <- tts

# Order values
meth_rsq_mean <- meth_rsq_mean[order(meth_rsq_mean)]


## Find mean CNA Rsq_adj for each tumor type
# Make empty vector
cna_rsq_mean <- c()

# For each tumor type
for(i in tts){
  
  # Subset model data for the tumor type, find mean Rsq_adj, add to vector
  cna_rsq_mean <- c(cna_rsq_mean, mean(models$lm_logE_logC_r2_adj[models$TumorType == i]))
  
}

# Set vector names
names(cna_rsq_mean) <- tts

# Order values
cna_rsq_mean <- cna_rsq_mean[order(cna_rsq_mean)]


## Find mean combined Rsq_adj for each tumor type
# Make empty vector
both_rsq_mean <- c()

# For each tumor type
for(i in tts){
  
  # Subset model data for the tumor type, find mean Rsq_adj, add to vector
  both_rsq_mean <- c(both_rsq_mean, mean(models$methGAM_cnaLM_logE_logC_M_r2_adj[models$TumorType == i]))
  
}

# Set vector names
names(both_rsq_mean) <- tts

# Order values
both_rsq_mean <- both_rsq_mean[order(both_rsq_mean)]


## Find mean expression variance
# Make empty vector
expressionMeanVariance <- c()

# For each tumor type
for(i in tts){
  
  # Subset model data for the tumor type, find mean gene expresion variance, add to vector
  expressionMeanVariance <- c(expressionMeanVariance, mean(models$expression_variance[models$TumorType == i]))
  
}

# Set vector names
names(expressionMeanVariance) <- tts

# Order values
expressionMeanVariance <- expressionMeanVariance[order(expressionMeanVariance)]


## Find mean CNA variance
# Make empty vector
cnaMeanVariance <- c()

# For each tumor type
for(i in tts){
  
  # Subset model data for the tumor type, find mean gene expresion variance, add to vector
  cnaMeanVariance <- c(cnaMeanVariance, mean(models$cna_variance[models$TumorType == i]))
  
}

# Set vector names
names(cnaMeanVariance) <- tts

# Order values
cnaMeanVariance <- cnaMeanVariance[order(cnaMeanVariance)]

## Find mean methylation variance
# Make empty vector
methMeanVariance <- c()

# For each tumor type
for(i in tts){
  
  # Subset model data for the tumor type, find mean gene expresion variance, add to vector
  methMeanVariance <- c(methMeanVariance, mean(models$meth_mean_variance[models$TumorType == i]))
  
  
}

# Set vector names
names(methMeanVariance) <- tts

# Order values
methMeanVariance <- methMeanVariance[order(methMeanVariance)]


## Genomic instability index
# Make empty vector
gii <- c()

# For each tumor type
for(i in tts){
  
  # Find the mean GII in the tumor type, add to the vector
  gii <- c(gii, mean(patientData$GII_ploidyCorrected[patientData$TumorType == i & patientData$inAll == 1], na.rm = TRUE))
  
}

# Set vector names
names(gii) <- tts

# Order values
gii <- gii[order(gii)]


## Create data frame with relevant data

## CNA Rsq_adj
# Order alphabetically by name of tumor type
avg_cna <- cna_rsq_mean[order(names(cna_rsq_mean))]

# Get names of tumor types
nms <- names(avg_cna)

# Set names of vector to null
names(avg_cna) <- NULL


## Methylation Rsq_adj
# Order alphabetically by name of tumor type
avg_meth <- meth_rsq_mean[order(names(meth_rsq_mean))]

# Set names of vector to null
names(avg_meth) <- NULL


## Combined Rsq_adj
# Order alphabetically by name of tumor type
avg_both <- both_rsq_mean[order(names(both_rsq_mean))]

# Set names of vector to null
names(avg_both) <- NULL


## CNA variance data
# Order alphabetically by name of tumor type
cnaVar <- cnaMeanVariance[order(names(cnaMeanVariance))]

# Set names of vector to null
names(cnaVar) <- NULL


## Methylation variance data
# Order alphabetically by name of tumor type
methVar <- methMeanVariance[order(names(methMeanVariance))]

# Set names of vector to null
names(methVar) <- NULL


## Gene expression variance data
# Order alphabetically by name of tumor type
exprVar <- expressionMeanVariance[order(names(expressionMeanVariance))]

# Set names of vector to null
names(exprVar) <- NULL

## Genomic instability index
# Order alphabetically by name of tumor type
gii_corr <- gii[order(names(gii))]

# Set names of vector to null
names(gii_corr) <- NULL


# Make aggregated data frame
longDF <- data.frame(TumorType = nms,
                     AvgMeth = avg_meth,
                     AvgCNA = avg_cna,
                     AvgBoth = avg_both,
                     VarMeth = methVar,
                     VarCNA  = cnaVar,
                     VarExpr = exprVar,
                     giiCor = gii_corr,
                     SampleNumbers = allSampleNumbers)


##################################################
##################################################
## Figures
##################################################
##################################################

##################################################
## CNA violin plot
# Start figure
pdf(paste0("./Output/Figures/Violin_CNA_Rsquare_adj_All_",
           w_BP_short,
           "kWindow.pdf"),
    pointsize = 7,
    height = 4.5,
    width = 2.25,
    useDingbats = FALSE)

  # Factor tumor type according to order of CNA Rsq_adj
  models$TumorType <- factor(models$TumorType, levels = names(cna_rsq_mean))
  
  # Make violin plot
  vioplot(lm_logE_logC_r2_adj ~ TumorType,
             horizontal = TRUE,
             data = models, 
             main = "Gene expression ~ CNA",
             font.main = 1,
             lty = 1,
             drawRect = FALSE,
             border = FALSE,
             ylim = c(0, 1),
             xlab = parse(text = '~italic(R[Adjusted]^"2")'),
             ylab = "",
             col = cnaColor,
             cex = 0.75,
             las = 1)
  
  # Add mean Rsq_adj dots
  stripchart(cna_rsq_mean ~ factor(names(cna_rsq_mean), levels = names(cna_rsq_mean)),
             add = TRUE,
             pch = 16,
             cex = 1.5,
             col = "white")

# End figure
dev.off()


##################################################
## Methylation violin plot
# Start figure
pdf(paste0("./Output/Figures/Violin_Meth_Rsquare_adj_All_",
           w_BP_short,
           "kWindow.pdf"),
    pointsize = 7,
    height = 4.5,
    width = 2.25,
    useDingbats = FALSE)
  
  # Factor tumor type according to order of methylation Rsq_adj
  models$TumorType <- factor(models$TumorType, levels = names(meth_rsq_mean))
  
  # Make violin plot
  vioplot(gam_logE_M_r2_adj ~ TumorType,
          horizontal = TRUE,
          data = models, 
          main = "Gene expression ~ Methylation",
          font.main = 1,
          lty = 1,
          drawRect = FALSE,
          border = FALSE,
          ylim = c(0, 1),
          xlab = parse(text = '~italic(R[Adjusted]^"2")'),
          ylab = "",
          col = methylationColor,
          cex = 0.75,
          las = 1)
  
  # Add mean Rsq_adj dots
  stripchart(meth_rsq_mean ~ factor(names(meth_rsq_mean), levels = names(meth_rsq_mean)),
             add = TRUE,
             pch = 16,
             cex = 1.5,
             col = "white")

# End figure
dev.off()


##################################################
## Combined violin plot
# Start figure
pdf(paste0("./Output/Figures/Violin_combined_Rsquare_adj_All_",
           w_BP_short,
           "kWindow.pdf"),
    pointsize = 7,
    height = 4.5,
    width = 2.25,
    useDingbats = FALSE)
  
  # Factor tumor type according to order of combined Rsq_adj
  models$TumorType <- factor(models$TumorType, levels = names(both_rsq_mean))
  
  # Make violin plot
  vioplot(methGAM_cnaLM_logE_logC_M_r2_adj ~ TumorType,
          horizontal = TRUE,
          data = models, 
          main = "Gene expression ~ Combined",
          font.main = 1,
          lty = 1,
          drawRect = FALSE,
          border = FALSE,
          ylim = c(0, 1),
          xlab = parse(text = '~italic(R[Adjusted]^"2")'),
          ylab = "",
          col = combinedColor,
          cex = 0.75,
          las = 1)
  
  # Add mean Rsq_adj dots
  stripchart(both_rsq_mean ~ factor(names(both_rsq_mean), levels = names(both_rsq_mean)),
             add = TRUE,
             pch = 16,
             cex = 1.5,
             col = "white")

# End figure
dev.off()


##################################################
## Methylation mean Rsq_adj vs. CNA mean Rsq_adj scatterplot

# Fit linear model
fit <- lm(longDF$AvgMeth ~ longDF$AvgCNA)

# Start figure
pdf(file = paste0("./Output/Figures/MethMeanRsqAdj_vs_CNAMeanRsqAdj_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)
  
  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))
  
  # Set color palette
  palette(cols)
  
  # Make plot
  plot(longDF$AvgMeth ~ longDF$AvgCNA,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.4),
       xlim = c(0, 0.15),
       cex = 1.5,
       las = 1,
       bty = "L",
       xlab = parse(text = paste('CNA: Mean ~italic(R[Adjusted]^"2")')),
       main = "Methylation ~ CNA",
       font.main = 1,
       ylab = parse(text = paste('Methylation: Mean ~italic(R[Adjusted]^"2")')))
  
  ## Rsq annotation
  # Get Rsq value
  rsq <- round(summary(fit)$r.squared, 2)
  
  # Add annotation
  text(x = 0.03, y = 0.4, pos = 4, labels = bquote(italic(R^2)==.(rsq)))
  
  
  ## P-value annotation
  # If P less than 0.001 
  if(lmp(fit) < 0.001){p <- " < 0.001"}
  
  # If P greater than 0.001
  if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
  
  # Add annotation
  text(x = 0.03, y = 0.37, pos = 4, labels = bquote(italic(P)~.(p)))
  
  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure 
dev.off()


##################################################
## GII vs CNA mean Rsq_adj scatterplot

# Fit linear model
fit <- lm(longDF$AvgCNA ~ longDF$giiCor)

# Predict data
sef <- predict(fit, se.fit = TRUE, interval = "confidence")

# Make data frame with predicted values and 95% confidence interval
d <- cbind(pred = sef$fit[, "fit"],
           se_up = sef$fit[, "lwr"],
           se_dn = sef$fit[, "upr"],
           vals = longDF$giiCor)

# Order data
d <- d[order(d[, "vals"]), ]

# Start figure
pdf(file = paste0("./Output/Figures/CNAMeanRsqAdj_vs_GII_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)
  
  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))
  
  # Set color palette
  palette(cols)

  # Make plot
  plot(longDF$AvgCNA ~ longDF$giiCor,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.15),
       xlim = c(0, 1),
       cex = 1.5,
       las = 1,
       bty = "L",
       xlab = "Mean GII",
       main = "CNA ~ GII",
       font.main = 1,
       ylab = parse(text = paste('CNA: Average ~italic(R[Adjusted]^"2")')))
  
  # Add 95% confidence interval
  polygon(x = c(d[, "vals"],
                d[order(d[, "vals"], decreasing = TRUE), "vals"]),
          y = c(d[, "se_up"],
                d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Add predicted line
  lines(pred ~ vals,
        data = d,
        lwd = 2,
        col = "black")
  
  ## Rsq annotation
  # Get Rsq value
  rsq <- round(summary(fit)$r.squared, 2)
  
  # Add annotation
  text(x = 0.15, y = 0.15, pos = 4, labels = bquote(italic(R^2)==.(rsq)))
  
  
  ## P-value annotation
  # If P less than 0.001 
  if(lmp(fit) < 0.001){p <- " < 0.001"}
  
  # If P greater than 0.001
  if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
  
  # Add annotation
  text(x = 0.15, y = 0.14, pos = 4, labels = bquote(italic(P)~.(p)))

  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


##################################################
## CNA mean Rsq_adj vs. CNA variance scatterplot

# Fit linear model
fit <- lm(longDF$AvgCNA ~ longDF$VarCNA)

# Predict data
sef <- predict(fit, se.fit = TRUE, interval = "confidence")

# Make data frame with predicted values and 95% confidence interval
d <- cbind(pred = sef$fit[, "fit"],
           se_up = sef$fit[, "lwr"],
           se_dn = sef$fit[, "upr"],
           vals = longDF$VarCNA)

# Order data
d <- d[order(d[, "vals"]), ]

# Start figure
pdf(file = paste0("./Output/Figures/CNAMeanRsqAdj_vs_CNA-variance_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))
  
  # Set color palette
  palette(cols)

  # Make plot
  plot(longDF$AvgCNA ~ longDF$VarCNA,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.16),
       xlim = c(0, 0.4),
       cex = 1.5,
       bty = "L",
       xlab = "Mean CNA variance",
       main = "CNA Rsq ~ CNA variance",
       font.main = 1,
       ylab = parse(text = paste('CNA: Average ~italic(R[Adjusted]^"2")')))

  # Add 95% confidence interval
  polygon(x = c(d[, "vals"],
                d[order(d[, "vals"], decreasing = TRUE), "vals"]),
          y = c(d[, "se_up"],
                d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Add predicted line
  lines(pred ~ vals,
        data = d,
        lwd = 2,
        col = "black")
  
  ## Rsq annotation
  # Get Rsq value
  rsq <- round(summary(fit)$r.squared, 2)
  
  # Add annotation
  text(x = 0.05, y = 0.16, pos = 4, labels = bquote(italic(R^2)==.(rsq)))
  
  
  ## P-value annotation
  # If P less than 0.001 
  if(lmp(fit) < 0.001){p <- " < 0.001"}
  
  # If P greater than 0.001
  if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
  
  # Add annotation
  text(x = 0.05, y = 0.15, pos = 4, labels = bquote(italic(P)~.(p)))
  
  
  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()

##################################################
## CNA mean Rsq_adj vs. expression variance scatterplot

# Fit linear model
fit <- lm(longDF$AvgCNA ~ longDF$VarExpr)

# Predict data
sef <- predict(fit, se.fit = TRUE, interval = "confidence")

# Make data frame with predicted values and 95% confidence interval
d <- cbind(pred = sef$fit[, "fit"],
           se_up = sef$fit[, "lwr"],
           se_dn = sef$fit[, "upr"],
           vals = longDF$VarExpr)

# Order data
d <- d[order(d[, "vals"]), ]

# Start figure
pdf(file = paste0("./Output/Figures/CNAMeanRsqAdj_vs_expression-variance_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))
  
  # Set color palette
  palette(cols)

  # Make plot
  plot(longDF$AvgCNA ~ longDF$VarExpr,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.15),
       xlim = c(0.8, 2),
       cex = 1.5,
       bty = "L",
       xlab = "Mean expression variance",
       main = "CNA Rsq ~ Expression variance",
       font.main = 1,
       ylab = parse(text = paste('CNA: Average ~italic(R[Adjusted]^"2")')))

  # Add 95% confidence interval
  polygon(x = c(d[, "vals"],
                d[order(d[, "vals"], decreasing = TRUE), "vals"]),
          y = c(d[, "se_up"],
                d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Add predicted line
  lines(pred ~ vals,
        data = d,
        lwd = 2,
        col = "black")
  
  ## Rsq annotation
  # Get Rsq value
  rsq <- round(summary(fit)$r.squared, 2)
  
  # Add annotation
  text(x = 1, y = 0.14, pos = 4, labels = bquote(italic(R^2)==.(rsq)))
  
  
  ## P-value annotation
  # If P less than 0.001 
  if(lmp(fit) < 0.001){p <- " < 0.001"}
  
  # If P greater than 0.001
  if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
  
  # Add annotation
  text(x = 1, y = 0.13, pos = 4, labels = bquote(italic(P)~.(p)))

  
  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


##################################################
## Meth mean Rsq_adj vs. meth variance scatterplot

# Fit linear model
fit <- lm(longDF$AvgMeth ~ longDF$VarMeth)

# Predict data
sef <- predict(fit, se.fit = TRUE, interval = "confidence")

# Make data frame with predicted values and 95% confidence interval
d <- cbind(pred = sef$fit[, "fit"],
           se_up = sef$fit[, "lwr"],
           se_dn = sef$fit[, "upr"],
           vals = longDF$VarMeth)

# Order data
d <- d[order(d[, "vals"]), ]

# Start figure
pdf(file = paste0("./Output/Figures/MethMeanRsqAdj_vs_meth-variance_",
                 w_BP_short,
                 "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))
  
  # Set color palette
  palette(cols)
  
  # Make plot
  plot(longDF$AvgMeth ~ longDF$VarMeth,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.4),
       xlim = c(0.04, 0.2),
       cex = 1.5,
       bty = "L",
       xlab = "Mean MethSig variance",
       main = "Methylation Rsq ~ MethSig variance",
       font.main = 1,
       ylab = parse(text = paste('Methylation: Average ~italic(R[Adjusted]^"2")')))
  
  # Add 95% confidence interval
  polygon(x = c(d[, "vals"],
                d[order(d[, "vals"], decreasing = TRUE), "vals"]),
          y = c(d[, "se_up"],
                d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Add predicted line
  lines(pred ~ vals,
        data = d,
        lwd = 2,
        col = "black")
  
  ## Rsq annotation
  # Get Rsq value
  rsq <- round(summary(fit)$r.squared, 2)
  
  # Add annotation
  text(x = 0.06, y = 0.4, pos = 4, labels = bquote(italic(R^2)==.(rsq)))
  
  
  ## P-value annotation
  # If P less than 0.001 
  if(lmp(fit) < 0.001){p <- " < 0.001"}
  
  # If P greater than 0.001
  if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
  
  # Add annotation
  text(x = 0.06, y = 0.37, pos = 4, labels = bquote(italic(P)~.(p)))

  
  # Add legend  
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()

##################################################
## Meth mean Rsq_adj vs. expression variance scatterplot

# Fit linear model
fit <- lm(longDF$AvgMeth ~ longDF$VarExpr)

# Start figure
pdf(file = paste0("./Output/Figures/MethMeanRsqAdj_vs_expression-variance_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)
  
  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))
  
  # Set color palette
  palette(cols)
  
  # Make plot
  plot(longDF$AvgMeth ~ longDF$VarExpr,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.4),
       xlim = c(0.8, 2),
       cex = 1.5,
       bty = "L",
       xlab = "Mean expression variance",
       main = "Methylation Rsq ~ Expression variance",
       font.main = 1,
       ylab = parse(text = paste('Methylation: Average ~italic(R[Adjusted]^"2")')))
  
  
  ## Rsq annotation
  # Get Rsq value
  rsq <- round(summary(fit)$r.squared, 2)
  
  # Add annotation
  text(x = 1, y = 0.4, pos = 4, labels = bquote(italic(R^2)==.(rsq)))
  
  
  ## P-value annotation
  # If P less than 0.001 
  if(lmp(fit) < 0.001){p <- " < 0.001"}
  
  # If P greater than 0.001
  if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
  
  # Add annotation
  text(x = 1, y = 0.37, pos = 4, labels = bquote(italic(P)~.(p)))
  
  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


##################################################
## Importance of sample numbers for Rsq_adj

## Empty vectors for models run with all samples
# CNA
avg_cna_allSamples <- c()

# Methylation
avg_meth_allSamples <- c()

# Combined model
avg_cmb_allSamples <- c()

# For each tumor type (including those with less than the minimum required number of samples)
for(i in names(allSampleNumbersIncl)){

  # Read model data for the tumor type, for models run with all samples
  currentData <- read.table(file = paste0("./Output/AllxSamplesxModels_", w_BP_short, "kWindow_", i, ".txt"))

  # Get mean Rsq_adj for CNA
  avg_cna_allSamples <- c(avg_cna_allSamples, mean(currentData$lm_logE_logC_r2_adj))

  # Get mean Rsq_adj for methylation
  avg_meth_allSamples <- c(avg_meth_allSamples, mean(currentData$gam_logE_M_r2_adj))

  # Get mean Rsq_adj for the combined model
  avg_cmb_allSamples <- c(avg_cmb_allSamples, mean(currentData$methGAM_cnaLM_logE_logC_M_r2_adj))
}

# Create a data frame with model data from above
allSamplesData <- data.frame(TumorType = names(allSampleNumbersIncl),
                             SampleNumbers = allSampleNumbersIncl,
                             avgMethAllSamples = avg_meth_allSamples,
                             avgCNAAllSamples = avg_cna_allSamples,
                             avgBothAllSamples = avg_cmb_allSamples)

# Colors for all tumor types
cols_allTTs <- c(rep(cols, times = ceiling(length(tts)/length(cols)))[1:length(tts)],
                 rep(c("#FFF302", "#BFBFBE"), times = length(allSampleNumbersIncl)))[1:length(allSampleNumbersIncl)]

# Point types for all tumor types
pointTypes_allTTs <- c(rep(c(21:25), each = ceiling(length(tts)/length(cols)))[1:length(tts)],
                       rep(c(21:25), times = length(allSampleNumbersIncl)))[1:length(allSampleNumbersIncl)]


# Reorder data (first tumor types with enough samples, then tumor types with too few samples)
allSamplesData <- rbind(allSamplesData[names(allSampleNumbersIncl) %in% tts, ],
                         allSamplesData[!names(allSampleNumbersIncl) %in% tts, ])

# Factor so that tumor types are plotted in correct order
allSamplesData$TumorType <- factor(allSamplesData$TumorType, levels = as.factor(allSamplesData$TumorType))

## Methylation: Corrected for sample size
# Start figure
pdf(file = paste0("./Output/Figures/MethMeanRsqAdj_vs_SampleNumber_",
                  model_sampleSize,
                  "Samples",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))

  # Set color palette
  palette(cols)

  # Make plot
  plot(longDF$AvgMeth ~ longDF$SampleNumbers,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.4),
       xlim = c(0, 800),
       cex = 1.5,
       bty = "L",
       xlab = "Sample number",
       main = "Methylation Rsq ~  Sample number (100)",
       font.main = 1,
       ylab = parse(text = paste('Methylation: Average ~italic(R[Adjusted]^"2")')))

  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


# Methylation: all samples included
pdf(file = paste0("./Output/Figures/MethMeanRsqAdj_vs_SampleNumber_AllSamples_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))

  # Set color palette
  palette(cols_allTTs)

  # Make plot
  plot(allSamplesData$avgMethAllSamples ~ allSamplesData$SampleNumbers,
       bg = allSamplesData$TumorType,
       col = "black",
       pch = pointTypes_allTTs,
       ylim = c(0, 0.4),
       xlim = c(0, 800),
       cex = 1.5,
       bty = "L",
       xlab = "Sample number",
       main = "Methylation Rsq ~  Sample number (All)",
       font.main = 1,
       ylab = parse(text = paste('Methylation: Average ~italic(R[Adjusted]^"2")')))


  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes_allTTs,
         legend = levels(allSamplesData$TumorType),
         pt.bg = cols_allTTs,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


## CNA: Corrected for sample size
# Start figure
pdf(file = paste0("./Output/Figures/CNAMeanRsqAdj_vs_SampleNumber_",
                  model_sampleSize,
                  "Samples",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))

  # Set color palette
  palette(cols)

  # Make plot
  plot(longDF$AvgCNA ~ longDF$SampleNumbers,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.15),
       xlim = c(0, 800),
       cex = 1.5,
       bty = "L",
       xlab = "Sample number",
       main = "CNA Rsq ~  Sample number (100)",
       font.main = 1,
       ylab = parse(text = paste('CNA: Average ~italic(R[Adjusted]^"2")')))


  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


## CNA: Not corrected for sample size
pdf(file = paste0("./Output/Figures/CNAMeanRsqAdj_vs_SampleNumber_AllSamples_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))

  # Set color palette
  palette(cols_allTTs)

  # Make plot
  plot(allSamplesData$avgCNAAllSamples ~ allSamplesData$SampleNumbers,
       bg = allSamplesData$TumorType,
       col = "black",
       pch = pointTypes_allTTs,
       ylim = c(0, 0.15),
       xlim = c(0, 800),
       cex = 1.5,
       bty = "L",
       xlab = "Sample number",
       main = "CNA Rsq ~  Sample number (All)",
       font.main = 1,
       ylab = parse(text = paste('CNA: Average ~italic(R[Adjusted]^"2")')))


  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes_allTTs,
         legend = levels(allSamplesData$TumorType),
         pt.bg = cols_allTTs,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()

## Combined: Corrected for sample size
# Start figure
pdf(file = paste0("./Output/Figures/CombinedMeanRsqAdj_vs_SampleNumber_",
                  model_sampleSize,
                  "Samples",
                  w_BP_short,
                  "kWindow.pdf"), width = 4, height = 3, pointsize = 8, useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))

  # Set color palette
  palette(cols)

  # Make plot
  plot(longDF$AvgBoth ~ longDF$SampleNumbers,
       bg = longDF$TumorType,
       col = "black",
       pch = pointTypes,
       ylim = c(0, 0.5),
       xlim = c(0, 800),
       cex = 1.5,
       bty = "L",
       xlab = "Sample number",
       main = "Combined Rsq ~  Sample number (100)",
       font.main = 1,
       ylab = parse(text = paste('Combined: Average ~italic(R[Adjusted]^"2")')))

  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes,
         legend = levels(longDF$TumorType),
         pt.bg = cols,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


## Combined: Not corrected for sample size
# Start figure
pdf(file = paste0("./Output/Figures/CombinedMeanRsqAdj_vs_SampleNumber_AllSamples_",
                  w_BP_short,
                  "kWindow.pdf"),
    width = 4,
    height = 3,
    pointsize = 8,
    useDingbats = FALSE)

  # Set plotting parameters
  par(mar = c(4, 5, 4, 9))

  # Set color palette
  palette(cols_allTTs)

  # Make plot
  plot(allSamplesData$avgBothAllSamples ~ allSamplesData$SampleNumbers,
       bg = allSamplesData$TumorType,
       col = "black",
       pch = pointTypes_allTTs,
       ylim = c(0, 0.5),
       xlim = c(0, 800),
       cex = 1.5,
       bty = "L",
       xlab = "Sample number",
       main = "Combined Rsq ~  Sample number (All)",
       font.main = 1,
       ylab = parse(text = paste('Combined: Average ~italic(R[Adjusted]^"2")')))


  # Add legend
  legend("right",
         inset = c(-0.5, 0),
         xpd = TRUE,
         pch = pointTypes_allTTs,
         legend = levels(allSamplesData$TumorType),
         pt.bg = cols_allTTs,
         col = "black",
         horiz = FALSE,
         ncol = 2)

# End figure
dev.off()


##################################################
## 2D histogram for methylation Rsq_adj vs CNA Rsq_adj at a per gene level

# For each tumor type
for(i in tts){
  
  # Subset data
  currentData <- models[models$TumorType == i, ]
  
  # Set the number of bins
  bins <- 50
  
  # Make a data frame of only zeros
  densData <- data.frame(matrix(data = 0,
                                nrow = bins,
                                ncol = bins))
  
  # For each gene
  for(j in 1:nrow(currentData)){
    
    # Find the nearest bin on the x-axis (CNA data)
    x <- (round(currentData$lm_logE_logC_r2_adj[j], digits = 2) * bins) + 1
    
    # Find the nearest bin on the y-axis (methylation data)
    y <- (round(currentData$gam_logE_M_r2_adj[j], digits = 2) * bins) + 1
    
    # Add one count to the appropriate bin
    densData[y, x] <- densData[y, x] + 1
  }
  
  # Flip axis
  densData <- densData[rev(rownames(densData)), ]
  
  # Log scale count data
  densData <- log2(densData + 1)
  
  # Start figure
  pdf(paste0("./Output/Figures/CNA_vs_Meth_densityPlot_",
             i,
             "_",
             w_BP_short,
             "kWindow.pdf"),
      width = 1.9,
      height = 1.2) 
  
    # Make heatmap
    p <- Heatmap(as.matrix(densData),
            name = paste0("hm", i),
            col = colorRamp2(c(0, 4, 8), c(adjustcolor("white", alpha = 0.00),
                                           adjustcolor("#3E9293", alpha = 1),
                                           adjustcolor("#337C7D", alpha = 1))),
            row_title = parse(text = paste('Methylation: ~italic(R[Adjusted]^"2")')),
            column_title = parse(text = paste('CNA: ~italic(R[Adjusted]^"2")')),
            column_title_side = "bottom",
            show_row_names = FALSE,
            show_column_names = FALSE,
            row_title_gp = gpar(fontsize = 7),
            column_title_gp = gpar(fontsize = 7),
            cluster_rows = FALSE,
            cluster_columns = FALSE)
    
    # Print heatmap
    print(p)
    
    # Add border
    decorate_heatmap_body(paste0("hm", i), {
      grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))
    })
  
  # End figure
  dev.off()
}

##################################################
## Is there a pattern of mutual exclusivity in whether genes are highly regulated by CNA and/or methylation?

# For cut-off levels between 0.5 and 0.9, in increments of 0.05
for(co in seq(from = 0.5, to = 0.9, by = 0.05)){
  
  # Create empty data frame for test of mutual exclusivity
  allOut <- data.frame(TumorType = tts,
                       p = NA,
                       OR = NA,
                       Trend = "NotSignificant",
                       stringsAsFactors = FALSE)
  
  # For each tumor type
  for(i in tts){
    
    # Subset the models
    currentData <- models[models$TumorType == i, ]
    
    # Find the Rsq_adj value at the specified percentile cut-off for methylation
    cutOffMeth <- quantile(currentData$gam_logE_M_r2_adj, co)
    
    # Find the Rsq_adj value at the specified percentile cut-off for CNA
    cutOffCNA <- quantile(currentData$lm_logE_logC_r2_adj, co)
    
    # Number of meth-high, CNA-low genes
    methOnlyPos <- nrow(currentData[(currentData$gam_logE_M_r2_adj > cutOffMeth & currentData$lm_logE_logC_r2_adj < cutOffCNA), ])
    
    # Number of meth-low, CNA-high genes
    cnaOnlyPos <- nrow(currentData[(currentData$gam_logE_M_r2_adj < cutOffMeth & currentData$lm_logE_logC_r2_adj > cutOffCNA), ])
    
    # Number of meth-low, CNA-low genes
    bothNeg <- nrow(currentData[(currentData$gam_logE_M_r2_adj < cutOffMeth & currentData$lm_logE_logC_r2_adj < cutOffCNA), ])
    
    # Number of meth-high, CNA-high genes
    bothPos <-  nrow(currentData[(currentData$gam_logE_M_r2_adj > cutOffMeth & currentData$lm_logE_logC_r2_adj > cutOffCNA), ])
  
    # Perform fisher's exact test on 2x2 contingency table
    out <- fisher.test(x = matrix(data = c(methOnlyPos, bothPos,
                                          bothNeg, cnaOnlyPos),
                                  byrow = TRUE,
                                  ncol = 2))
    
    # Enter p-value into data frame
    allOut$p[allOut$TumorType == i] <- out$p.value
    
    # Enter odds ratio into data frame
    allOut$OR[allOut$TumorType == i] <- out[["estimate"]][["odds ratio"]]
    
    # If odds-ratio greater than 1 and p-value less than 0.05
    if(out[["estimate"]][["odds ratio"]] > 1 & out$p.value < 0.05){
      
      # Trend for mutual exclusivity, enter into data frame
      allOut$Trend[allOut$TumorType == i] <- "MutuallyExclusive"
    
    # If odds-ratio less than 1 and p-value less than 0.05
    } else if(out[["estimate"]][["odds ratio"]] < 1 & out$p.value < 0.05){
      
      # Trend for mutual inclusivity, enter into data frame
      allOut$Trend[allOut$TumorType == i] <- "MutuallyInclusive"
      
    }
  }
  
  # Order by odds ratio
  allOut <- allOut[order(allOut$OR),]
  
  # Factor so that tumor types are plotted in correct order
  allOut$TumorType <- factor(allOut$TumorType, levels = allOut$TumorType)
  
  # Start figure
  pdf(paste0("./Output/Figures/GAM_MutualExclusivityOverview_",
             w_BP_short,
             "kWindow_",
             co,
             "cutoff.pdf"),
      width = 2.25,
      height = 3.5,
      pointsize = 7,
      useDingbats = FALSE)
  
    # Color palette for the figure
    colstwo <- c("#1E7AB2", # Blue, mutually exclusive
                 "#BC3839", # Red, mutually inclusive
                  "#000000") # Black, not significant
    
    # Make chart
    stripchart(OR ~ TumorType,
           data = allOut,
           vertical = FALSE,
           las = 1,
           pch = 16,
           main = "",
           xlab = "Odds ratio")
    
    # Dotted line to show odds ratio of 1
    segments(x0 = 1, x1 = 1, y0 = -0.5, y1 = 33, lty = 2, col = "darkgrey")
    
    # Add dots so that they are placed on top of the dotted line
    stripchart(OR ~ TumorType,
               data = allOut,
               add = TRUE,
               col = colstwo[factor(allOut$Trend)],
               vertical = FALSE,
               las = 1,
               pch = 16,
               xlab = "Odds ratio")
    
    # Add legend
    legend("bottomright",
           pch = 16,
           legend = levels(factor(c("Mutually exclusive", "Mutually inclusive", "Not significant"))),
           col = colstwo,
           horiz = FALSE,
           bty = "n")
    
  # Ennd figure
  dev.off()
}
 

##################################################
## To what extent is CNA_Rsq_adj determined by CNA variance?
# Start figure
pdf(paste0("./Output/Figures/CNA-RsqAdj_vs_CNA-Variance_",
           w_BP_short,
           "kWindow.pdf"),
    pointsize = 9,
    width = 8.27,
    height = 11.69,
    useDingbats = FALSE)
  
  # Set plotting parameters
  par(mfrow = c(ceiling(length(tts)/3), 3), mar = c(4, 5, 2, 0.4))
  
  # Make a data frame to save mean CN variance in MethHigh and MethLow groups
  cnaVarMeans <- data.frame(TumorType = c(rep(tts, each = 3)),
                            class = rep(c("All", "MethHigh", "MethLow")),
                            values = NA)
  
  # For each tumor type + pan-cancer
  for(i in c(tts, "Pan-cancer")){
    
    # If not pan-cancer
    if(i != "Pan-cancer"){
      
      # Subset data for current tumor type
      currentData <- models[models$TumorType == i, ]
      
      # Subset models with high E-M association
      currentDataMethHigh <- currentData[currentData$gam_logE_M_r2_adj > quantile(currentData$gam_logE_M_r2_adj, highAssociationCutOff), ]
      
      # Subset models with low E-M association
      currentDataMethLow <- currentData[!currentData$gam_logE_M_r2_adj > quantile(currentData$gam_logE_M_r2_adj, highAssociationCutOff), ]
    }
    
    # If pan-cancer
    if(i == "Pan-cancer"){
      
      # Set currentData as all models
      currentData <- models
      
      # Make an empty data frame for meth-high genes
      currentDataMethHigh <- currentData[NULL, ]
      
      # Make an empty data frame for meth-low genes
      currentDataMethLow <- currentData[NULL, ]
      
      # Get meth-high/meth-low genes for each tumor type
      for(j in tts){
        
        # Subset data for the tumor type
        currentDataTT <- currentData[currentData$TumorType == j, ]
        
        # Subset meth-high genes for the tumor type
        currentDataMethHighTT <- currentDataTT[currentDataTT$gam_logE_M_r2_adj > quantile(currentDataTT$gam_logE_M_r2_adj, highAssociationCutOff), ]
        
        # Bind meth-high genes to the pan-cancer data frame
        currentDataMethHigh <- rbind(currentDataMethHigh, currentDataMethHighTT)
        
        # Subset meth-low genes for the tumor type
        currentDataMethLowTT <- currentDataTT[!currentDataTT$gam_logE_M_r2_adj > quantile(currentDataTT$gam_logE_M_r2_adj, highAssociationCutOff), ]
        
        # Bind meth-low genes to the pan-cancer data frame
        currentDataMethLow <- rbind(currentDataMethLow, currentDataMethLowTT)
        
      }
    }
    
    ## All genes
    # Fit linear model model (Copy number Rsq_adj ~ copy number variance)
    fit <- lm(currentData$lm_logE_logC_r2_adj ~ currentData$cna_variance)
    
    # Predict data
    sef <- predict(fit, se.fit = TRUE, interval = "confidence")
    
    # Make data frame for predicted data and 95% confidence interval
    d <- cbind(pred = sef$fit[, "fit"],
               se_up = sef$fit[, "lwr"],
               se_dn = sef$fit[, "upr"],
               vals = currentData$cna_variance)
    
    # Order data
    d <- d[order(d[, "vals"]), ]
    
    ## Meth-high genes
    # Fit linear model model (Copy number Rsq_adj ~ copy number variance)
    fitMethHigh <- lm(currentDataMethHigh$lm_logE_logC_r2_adj ~ currentDataMethHigh$cna_variance)
    
    # Predict data
    sefMethHigh <- predict(fitMethHigh, se.fit = TRUE, interval = "confidence")
    
    # Make data frame for predicted data and 95% confidence interval
    dMethHigh <- cbind(pred = sefMethHigh$fit[, "fit"],
               se_up = sefMethHigh$fit[, "lwr"],
               se_dn = sefMethHigh$fit[, "upr"],
               vals = currentDataMethHigh$cna_variance)
    
    # Order data
    dMethHigh <- dMethHigh[order(dMethHigh[, "vals"]), ]
    
    
    ## Meth-low genes
    # Fit linear model model (Copy number Rsq_adj ~ copy number variance)
    fitMethLow <- lm(currentDataMethLow$lm_logE_logC_r2_adj ~ currentDataMethLow$cna_variance)
    
    # Predict data
    sefMethLow <- predict(fitMethLow, se.fit = TRUE, interval = "confidence")
    
    # Make data frame for predicted data and 95% confidence interval
    dMethLow <- cbind(pred = sefMethLow$fit[, "fit"],
                       se_up = sefMethLow$fit[, "lwr"],
                       se_dn = sefMethLow$fit[, "upr"],
                       vals = currentDataMethLow$cna_variance)
    
    # Order data
    dMethLow <- dMethLow[order(dMethLow[, "vals"]), ]
    
    # If not-pancancer
    if(i != "Pan-cancer"){
      # Enter mean CN variance for all genes into cnaVarMeans data frame
      cnaVarMeans$values[cnaVarMeans$TumorType == i & cnaVarMeans$class == "All"] <- mean(currentData$cna_variance)
      
      # Enter mean CN variance for meth-high genes into cnaVarMeans data frame
      cnaVarMeans$values[cnaVarMeans$TumorType == i & cnaVarMeans$class == "MethHigh"] <- mean(currentDataMethHigh$cna_variance)
      
      # Enter mean CN variance for meth-low genes into cnaVarMeans data frame
      cnaVarMeans$values[cnaVarMeans$TumorType == i & cnaVarMeans$class == "MethLow"] <- mean(currentDataMethLow$cna_variance)
    }
    
    # Make plot
    plot(currentData$lm_logE_logC_r2_adj ~ currentData$cna_variance,
         pch = 16,
         font.main = 1,
         col = adjustcolor(cnaColor, alpha.f = 0.25),
         xlab = "",
         ylab = "",
         las = 1,
         main = i,
         cex = 0.5,
         lwd = 0.75,
         bty = "n",
         ylim = c(0,1),
         xlim = c(0, 2))
    
    ## Add annotation dots to top of plot
    # Data frame for dot locations
    pointsDF <- data.frame(x = c(0.1, 0.8, 1.5), y = 0.875)
    
    # Plot points
    points(y ~ x,
           data = pointsDF,
           pch = 16,
           cex = 2,
           col = c("black", methylationColor, cols[4]))
    
    ## Linear regression for all genes
    # Plot 95% confidence interval
    polygon(x = c(d[, "vals"],
                  d[order(d[, "vals"], decreasing = TRUE), "vals"]),
            y = c(d[, "se_up"],
                  d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
            col = adjustcolor("black", alpha.f = 0.15),
            border = NA)
    
    # Plot predicted line
    lines(pred ~ vals,
          data = d,
          lwd = 2,
          col = "black")
    
    # Find Rsq for model
    rsq <-round(summary(fit)$r.squared, 2)
    
    # Annotate Rsq
    text(x = 0.125,
         y = 0.95,
         pos = 4,
         labels = bquote(italic(R^2)==.(rsq)))
    
    # Find p-value for model
    if(lmp(fit) < 0.001){p <- " < 0.001"}
    if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
    
    # Add p-value
    text(x = 0.125,
         y = 0.8,
         pos = 4,
         labels = bquote(italic(P)~.(p)))
    

    ## Meth-high genes linear regression
    # Plot 95% confidence interval
    polygon(x = c(dMethHigh[, "vals"],
                  dMethHigh[order(dMethHigh[, "vals"], decreasing = TRUE), "vals"]),
            y = c(dMethHigh[, "se_up"],
                  dMethHigh[order(dMethHigh[, "vals"], decreasing = TRUE), "se_dn"]),
            col = adjustcolor("black", alpha.f = 0.15),
            border = NA)
    
    # Plot predicted line
    lines(pred ~ vals,
          data = dMethHigh,
          lwd = 2,
          col = methylationColor)
    
    # Find Rsq for model
    rsq <-round(summary(fitMethHigh)$r.squared, 2)
    
    # Annotate Rsq
    text(x = 0.825,
         y = 0.95,
         pos = 4,
         labels = bquote(italic(R^2)==.(rsq)))
    
    # Find p-value for model 
    if(lmp(fitMethHigh) < 0.001){p <- " < 0.001"}
    if(!lmp(fitMethHigh) < 0.001){p <- paste0(" = ", round(lmp(fitMethHigh), 3))}
    
    # Add p-value
    text(x = 0.825,
         y = 0.8,
         pos = 4,
         labels = bquote(italic(P)~.(p)))
    

    ## Meth-low genes linear regression  
    # Plot 95% confidence interval
    polygon(x = c(dMethLow[, "vals"],
                  dMethLow[order(dMethLow[, "vals"], decreasing = TRUE), "vals"]),
            y = c(dMethLow[, "se_up"],
                  dMethLow[order(dMethLow[, "vals"], decreasing = TRUE), "se_dn"]),
            col = adjustcolor("black", alpha.f = 0.15),
            border = NA)
    
    # Plot predicted line
    lines(pred ~ vals,
          data = dMethLow,
          lwd = 2,
          col = cols[4])
    
    # Find Rsq for model
    rsq <-round(summary(fitMethLow)$r.squared, 2)
    
    # Annotate Rsq
    text(x = 1.525,
         y = 0.95,
         pos = 4,
         labels = bquote(italic(R^2)==.(rsq)))
    
    # Find p-value for model
    if(lmp(fitMethLow) < 0.001){p <- " < 0.001"}
    if(!lmp(fitMethLow) < 0.001){p <- paste0(" = ", round(lmp(fitMethLow), 3))}
    
    # Add p-value
    text(x = 1.525,
         y = 0.8,
         pos = 4,
         labels = bquote(italic(P)~.(p)))
    
  }

# End figure
dev.off()

##################################################
## To what extent is CNA_Rsq_adj determined by CNA variance? - minified

# For each tumor type + pan-cancer
for(i in c(tts, "Pan-cancer")){
  
  # Start figure
  pdf(paste0("./Output/Figures/MinifiedCNA-Variance_vs_CNA-Rsquare_",
             i,
             "_",
             w_BP_short,
             "kWindow.pdf"),
      pointsize = 8,
      width = 1.7,
      height = 2,
      useDingbats = FALSE)
  
  # If not pan-cancer
  if(i != "Pan-cancer"){
    
    # Subset data for current tumor type
    currentData <- models[models$TumorType == i, ]
    
    # Subset models with high E-M association
    currentDataMethHigh <- currentData[currentData$gam_logE_M_r2_adj > quantile(currentData$gam_logE_M_r2_adj, highAssociationCutOff), ]
    
    # Subset models with low E-M association
    currentDataMethLow <- currentData[!currentData$gam_logE_M_r2_adj > quantile(currentData$gam_logE_M_r2_adj, highAssociationCutOff), ]
    
  }
  
  # If pan-cancer
  if(i == "Pan-cancer"){
    
    # Set currentData as all models
    currentData <- models

    # Make an empty data frame for meth-high genes
    currentDataMethHigh <- currentData[NULL, ]
    
    # Make an empty data frame for meth-low genes
    currentDataMethLow <- currentData[NULL, ]
    
    # Get meth-high/meth-low genes for each tumor type
    for(j in tts){
      
      # Subset data for the tumor type
      currentDataTT <- currentData[currentData$TumorType == j, ]
      
      # Subset meth-high genes for the tumor type
      currentDataMethHighTT <- currentDataTT[currentDataTT$gam_logE_M_r2_adj > quantile(currentDataTT$gam_logE_M_r2_adj, highAssociationCutOff), ]
      
      # Bind meth-high genes to the pan-cancer data frame
      currentDataMethHigh <- rbind(currentDataMethHigh, currentDataMethHighTT)
      
      # Subset meth-low genes for the tumor type
      currentDataMethLowTT <- currentDataTT[!currentDataTT$gam_logE_M_r2_adj > quantile(currentDataTT$gam_logE_M_r2_adj, highAssociationCutOff), ]
      
      # Bind meth-low genes to the pan-cancer data frame
      currentDataMethLow <- rbind(currentDataMethLow, currentDataMethLowTT)
      
    }
  }
  
  
  ## All genes
  # Fit linear model model (Copy number Rsq_adj ~ copy number variance)
  fit <- lm(currentData$lm_logE_logC_r2_adj ~ currentData$cna_variance)
  
  # Predict data
  sef <- predict(fit, se.fit = TRUE, interval = "confidence")
  
  # Make data frame for predicted data and 95% confidence interval
  d <- cbind(pred = sef$fit[, "fit"],
             se_up = sef$fit[, "lwr"],
             se_dn = sef$fit[, "upr"],
             vals = currentData$cna_variance)
  
  # Order data
  d <- d[order(d[, "vals"]), ]
  
  
  ## Meth-high genes
  # Fit linear model model (Copy number Rsq_adj ~ copy number variance)
  fitMethHigh <- lm(currentDataMethHigh$lm_logE_logC_r2_adj ~ currentDataMethHigh$cna_variance)
  
  # Predict data
  sefMethHigh <- predict(fitMethHigh, se.fit = TRUE, interval = "confidence")
  
  # Make data frame for predicted data and 95% confidence interval
  dMethHigh <- cbind(pred = sefMethHigh$fit[, "fit"],
                     se_up = sefMethHigh$fit[, "lwr"],
                     se_dn = sefMethHigh$fit[, "upr"],
                     vals = currentDataMethHigh$cna_variance)
  
  # Order data
  dMethHigh <- dMethHigh[order(dMethHigh[, "vals"]), ]
  
  
  ## Meth-low genes
  # Fit linear model model (Copy number Rsq_adj ~ copy number variance)
  fitMethLow <- lm(currentDataMethLow$lm_logE_logC_r2_adj ~ currentDataMethLow$cna_variance)
  
  # Predict data
  sefMethLow <- predict(fitMethLow, se.fit = TRUE, interval = "confidence")
  
  # Make data frame for predicted data and 95% confidence interval
  dMethLow <- cbind(pred = sefMethLow$fit[, "fit"],
                    se_up = sefMethLow$fit[, "lwr"],
                    se_dn = sefMethLow$fit[, "upr"],
                    vals = currentDataMethLow$cna_variance)
  
  # Order data
  dMethLow <- dMethLow[order(dMethLow[, "vals"]), ]
  
  # Make plots
  plot(currentData$lm_logE_logC_r2_adj[1] ~ currentData$cna_variance[1],
       pch = 16,
       font.main = 1,
       col = "white",
       xlab = "",
       ylab = "",
       las = 1,
       main = "",
       cex = 0.5,
       lwd = 0.75,
       ylim = c(0,0.5),
       xlim = c(0, 1))
  
  ## All genes
  # Plot 95% confidence interval
  polygon(x = c(d[, "vals"],
                d[order(d[, "vals"], decreasing = TRUE), "vals"]),
          y = c(d[, "se_up"],
                d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Plot predicted line
  lines(pred ~ vals,
        data = d,
        lwd = 2,
        col = "black")
  
  
  ## Meth-high genes
  # Plot 95% confidence interval
  polygon(x = c(dMethHigh[, "vals"],
                dMethHigh[order(dMethHigh[, "vals"], decreasing = TRUE), "vals"]),
          y = c(dMethHigh[, "se_up"],
                dMethHigh[order(dMethHigh[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Plot predicted line
  lines(pred ~ vals,
        data = dMethHigh,
        lwd = 2,
        col = methylationColor)
  
  ## Meth-low genes
  # Plot 95% confidence interval
  polygon(x = c(dMethLow[, "vals"],
                dMethLow[order(dMethLow[, "vals"], decreasing = TRUE), "vals"]),
          y = c(dMethLow[, "se_up"],
                dMethLow[order(dMethLow[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Plot predicted line
  lines(pred ~ vals,
        data = dMethLow,
        lwd = 2,
        col = cols[4])
  
  # End figure
  dev.off()
  
}


##################################################
## To what extent is Methylation_Rsq_adj determined by methylation variance?

# Start figure
pdf(paste0("./Output/Figures/Meth-Variance_vs_Meth-Rsquare_",
           w_BP_short,
           "kWindow.pdf"),
    pointsize = 9,
    width = 8.27,
    height = 11.69,
    useDingbats = FALSE)
  
  # Set plotting parameters
  par(mfrow = c(ceiling(length(tts)/3), 3), mar = c(4, 5, 2, 0.4))
  
  # Make a data frame to save mean MethSig variance in CNAHigh and CNALow groups
  methVarMeans <- data.frame(TumorType = c(rep(tts, each = 3)),
                             class = rep(c("All", "cnaHigh", "cnaLow")),
                             values = NA)
  
  # For each tumor type + pan-cancer
  for(i in c(tts, "Pan-cancer")){
    
    # If not pan-cancer
    if(i != "Pan-cancer"){
      
      # Subset data for current tumor type
      currentData <- models[models$TumorType == i, ]
      
      # Subset models with high E-C association
      currentDataCNAHigh <- currentData[currentData$lm_logE_logC_r2_adj > quantile(currentData$lm_logE_logC_r2_adj, highAssociationCutOff), ]
      
      # Subset models with low E-C association
      currentDataCNALow <- currentData[!currentData$lm_logE_logC_r2_adj > quantile(currentData$lm_logE_logC_r2_adj, highAssociationCutOff), ]
      
    }
    
    # If pan-cancer 
    if(i == "Pan-cancer"){
      
      # Set currentData as all models
      currentData <- models
      
      # Make an empty data frame for CNA-high genes
      currentDataCNAHigh <- currentData[NULL, ]
      
      # Make an empty data frame for CNA-low genes
      currentDataCNALow <- currentData[NULL, ]
      
      # Get CNA-high/CNA-low genes for each tumor type
      for(j in tts){
        
        # Subset data for the tumor type
        currentDataTT <- currentData[currentData$TumorType == j, ]
        
        # Subset CNA-high genes for the tumor type
        currentDataCNAHighTT <- currentDataTT[currentDataTT$lm_logE_logC_r2_adj > quantile(currentDataTT$lm_logE_logC_r2_adj, highAssociationCutOff), ]
        
        # Bind CNA-high genes to the pan-cancer data frame
        currentDataCNAHigh <- rbind(currentDataCNAHigh, currentDataCNAHighTT)
        
        # Subset CNA-low genes for the tumor type
        currentDataCNALowTT <- currentDataTT[!currentDataTT$lm_logE_logC_r2_adj > quantile(currentDataTT$lm_logE_logC_r2_adj, highAssociationCutOff), ]
        
        # Bind CNA-low genes to the pan-cancer data frame
        currentDataCNALow <- rbind(currentDataCNALow, currentDataCNALowTT)
        
      }
    }
    
    ## All genes
    # Fit linear model model (methylation Rsq_adj ~ MethSig variance)
    fit <- lm(currentData$gam_logE_M_r2_adj ~ currentData$meth_mean_variance)
    
    # Predict data
    sef <- predict(fit, se.fit = TRUE, interval = "confidence")
    
    # Make data frame for predicted data and 95% confidence interval
    d <- cbind(pred = sef$fit[, "fit"],
               se_up = sef$fit[, "lwr"],
               se_dn = sef$fit[, "upr"],
               vals = currentData$meth_mean_variance)
    
    # Order data
    d <- d[order(d[, "vals"]), ]
    
    ## CNA-high genes
    # Fit linear model model (methylation Rsq_adj ~ MethSig variance)
    fitCNAHigh <- lm(currentDataCNAHigh$gam_logE_M_r2_adj ~ currentDataCNAHigh$meth_mean_variance)
    
    # Predict data
    sefCNAHigh <- predict(fitCNAHigh, se.fit = TRUE, interval = "confidence")
    
    # Make data frame for predicted data and 95% confidence interval
    dCNAHigh <- cbind(pred = sefCNAHigh$fit[, "fit"],
                       se_up = sefCNAHigh$fit[, "lwr"],
                       se_dn = sefCNAHigh$fit[, "upr"],
                       vals = currentDataCNAHigh$meth_mean_variance)
    
    # Order data
    dCNAHigh <- dCNAHigh[order(dCNAHigh[, "vals"]), ]
    
    
    ## CNA-low genes
    # Fit linear model model (methylation Rsq_adj ~ MethSig variance)
    fitCNALow <- lm(currentDataCNALow$gam_logE_M_r2_adj ~ currentDataCNALow$meth_mean_variance)
    
    # Make data frame for predicted data and 95% confidence interval
    sefCNALow <- predict(fitCNALow, se.fit = TRUE, interval = "confidence")
    
    # Make data frame for predicted data and 95% confidence interval
    dCNALow <- cbind(pred = sefCNALow$fit[, "fit"],
                      se_up = sefCNALow$fit[, "lwr"],
                      se_dn = sefCNALow$fit[, "upr"],
                      vals = currentDataCNALow$meth_mean_variance)
    
    # Order data
    dCNALow <- dCNALow[order(dCNALow[, "vals"]), ]
    
    # If not pan-cancer
    if(i != "Pan-cancer"){
      
      # Enter mean MethSig variance for all genes into methVarMeans data frame
      methVarMeans$values[methVarMeans$TumorType == i & methVarMeans$class == "All"] <- mean(currentData$meth_mean_variance)
      
      # Enter mean MethSig variance for CNA-high genes into methVarMeans data frame
      methVarMeans$values[methVarMeans$TumorType == i & methVarMeans$class == "cnaHigh"] <- mean(currentDataCNAHigh$meth_mean_variance)
      
      # Enter mean MethSig variance for CNA-low genes into methVarMeans data frame
      methVarMeans$values[methVarMeans$TumorType == i & methVarMeans$class == "cnaLow"] <- mean(currentDataCNALow$meth_mean_variance)
    }
    
    # Make plot
    plot(currentData$gam_logE_M_r2_adj ~ currentData$meth_mean_variance,
         pch = 16,
         font.main = 1,
         col = adjustcolor(methylationColor, alpha.f = 0.25),
         xlab = "",
         ylab = "",
         las = 1,
         main = i,
         cex = 0.5,
         lwd = 0.75,
         bty = "n",
         ylim = c(0,1),
         xlim = c(0, 2))
    
    ## Add annotation dots to top of plot
    # Data frame for dot locations
    pointsDF <- data.frame(x = c(0.1, 0.8, 1.5), y = 0.875)
    
    # Plot points
    points(y ~ x,
          data = pointsDF,
          pch = 16,
          cex = 2,
          col = c("black", cnaColor, cols[4]))
    
    ## Linear regression for all genes
    # Plot 95% confidence interval
    polygon(x = c(d[, "vals"],
                  d[order(d[, "vals"], decreasing = TRUE), "vals"]),
            y = c(d[, "se_up"],
                  d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
            col = adjustcolor("black", alpha.f = 0.15),
            border = NA)
    
    # Plot predicted line
    lines(pred ~ vals,
          data = d,
          lwd = 2,
          col = "black")
    
    # Find Rsq for model
    rsq <- round(summary(fit)$r.squared, 2)
    
    # Annotate Rsq
    text(x = 0.125,
        y = 0.95,
        pos = 4,
        labels = bquote(italic(R^2)==.(rsq)))
    
    # Find p-value for model 
    if(lmp(fit) < 0.001){p <- " < 0.001"}
    if(!lmp(fit) < 0.001){p <- paste0(" = ", round(lmp(fit), 3))}
    
    # Add p-value
    text(x = 0.125,
         y = 0.8,
         pos = 4,
         labels = bquote(italic(P)~.(p)))
    
    
    ## Linear regression for CNA-high genes
    # Plot 95% confidence interval
    polygon(x = c(dCNAHigh[, "vals"],
                  dCNAHigh[order(dCNAHigh[, "vals"], decreasing = TRUE), "vals"]),
            y = c(dCNAHigh[, "se_up"],
                  dCNAHigh[order(dCNAHigh[, "vals"], decreasing = TRUE), "se_dn"]),
            col = adjustcolor("black", alpha.f = 0.15),
            border = NA)
    
    # Plot predicted line
    lines(pred ~ vals,
          data = dCNAHigh,
          lwd = 2,
          col = cnaColor)
    
    # Find Rsq for model
    rsq <- round(summary(fitCNAHigh)$r.squared, 2)
    
    # Annotate Rsq
    text(x = 0.825,
         y = 0.95,
         pos = 4,
         labels = bquote(italic(R^2)==.(rsq)))
    
    # Find p-value for model 
    if(lmp(fitCNAHigh) < 0.001){p <- " < 0.001"}
    if(!lmp(fitCNAHigh) < 0.001){p <- paste0(" = ", round(lmp(fitCNAHigh), 3))}
    
    # Add p-value
    text(x = 0.825,
         y = 0.8,
         pos = 4,
         labels = bquote(italic(P)~.(p)))
    
    
    ## Linear regression for CNA-low genes
    # Plot 95% confidence interval
    polygon(x = c(dCNALow[, "vals"],
                  dCNALow[order(dCNALow[, "vals"], decreasing = TRUE), "vals"]),
            y = c(dCNALow[, "se_up"],
                  dCNALow[order(dCNALow[, "vals"], decreasing = TRUE), "se_dn"]),
            col = adjustcolor("black", alpha.f = 0.15),
            border = NA)
    
    # Plot predicted line
    lines(pred ~ vals,
          data = dCNALow,
          lwd = 2,
          col = cols[4])
    
    # Find Rsq for model
    rsq <- round(summary(fitCNALow)$r.squared, 2)
    
    # Annotate Rsq
    text(x = 1.525,
         y = 0.95,
         pos = 4,
         labels = bquote(italic(R^2)==.(rsq)))
    
    # Find p-value for model 
    if(lmp(fitCNALow) < 0.001){p <- " < 0.001"}
    if(!lmp(fitCNALow) < 0.001){p <- paste0(" = ", round(lmp(fitCNALow), 3))}
    
    # Add p-value
    text(x = 1.525,
         y = 0.8,
         pos = 4,
         labels = bquote(italic(P)~.(p)))
    
  }

# End figure
dev.off()


##################################################
## To what extent is Methylation_Rsq_adj determined by methylation variance? - minified

# For each tumor type + pan-cancer
for(i in c(tts, "Pan-cancer")){
  
  # Start figure
  pdf(paste0("./Output/Figures/MinifiedMeth-Variance_vs_Meth-Rsquare_",
             i,
             "_",
             w_BP_short,
             "kWindow.pdf"),
      pointsize = 8,
      width = 1.7,
      height = 2,
      useDingbats = FALSE)
  
  # If not pan-cancer
  if(i != "Pan-cancer"){
    
    # Subset data for current tumor type
    currentData <- models[models$TumorType == i, ]
    
    # Subset models with high E-C association
    currentDataCNAHigh <- currentData[currentData$lm_logE_logC_r2_adj > quantile(currentData$lm_logE_logC_r2_adj, highAssociationCutOff), ]
    
    # Subset models with low E-C association
    currentDataCNALow <- currentData[!currentData$lm_logE_logC_r2_adj > quantile(currentData$lm_logE_logC_r2_adj, highAssociationCutOff), ]
    
  }
  
  # If pan-cancer
  if(i == "Pan-cancer"){
    
    # Set currentData as all models
    currentData <- models
    
    # Make an empty data frame for CNA-high genes
    currentDataCNAHigh <- currentData[NULL, ]
    
    # Make an empty data frame for CNA-low genes
    currentDataCNALow <- currentData[NULL, ]
    
    # Get CNA-high/CNA-low genes for each tumor type
    for(j in tts){
      
      # Subset data for the tumor type
      currentDataTT <- currentData[currentData$TumorType == j, ]
      
      # Subset CNA-high genes for the tumor type
      currentDataCNAHighTT <- currentDataTT[currentDataTT$lm_logE_logC_r2_adj > quantile(currentDataTT$lm_logE_logC_r2_adj, highAssociationCutOff), ]
      
      # Bind CNA-high genes to the pan-cancer data frame
      currentDataCNAHigh <- rbind(currentDataCNAHigh, currentDataCNAHighTT)
      
      # Subset CNA-low genes for the tumor type
      currentDataCNALowTT <- currentDataTT[!currentDataTT$lm_logE_logC_r2_adj > quantile(currentDataTT$lm_logE_logC_r2_adj, highAssociationCutOff), ]
      
      # Bind CNA-low genes to the pan-cancer data frame
      currentDataCNALow <- rbind(currentDataCNALow, currentDataCNALowTT)
      
    }
  }
  
  
  ## All genes
  # Fit linear model model (methylation Rsq_adj ~ MethSig variance)
  fit <- lm(currentData$gam_logE_M_r2_adj ~ currentData$meth_mean_variance)
  
  # Predict data
  sef <- predict(fit, se.fit = TRUE, interval = "confidence")
  
  # Make data frame for predicted data and 95% confidence interval
  d <- cbind(pred = sef$fit[, "fit"],
             se_up = sef$fit[, "lwr"],
             se_dn = sef$fit[, "upr"],
             vals = currentData$meth_mean_variance)
  
  # Order data
  d <- d[order(d[, "vals"]), ]
  
  ## CNA-high genes
  # Fit linear model model (methylation Rsq_adj ~ MethSig variance)
  fitCNAHigh <- lm(currentDataCNAHigh$gam_logE_M_r2_adj ~ currentDataCNAHigh$meth_mean_variance)
  
  # Predict data
  sefCNAHigh <- predict(fitCNAHigh, se.fit = TRUE, interval = "confidence")
  
  # Make data frame for predicted data and 95% confidence interval
  dCNAHigh <- cbind(pred = sefCNAHigh$fit[, "fit"],
                    se_up = sefCNAHigh$fit[, "lwr"],
                    se_dn = sefCNAHigh$fit[, "upr"],
                    vals = currentDataCNAHigh$meth_mean_variance)
  
  # Order data
  dCNAHigh <- dCNAHigh[order(dCNAHigh[, "vals"]), ]
  
  
  ## CNA-low genes
  # Fit linear model model (methylation Rsq_adj ~ MethSig variance)
  fitCNALow <- lm(currentDataCNALow$gam_logE_M_r2_adj ~ currentDataCNALow$meth_mean_variance)
  
  # Predict data
  sefCNALow <- predict(fitCNALow, se.fit = TRUE, interval = "confidence")
  
  # Make data frame for predicted data and 95% confidence interval
  dCNALow <- cbind(pred = sefCNALow$fit[, "fit"],
                   se_up = sefCNALow$fit[, "lwr"],
                   se_dn = sefCNALow$fit[, "upr"],
                   vals = currentDataCNALow$meth_mean_variance)
  
  # Order data
  dCNALow <- dCNALow[order(dCNALow[, "vals"]), ]

  # Make plot
  plot(currentData$gam_logE_M_r2_adj[1] ~ currentData$meth_mean_variance[1],
       pch = 16,
       font.main = 1,
       col = "white",
       xlab = "",
       ylab = "",
       main = "",
       las = 1,
       cex = 0.5,
       lwd = 0.75,
       ylim = c(0,0.5),
       xlim = c(0, 1))
  
  ## All genes
  # Plot 95% confidence interval
  polygon(x = c(d[, "vals"],
                d[order(d[, "vals"], decreasing = TRUE), "vals"]),
          y = c(d[, "se_up"],
                d[order(d[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Plot predicted line
  lines(pred ~ vals,
        data = d,
        lwd = 2,
        col = "black")
  
  ## CNA-high genes
  # Plot 95% confidence interval
  polygon(x = c(dCNAHigh[, "vals"],
                dCNAHigh[order(dCNAHigh[, "vals"], decreasing = TRUE), "vals"]),
          y = c(dCNAHigh[, "se_up"],
                dCNAHigh[order(dCNAHigh[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Plot predicted line
  lines(pred ~ vals,
        data = dCNAHigh,
        lwd = 2,
        col = cnaColor)
  
  ## CNA-low genes
  # Plot 95% confidence interval
  polygon(x = c(dCNALow[, "vals"],
                dCNALow[order(dCNALow[, "vals"], decreasing = TRUE), "vals"]),
          y = c(dCNALow[, "se_up"],
                dCNALow[order(dCNALow[, "vals"], decreasing = TRUE), "se_dn"]),
          col = adjustcolor("black", alpha.f = 0.15),
          border = NA)
  
  # Plot predicted line
  lines(pred ~ vals,
        data = dCNALow,
        lwd = 2,
        col = cols[4])
  
  # End figure
  dev.off()
}



