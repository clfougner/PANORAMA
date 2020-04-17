##################################################
## Make an overview of the samples (file names and tumor types for the different samples and data types)
##################################################

## Set working directory
setwd("/path/to/PANORAMA/")

##################################################
## Load libraries
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
library(data.table)

if (!requireNamespace("readxl", quietly = TRUE))
  install.packages("readxl")
library("readxl")

##################################################
## Read "clinical" data
patientData <- read_excel("./Data/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1)

##################################################
## Methylation

# Read methylation 450k data
meth <- fread(file = "./Data/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv", verbose = TRUE)

# Get sample names in methylation data
methNames <- colnames(meth)
methNames <- methNames[-1]

# Create a data frame for methylation samples
methDF <- data.frame(PatientID = substr(methNames, start = 1, stop = 12),
                     SampleID = substr(methNames, start = 1, stop = 16),
                     AnalysisID = methNames,
                     SampleType = substr(methNames, start = 14, stop = 16),
                     TumorType = NA)

# Find the tumor type for each methylation sample
toRemove <- c()
for(i in 1:nrow(methDF)){
  tt <- patientData$type[patientData$bcr_patient_barcode == methDF[i, "PatientID"]]
  if(length(tt) == 1){methDF[i, "TumorType"] <- tt}
  if(length(tt) == 0){toRemove <- c(toRemove, i)}
}

# Remove samples with unknown tumor type (not found in the CDR)
methDF <- methDF[-toRemove, ]

# Remove non-tumor samples
methDF <- methDF[methDF$SampleType != "11A", ]
methDF <- methDF[methDF$SampleType != "11B", ]

# There are also too few OV tumors for meaningful analysis (n = 10), therefore we remove those tumors
methDF <- methDF[methDF$TumorType != "OV", ]

# Convert levels to character
methDF[] <- lapply(methDF, as.character)

##################################################
## Copy number

segs <- read.table(file = "./Data/ISAR_GISTIC.all_data_by_genes.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   sep = "\t")

cnaNames <- unique(colnames(segs[4:ncol(segs)]))

# Change periods to dashes in names
cnaNames <- gsub("\\.", "-", cnaNames)

# Create a data frame for copy number samples
cnaDF <- data.frame(PatientID = substr(cnaNames, start = 1, stop = 12),
                    SampleID = substr(cnaNames, start = 1, stop = 16),
                    AnalysisID = cnaNames,
                    SampleType = substr(cnaNames, start = 14, stop = 16),
                    TumorType = NA)

# Find the tumor type for each CNA sample
toRemove <- c()
for(i in 1:nrow(cnaDF)){
  tt <- patientData$type[patientData$bcr_patient_barcode == cnaDF[i, "PatientID"]]
  if(length(tt) == 1){cnaDF[i, "TumorType"] <- tt}
  if(length(tt) == 0){toRemove <- c(toRemove, i)}
}

# Remove samples with unknown tumor type (not found in the CDR)
cnaDF <- cnaDF[-toRemove, ]

# Remove OV samples (because too few samples with methylation data)
cnaDF <- cnaDF[cnaDF$TumorType != "OV", ]

# Convert levels to character
cnaDF[] <- lapply(cnaDF, as.character)

##################################################
## Gene expression

# Read data
exprs <- read.table(file = "./Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv", header = TRUE, stringsAsFactors = FALSE)

# Get names
exprs_names <- colnames(exprs)

# Change periods to dashes in names
exprs_names <- gsub("\\.", "-", exprs_names)

# Create a data frame for copy number samples
exprsDF <- data.frame(PatientID = substr(exprs_names, start = 1, stop = 12),
                    SampleID = substr(exprs_names, start = 1, stop = 16),
                    AnalysisID = exprs_names,
                    SampleType = substr(exprs_names, start = 14, stop = 16),
                    TumorType = NA)

# Find the tumor type for each expression sample
toRemove <- c()
for(i in 1:nrow(exprsDF)){
  tt <- patientData$type[patientData$bcr_patient_barcode == exprsDF[i, "PatientID"]]
  if(length(tt) == 1){exprsDF[i, "TumorType"] <- tt}
  if(length(tt) == 0){toRemove <- c(toRemove, i)}
}

# Remove samples with unknown tumor type (not found in the CDR)
exprsDF <- exprsDF[-toRemove, ]

# Remove OV samples (because too few samples with methylation data)
exprsDF <- exprsDF[exprsDF$TumorType != "OV", ]

# Convert levels to character
exprsDF[] <- lapply(exprsDF, as.character)

##################################################
## All levels

# Gather all data
allTumors <- rbind(cnaDF, exprsDF, methDF)

# Remove duplicated
allTumors <- allTumors[!duplicated(allTumors$SampleID), ]

# Create a data frame with each tumor sample as a separate row
allDF <- data.frame(PatientID = allTumors$PatientID,
                    SampleID = allTumors$SampleID,
                    SampleType = allTumors$SampleType,
                    TumorType = allTumors$TumorType,
                    ExprID = NA,
                    MethID = NA,
                    CNA_ID = NA,
                    inExpr = 0,
                    inMeth = 0,
                    inCNA = 0,
                    inExprInMeth = 0,
                    inExprInCNA = 0,
                    inMethInCNA = 0,
                    inAll = 0,
                  stringsAsFactors = FALSE)

# Populate the data frame
for(i in 1:nrow(allDF)){
  
  # Get current sample ID
  currentSample <- allDF$SampleID[i]
  
  # Methylation
  currentMeth <- methDF[methDF$SampleID == currentSample, ]
  if(nrow(currentMeth) == 1){
    allDF[i, "MethID"] <- currentMeth[, "AnalysisID"]
    allDF[i, "inMeth"] <- 1
  }
  
  # Expression
  currentExprs <- exprsDF[exprsDF$SampleID == currentSample, ]
  if(nrow(currentExprs) == 1){
    allDF[i, "ExprID"] <- currentExprs[, "AnalysisID"]
    allDF[i, "inExpr"] <- 1
  }
  
  # Copy number
  currentCNA <- cnaDF[cnaDF$SampleID == currentSample, ]
  if(nrow(currentCNA) == 1){
    allDF[i, "CNA_ID"] <- currentCNA[, "AnalysisID"]
    allDF[i, "inCNA"] <- 1
  }
  
  # Expression and methylation available
  if(allDF[i, "inExpr"] == 1 & allDF[i, "inMeth"] == 1){allDF[i, "inExprInMeth"] <- 1}
  
  # Expression and CNA available
  if(allDF[i, "inExpr"] == 1 & allDF[i, "inCNA"] == 1){allDF[i, "inExprInCNA"] <- 1}
  
  # Methylation and CNA available
  if(allDF[i, "inMeth"] == 1 & allDF[i, "inCNA"] == 1){allDF[i, "inMethInCNA"] <- 1}
    
  # All data levels available
  if(allDF[i, "inMeth"] == 1 & allDF[i, "inCNA"] == 1 & allDF[i, "inExpr"] == 1){allDF[i, "inAll"] <- 1}
}

# Save the data frame to file
write.table(x = allDF, file = "./Output/SampleOverview.txt", row.names = FALSE, quote = FALSE, sep = "\t")
