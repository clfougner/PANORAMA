##################################################
# Principal component analysis (PCA) of methylation data
##################################################

## Set the working directory
setwd("/path/to/PANORAMA/")

##################################################
## Load/install libraries
if (!requireNamespace("doMC", quietly = TRUE))
  install.packages("doMC")
library(doMC)

if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
library(data.table)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

if (!requireNamespace("impute", quietly = TRUE))
  BiocManager::install("impute")
library(impute)

if (!requireNamespace("gtools", quietly = TRUE))
  BiocManager::install("gtools")
library(gtools)

if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)

##################################################
## Parameters

# Register number of cores for running in parallel
# Note that large methylation files can make memory the limiting factor here, therefore certain functions should only be run on few cores
coresHigh <- 20
coresLow <- 3

# Number of base pairs around gene (window) for which a CpG should be classified as "belonging" to a gene
windowBP <- 50000
w_BP_short <- windowBP/1000

# Proportion of samples allowed to be NA for a given probe. Remove the probe if more than this number of probes is NA.
removeMissingProbes <- 0.1


##################################################
## Prepare annotations for Illumina 450k array data

# Create a data frame with annotation data
annot <- data.frame(Chrom = Locations@listData$chr,
                  Start = Locations@listData$pos,
                  End = Locations@listData$pos,
                  Probe = Locations@rownames,
                  row.names = Locations@rownames)

# Remove "Chr" so as to match the reference file with protein coding gene positions (./ReferenceFiles/GenePositionsProteinCoding.bed)
annot$Chrom <- gsub(x = annot$Chrom, pattern = "chr", replacement = "")

# Write probes to file (required for bedtools)
write.table(x = annot,
            file = "./ReferenceFiles/Illumina450kProbes.bed",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Load gene positions, downloaded from Ensembl Biomart
genePositions <- read.table(file = "./ReferenceFiles/GenePositionsProteinCoding.bed", header = TRUE)

# Add the window size to the gene start and end positions
genePositions$Start <- genePositions$Start - windowBP
genePositions$End <- genePositions$End + windowBP

# If any start positions are set to a negative coordinate, set as 1
genePositions$Start[genePositions$Start < 0] <- 1

# Write table with gene positions to file (for bedtools)
write.table(x = genePositions,
            file = paste0("./ReferenceFiles/GenePositions_", w_BP_short, "kWindow.bed"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Intersect annotation and probe positions with bedtools
system(paste0("bedtools intersect -b ./ReferenceFiles/Illumina450kProbes.bed -a ./ReferenceFiles/GenePositions_", w_BP_short, "kWindow.bed -wa -wb > ./ReferenceFiles/Illumina450k_GeneAnnot_", w_BP_short,"kWindow.txt"))

# Read output from bedtools intersect
intersectProbes <- read.table(file = paste0("./ReferenceFiles/Illumina450k_GeneAnnot_", w_BP_short, "kWindow.txt"), header = FALSE, stringsAsFactors = FALSE)

# Subset relevant data
intersectProbes <- data.frame(Gene = intersectProbes$V4, Probe = intersectProbes$V8, stringsAsFactors = FALSE)

# Remove non-cg probes (what are actually the non-cg probes?)
intersectProbes <- intersectProbes[grep(x = intersectProbes$Probe, pattern = "cg"), ]

# Create list of genes with associated probes
probeList <- list()

# For each gene
for(i in 1:length(unique(intersectProbes$Gene))){
  
  # Find the probes associated with the given gene
  probes <- intersectProbes$Probe[intersectProbes$Gene == unique(intersectProbes$Gene)[i]]
  
  # Add an entry to the list, containing all probes associated with the given gene
  probeList[[length(probeList) + 1]] <- probes
}

# Name each item in the list with the correct gene name
names(probeList) <- unique(intersectProbes$Gene)

# Save probes
saveRDS(probeList, file = paste0("./ReferenceFiles/Illumina450k_", w_BP_short, "kWindow_probeList.rds"))

##################################################
## Run a per-gene PCA on the methylation data (PCAs trained per tumor type)

# Restore probes
probeList <- readRDS(file = paste0("./ReferenceFiles/Illumina450k_", w_BP_short, "kWindow_probeList.rds"))

# Read methylation beta values
meth <- fread(file = "./Data/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv", verbose = TRUE)

# Remove probes with too many NA values
meth <- meth[rowSums(is.na(meth))/length(meth) < removeMissingProbes, ]

# Sample names in methylation data
methNames <- colnames(meth)
methNames <- methNames[-1]

# Read patientData 
patientData <- read.table(file = "./Output/SampleOverview.txt", header = TRUE, stringsAsFactors = FALSE)

# Get tumor types
tumorTypes <- unique(patientData$TumorType)

# Find all tumors in each tumor type and which have data for all three data types
tumorList <- list()
for(i in tumorTypes){
  methSamplesInTumorType <- patientData$MethID[patientData$TumorType == i & patientData$inAll == 1]
  tumorList[[length(tumorList) + 1]] <- methSamplesInTumorType
}

# Impute missing values
registerDoMC(coresLow)

foreach(i = 1:length(tumorList)) %dopar% {
  tums <- tumorList[[i]]
  currentMethImp <- impute.knn(data = as.matrix(meth[, ..tums]), k = 10)$data
  write.table(x = currentMethImp,
              file = paste0("./Output/TempFiles/MethylationImputedByTT_", tumorTypes[i], ".txt"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}


# Make list of CpGs
methCpG <- meth$V1

# Write CpGs to file
saveRDS(methCpG, file = paste0("./ReferenceFiles/methCpG.rds"))

# Restore CpGs
methCpG <- readRDS(file = "./ReferenceFiles/methCpG.rds")

# Remove methylation input from memory
rm(meth)

# Split genes into equal chunks
numGenes <- 1:length(probeList)
numGenesSplit <- ceiling(length(probeList)/coresHigh)
spl <- split(numGenes, ceiling(seq_along(numGenes)/numGenesSplit))

# Get files with imputed methylation values
fls <- list.files("./Output/TempFiles/", pattern = "MethylationImputedByTT_")

# Register appropriate number of cores to run the PCA
registerDoMC(coresHigh)

# i = for each file
for(i in 1:length(fls)){
  
  # Read methylation beta values
  currentMeth <- read.table(file = paste0("./Output/TempFiles/", fls[i]), header = TRUE)
  
  # Set row names to the CpG
  row.names(currentMeth) <- methCpG
  
  # Get the tumor type for the given file
  tt <- gsub(x = fls[i], pattern = "MethylationImputedByTT_", replacement = "")
  tt <- gsub(x = tt, pattern = ".txt", replacement = "")
  
  # j = for each chunk (one chunk per core)
  foreach(j = 1:length(spl)) %dopar% {
    
    # Create an empty list to save the variance captured by each principal component
    propVariance <- list()
    
    # Create an empty list to save the cumulative variance captured by each principal component
    propVariance_cumulative <- list()
   
    # Create an empty list to save all PCA data
    currentSplitPCA_all <- list()
    
    # k = for each gene
    for(k in 1:length(spl[[j]])){
      
      # Get the location in the probelist with CpGs for the given gene
      idx <- spl[[j]][k]
      
      # Find rows for the CpGs associated with the given gene
      currentRows <- methCpG %in% probeList[[idx]]
      
      # If there are no CpGs associated with the current gene
      if(nrow(currentMeth[currentRows, ]) == 0){
        
        # Set the PCA values and the variance to NA
        propVariance[[length(propVariance) + 1]] <- NA
        propVariance_cumulative[[length(propVariance_cumulative) + 1]] <- NA
        currentSplitPCA_all[[length(currentSplitPCA_all) + 1]] <- NA
      
      # If there is only one CpG associated with the given gene
      } else if(nrow(currentMeth[currentRows, ]) == 1){
        
        # Run a PCA (prcomp)
        currentPCA <- prcomp(x = t(currentMeth[currentRows, ]), scale = FALSE, center = TRUE)
        
        # Insert the values for all PCs into the list
        currentSplitPCA_all[[length(currentSplitPCA_all) + 1]] <- t(currentPCA$x)
        
        # Save the variance for the given principal components
        tmp <- summary(currentPCA)[["importance"]]["Proportion of Variance", ]
        names(tmp) <- "PC1"
        propVariance[[length(propVariance) + 1]] <- tmp
        
        # Save the cumulative variance for the given principal components
        tmp <- summary(currentPCA)[["importance"]]["Cumulative Proportion", ]
        names(tmp) <- "PC1"
        propVariance_cumulative[[length(propVariance_cumulative) + 1]] <- tmp
        
        # If there is more than one CpG for the given gene
        } else if(nrow(currentMeth[currentRows, ]) > 1){
        
        # Run a PCA (prcomp)
        currentPCA <- prcomp(x = t(currentMeth[currentRows, ]), scale = FALSE, center = TRUE)
          
        # Insert the values for all PCs into the list
        currentSplitPCA_all[[length(currentSplitPCA_all) + 1]] <- t(currentPCA$x)
        
        # Save the variance for the given principal components
        propVariance[[length(propVariance) + 1]] <- summary(currentPCA)[["importance"]]["Proportion of Variance", ]
        
        # Save the cumulative variance for the given principal components
        propVariance_cumulative[[length(propVariance_cumulative) + 1]] <- summary(currentPCA)[["importance"]]["Cumulative Proportion", ]
      }
    }

    # Asssign each member of the variance list with the correct gene name
    names(propVariance) <- names(probeList[spl[[j]]])
    
    # Asssign each member of the cumulative variance list with the correct gene name
    names(propVariance_cumulative) <- names(probeList[spl[[j]]])
    
    # Asssign each member of the full PCA list  with the correct gene name
    names(currentSplitPCA_all) <- names(probeList[spl[[j]]])
    
     # Write the list with all PCA values to files
    saveRDS(currentSplitPCA_all, file = paste0("./Output/TempFiles/MethylationList_AllPCs_", w_BP_short, "kWindow_", tt, "_", j, ".rds"))
    
    # Write the variances for the chunk to file
    saveRDS(propVariance, file = paste0("./Output/TempFiles/VarianceAllPCs_individual_", w_BP_short, "kWindow_", tt, "_", j, ".rds"))
    
    # Write the cumulative variances for the chunk to file
    saveRDS(propVariance_cumulative, file = paste0("./Output/TempFiles/VarianceAllPCs_cumulative_", w_BP_short, "kWindow_", tt, "_", j, ".rds"))
  }
}

## Collect the outputs from above into a single file per tumor type
foreach(i = 1:length(tumorTypes)) %dopar% {
    
  ## PC values
  
  # Find files
  filesAllPC <- list.files("./Output/TempFiles/", pattern = paste0("MethylationList_AllPCs_", w_BP_short, "kWindow_", tumorTypes[i], "_"))
  filesAllPC <- mixedsort(filesAllPC)
  
  # For each file
  for(j in 1:length(filesAllPC)){
    
    # Read the data 
    currentAllPCs <- readRDS(file = paste0("./Output/TempFiles/", filesAllPC[j]))
    
    # If first file
     if(j == 1){
       # Make a new object
       gatheredAllPCs <- currentAllPCs
     }
    
     # If last file
     if(j == length(filesAllPC)){
       # Bind to the aggregated object
       gatheredAllPCs <- c(gatheredAllPCs, currentAllPCs)
       # Write aggregated object to file
       saveRDS(gatheredAllPCs, paste0("./Output/MethylationAllPCValues_", w_BP_short, "kWindow_", tumorTypes[i], ".rds"))
     }
    
     # If not first or last file
     if(j != 1 & j != length(filesAllPC)){
       # Bind current data to the aggregated object
       gatheredAllPCs <- c(gatheredAllPCs, currentAllPCs)
     }
   }
   
   ## Proportion of variance
  
   # Find files
  filesPropVariance <- list.files("./Output/TempFiles/", pattern = paste0("VarianceAllPCs_individual_", w_BP_short, "kWindow_", tumorTypes[i], "_"))
  filesPropVariance <- mixedsort(filesPropVariance)
  
  # For each file
  for(j in 1:length(filesPropVariance)){
    
    # Read the data
    currentAllVariances <- readRDS(file = paste0("./Output/TempFiles/", filesPropVariance[j]))
    
    # If first file
    if(j == 1){
      # Make a new object
      gatheredAllVariances <- currentAllVariances
    }
    
    # If last file
    if(j == length(filesPropVariance)){
      # Bind to aggregated list
      gatheredAllVariances <- c(gatheredAllVariances, currentAllVariances)
      # Write file
      saveRDS(gatheredAllVariances, paste0("./Output/MethylationAllVariances_", w_BP_short, "kWindow_", tumorTypes[i], ".rds"))
    }
    
    # If not first or last file
    if(j != 1 & j != length(filesPropVariance)){
      # Bind to aggregated list
      gatheredAllVariances <- c(gatheredAllVariances, currentAllVariances)
    }
  }
  
  ## Proportion of variance - cumulative
  
  # Find files
  filesPropVariance_cumulative <- list.files("./Output/TempFiles/", pattern = paste0("VarianceAllPCs_cumulative_", w_BP_short, "kWindow_", tumorTypes[i], "_"))
  filesPropVariance_cumulative <- mixedsort(filesPropVariance_cumulative)
  
  # For each file
  for(j in 1:length(filesPropVariance_cumulative)){
    
    # Read the data
    currentAllVariances_cumulative <- readRDS(file = paste0("./Output/TempFiles/", filesPropVariance_cumulative[j]))
    
    # If first file
    if(j == 1){
      # Make a new object
      gatheredAllVariances_cumulative <- currentAllVariances_cumulative
    }
    
    # If last file
    if(j == length(filesPropVariance_cumulative)){
      # Bind to aggregated list
      gatheredAllVariances_cumulative <- c(gatheredAllVariances_cumulative, currentAllVariances_cumulative)
      # Write file
      saveRDS(gatheredAllVariances_cumulative, paste0("./Output/MethylationAllCumulativeVariances_", w_BP_short, "kWindow_", tumorTypes[i], ".rds"))
    }
    
    # If not first or last file
    if(j != 1 & j != length(filesPropVariance_cumulative)){
      # Bind to aggregated list
      gatheredAllVariances_cumulative <- c(gatheredAllVariances_cumulative, currentAllVariances_cumulative)
    }
  }
}
