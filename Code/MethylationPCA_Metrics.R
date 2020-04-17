##################################################
# Determine the variance captured by the principal component analyses
##################################################

## Set the working directory
setwd("/path/to/PANORAMA/")


##################################################
## Parameters
# Number of base pairs around gene (window) for which a CpG should be classified as "belonging" to a gene
windowBP <- 50000
w_BP_short <- windowBP/1000

# Minimum number of tumors for a tumor type to be included
minTums <- 100

# Number of PCs to plot
plotPCs <- 10

# Number of cores to use
crs <- 20
print(paste0(crs, " cores"))


##################################################
## Load and install packaes
if (!requireNamespace("gtools", quietly = TRUE))
  BiocManager::install("gtools")
library(gtools)

if (!requireNamespace("doMC", quietly = TRUE))
  install.packages("doMC")
library(doMC)
registerDoMC(crs)

##################################################
## Prepare sample overview data

# Load patient overview
patientData <- read.table(file = "./Output/SampleOverview.txt",
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE)

# Empty vector for tumor types with minimum number of tumors
ttsMinTums <- c()

# Find number of tumors in each tumor type with all data required
for(i in unique(patientData$TumorType)){
  ttsMinTums <- c(ttsMinTums, nrow(patientData[patientData$TumorType == i & patientData$inAll == 1, ]))
}

# Name the vector
names(ttsMinTums) <- unique(patientData$TumorType)

# Remove tumor types with too few samples
ttsMinTums <- ttsMinTums[ttsMinTums >= minTums]
 

##################################################
## Convert data from list to long format

# Files with cumulative variance captured by each principal component for each gene
cumulativeVarianceFiles <- list.files("./Output/", pattern =  paste0("MethylationAllCumulativeVariances_", w_BP_short, "kWindow_"))

# For each file
foreach (i = 1:length(cumulativeVarianceFiles)) %dopar% {
  
  # Read the input data
  currentData <- readRDS(paste0("./Output/", cumulativeVarianceFiles[i]))
  
  # Create an empty data frame to place the data in to
  currentLong <- data.frame(matrix(data = NA, ncol = 3, nrow = sum(lengths(currentData)), dimnames = list(NULL, c("Gene", "PC", "CumulativeVariance"))))
  
  # Make vector of gene names
  geneNames <- c()
  for(j in 1:length(currentData)){
    geneNames <- c(geneNames, rep(names(currentData[j]), times = lengths(currentData[j])))
  }
  
  # Set gene names
  currentLong$Gene <- geneNames
  
  # Enter data for each gene
  for(j in 1:length(currentData)){
    currentLong$CumulativeVariance[currentLong$Gene == names(currentData[j])] <- currentData[[j]]
    if(!is.null(names(currentData[[j]][1]))){
      currentLong$PC[currentLong$Gene == names(currentData[j])] <- names(currentData[[j]])  
    }
  }
  
  # Get the tumor type for the current file
  tt <- gsub(x = cumulativeVarianceFiles[i], pattern = paste0("MethylationAllCumulativeVariances_", w_BP_short, "kWindow_"), replacement = "")
  
  # Save the long data frame to file
  saveRDS(currentLong, paste0("./Output/TempFiles/Methylation_long_AllCumulativeVariances_", w_BP_short, "kWindow_", tt))
}

##################################################
## Calculate the median variances

# Get the files output from the previous sections
files <- paste0("./Output/TempFiles/Methylation_long_AllCumulativeVariances_", w_BP_short, "kWindow_", names(ttsMinTums), ".rds")

# Aggregate files into one
varsFirst <- readRDS(file = files[1])

# Set the principal components as factors in the corrrect order
varsFirst$PC <- factor(varsFirst$PC, levels = mixedsort(unique(varsFirst$PC)))

# Subset to only the wanted PCs
varsFirst <- varsFirst[varsFirst$PC %in% paste0("PC", 1:plotPCs), ]
varsFirst <- droplevels(varsFirst)

# Add tumor type as a new column
varsFirst$TumorType <- names(ttsMinTums)[1]

# For each remaining tumor type
for(i in 2:length(ttsMinTums)){
  
  # Read variance file
  currentVars <- readRDS(file = files[i])
  
  # Set the principal components as factors in the corrrect order
  currentVars$PC <- factor(currentVars$PC, levels = mixedsort(unique(currentVars$PC)))
  
  # Subset to only the wanted PCs
  currentVars <- currentVars[currentVars$PC %in% paste0("PC", 1:plotPCs), ]
  currentVars <- droplevels(currentVars)
  
  # Add tumor type to TT column
  currentVars$TumorType <- names(ttsMinTums)[i]
  
  # Gather data
  varsFirst <- rbind(varsFirst, currentVars)
}

# Add an item for the pan-cancer medians
ttsPlusAll <- c(names(ttsMinTums), "All")

# Make a data frame for the medians
mediansDF <- data.frame(matrix(nrow = length(ttsPlusAll) * plotPCs,
                               ncol = 3,
                               dimnames = list(NULL, c("PC", "MedianVar", "TumorType"))))

# Enter PC values into the data frame
mediansDF$PC <- paste0("PC", 1:plotPCs)

# Enter Tumor type values into the data frame
mediansDF$TumorType <- rep(ttsPlusAll, each = plotPCs)

# For each tumor type
for(i in ttsPlusAll){
  # For each principal component
  for(j in unique(mediansDF$PC)){
    
    # If the "tumor type" is all, calculate medians from mediansDF (i.e. median of medians)
    if(i == "All"){
      currentMedian <- median(mediansDF$MedianVar[mediansDF$PC == j], na.rm = TRUE)
    }
    
    # If the tumor type is not all, calulate median variance captured by all genes
    if(i != "All"){
      currentMedian <- median(varsFirst$CumulativeVariance[varsFirst$TumorType == i & varsFirst$PC == j])  
    }
    
    # Enter into medians data frame
    mediansDF$MedianVar[mediansDF$TumorType == i & mediansDF$PC == j] <- currentMedian
    
  }
}

# Set PC column as a factor
mediansDF$PC <- factor(mediansDF$PC, levels = mixedsort(unique(mediansDF$PC)))

# Save medians to file
write.table(mediansDF, file = paste0("./Output/VarianceCaptured_", w_BP_short, "kWindow.txt"), sep = "\t", quote = FALSE)


##################################################
## Make plots

pdf(paste0("./Output/Figures/Methylation_PCA_VariancesCumulative_boxplot_", w_BP_short, "kWindow.pdf"),
    pointsize = 7,
    height = 11,
    width = 8,
    useDingbats = FALSE)
par(mfrow = c(8, 3), mai = c(0.2, 0.2, 0.4, 0.2)) # bottom, left, top, right

# For each tumor type
for(i in ttsPlusAll){

  # If not "All"
  if(i != "All"){

  # Get the input data for the given tumor type
  currentData <- varsFirst[varsFirst$TumorType == i, ]

  # Make boxplot
  boxplot(CumulativeVariance ~ PC,
          data = currentData,
          boxwex = 0.75,
          las = 1,
          ylim = c(0, 1),
          frame = FALSE,
          outline = FALSE,
          main = i,
          font = 1,
          col = "#4C856E") 
  }
  
  # If "All" (pan-cancer medians)
  if(i == "All"){
    
    # Get the input data for the given tumor type
    currentData <- mediansDF[mediansDF$TumorType != i, ]
    
    # Make boxplot
    boxplot(MedianVar ~ PC,
            data = currentData,
            boxwex = 0.75,
            las = 1,
            ylim = c(0, 1),
            frame = FALSE,
            outline = TRUE,
            main = "Medians",
            font = 1,
            col = "#4C856E",
            pch = 19)
  }
}

dev.off()
