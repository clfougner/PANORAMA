##################################################
# Calculate genomic instability index in each tumor type
##################################################

## Set the working directory
setwd("/path/to/PAMORAMA/")

##################################################
## Load data
# Read ploidy data
ploidyData <- read.table(file = "./Data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         sep = "\t")

# Read segments
segs <- read.table(file = "./Data/ISAR_corrected.PANCAN_Genome_Wide_SNP_6_whitelisted.seg",
                   header = TRUE,
                   stringsAsFactors = FALSE)

# Read sample overview 
patientData <- read.table(file = "./Output/SampleOverview.txt",
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          na.strings = "<NA>")

##################################################
## Format data

# Segment data is in the following format: log2(copy-number/2) (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/). Revert to copy number 
segs$Segment_Mean <- 2^(segs$Segment_Mean) * 2

# Add a column for GII to patientData
patientData$GII <- NA

# Add a column for ploidy corrected GII to patientData
patientData$GII_ploidyCorrected <- NA

# For i = each sample with CNA data
for(i in patientData$CNA_ID[patientData$inCNA == 1]){
  
  # Get the segments for the current sample
  currentSegs <- segs[segs$Sample == i, ]
  
  if(nrow(currentSegs) == 0){next}
  
  # Get the ploidy for the current sample (nearest integer value)
  currentPloidy <- round(ploidyData[ploidyData$sample == i, "ploidy"])
  
  # If there is no value available for ploidy, set ploidy to 2
  if(is.na(currentPloidy)){
    currentPloidy <- 2
  }
    
  # Find the number of nucleotides in all segments
  currentSegsNucleotides <- sum(currentSegs$End - currentSegs$Start)
  

  # Create an empty list for the lengths of all copy number aberrant segments (relative to ploidy)
  currentAberrantPloidyCorrected <- c()
  # For j = each segment
  for(j in 1:nrow(currentSegs)){
    # If the nearest integer value for segment mean is different from the nearest integer value for ploidy
    if(round(currentSegs[j, "Segment_Mean"]) != currentPloidy){
      # Place the current segment length into the aggregated list
      currentAberrantPloidyCorrected <- c(currentAberrantPloidyCorrected,
                                          (currentSegs[j, "End"] - currentSegs[j, "Start"]))
    }
  }
    # Calculate the ploidy corrected GII and place into the patientData dataframe
  patientData[patientData$CNA_ID == i, "GII_ploidyCorrected"] <- round(sum(currentAberrantPloidyCorrected)/currentSegsNucleotides, digits = 3) 
  
  
  # Create an emply list for all copy number aberrant segments (copy number != 2; not ploidy corrected)
  currentAberrant <- c()
  # For j = each segment
  for(j in 1:nrow(currentSegs)){
    # If the nearest integer value for segment mean is not 2
    if(round(currentSegs[j, "Segment_Mean"]) != 2){
      # Place the current segment length into the aggregated list
      currentAberrant <- c(currentAberrant,
                          (currentSegs[j, "End"] - currentSegs[j, "Start"]))
    }
  }
  # Calculate the GII (not ploidy-corrected) and place into the patientData dataframe
  patientData[patientData$CNA_ID == i, "GII"] <- round(sum(currentAberrant)/currentSegsNucleotides, digits = 3)
}


# Save the data frame to file
write.table(x = patientData,
            file = "./Output/SampleOverview.txt",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

