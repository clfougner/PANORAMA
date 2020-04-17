##################################################
# Format gene expression data, and remove genes with no/low variance
##################################################

## Set the working directory
setwd("/path/to/PANORAMA/")

##################################################
## Load/install packages
if (!requireNamespace("ks", quietly = TRUE))
  install.packages("ks")
library(ks)


##################################################
## Load and format data

# Read the gene expression data
exprs <- read.table(file = "./Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv",
                    header = TRUE,
                    stringsAsFactors = FALSE)

# Make list of gene IDs
exprs_gi <- exprs$gene_id

# Remove EntrezIDs, keep gene names
exprs_gi <- gsub(x = exprs_gi, pattern = "\\|\\d+", replacement = "")

# Remove duplicates: this removes rows with no gene name and one row of SLC35E2, which for some reason is duplicated
exprs_new <- exprs[!duplicated(exprs_gi), ]
exprs_gi <- exprs_gi[!duplicated(exprs_gi)]

# Remove gene IDs from the expression data
exprs_new <- exprs_new[, 2:ncol(exprs_new)]

# Set row names to gene names
rownames(exprs_new) <- exprs_gi

# Remove the remaining row with an unknown gene name
exprs_new <- exprs_new[!rownames(exprs_new) == "?", ]

# Log2 transform expression data
exprs_new <- log2(exprs_new + 1)

# Get the standard deviation for each gene
exprs_new <- transform(exprs_new, SD = apply(exprs_new, 1, sd, na.rm = TRUE))

# Generate the data used for the density plot (for standard deviation)
d <- density(exprs_new[, "SD"])

# Find the inflection point for the density plot
h <- hns(exprs_new[, "SD"], deriv.order = 1)

# Generate the plot
pdf("./Output/Figures/FilterExpressionGenes_StandardDeviation.pdf", width = 2.5, height = 2.5, pointsize = 8)
  
plot(d,
    xlab = "Standard deviation",
    y = "Density",
    frame = FALSE,
    main = "",
    las = 1,
    ylim = c(0, 1))

  # Line for the inflection point
  abline(v = h)
  text(x = 4, y = 0.8, labels = paste0("x = ", round(h, digits = 4)))
  
dev.off()

# Remove genes with a standard deviation less than the inflection point
exprs_new <- exprs_new[exprs_new[, "SD"] > h, ]

# Remove standard deviation column
exprs_new <- exprs_new[, 1:(ncol(exprs_new) - 1)]

# Change periods to dashes in sample names
colnames(exprs_new) <- gsub("\\.", "-", colnames(exprs_new))

# Save to file
write.table(exprs_new,
            file = "./Output/GeneExpressionFormatted_sdFiltered.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)
