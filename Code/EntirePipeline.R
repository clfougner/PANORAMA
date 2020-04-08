##################################################
## Entire pipeline
##################################################

# First, set the working directory to /path/to/GeneExpressionModulation/

# 1) Get an overview of the available data for each data type in the TCGA cohort
source("./Code/SampleOverview.R")

# 2) Format gene expression data, filter out genes with low/no variance
source("./Code/FormatExpression.R")

# 3) Calculate GII for each sample
source("./Code/Calculate_GII.R") 

# 4) Methylation data to principal components
source("./Code/MethylationPCA.r")

# 5) Plot the variance captured by the principal components
source("./Code/MethylationPCA_Metrics.R")

# 6) Make models
source("./Code/Models.R")

# 7) Gather models into .rds objects
source("./Code/GatherModels.R")

# 8) Run models with all samples
source("./Code/Models.R")

# 9) Visualize results from models
source("./Code/PlotModels.R")