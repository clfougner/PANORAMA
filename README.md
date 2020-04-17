# A pan-cancer atlas of transcriptional dependence on DNA methylation and copy number aberrations

Fougner, C., Höglander, E.K., Lien, T.G., Sørlie, T.S., Nord, S. & Lingjærde, O.C. A pan-cancer atlas of transcriptional dependence on DNA methylation and copy number aberrations. bioRxiv (2020).

## Data
The [TCGA pan-cancer cohort](https://www.cell.com/pb-assets/consortium/pancanceratlas/pancani3/index.html) is used in this study. To download the necessary data, first install the [GDC data transfer tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool). Once installed, clone this repository and set the path to `/path/to/GeneExpressionModulation/` and run the following to download:

```bash
cd ./Data/
gdc-client download -m TCGAmanifest.txt
```

If for some reason the files can't be downloaded with the GDC client, the individual files can be downloaded from [here](https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin) and [here](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018). This can be done using wget by running the script  `sh ./Data/DownloadTCGA.sh`. This is however not recommended, and takes about ten times longer than using the GDC client.

The files that need to be in `./Data/`to run the analyses are:
```
# Clinical data
TCGA-CDR-SupplementalTableS1.xlsx

# Methylation beta values
jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv

# Copy number segments
ISAR_corrected.PANCAN_Genome_Wide_SNP_6_whitelisted.seg

# Gene-centric copy number values
ISAR_GISTIC.all_data_by_genes.txt

# Gene expression data
EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv

# ABSOLUTE calls
TCGA_mastercalls.abs_tables_JSedit.fixed.txt
```

## Reference files
Reference files are saved under `./ReferenceFiles/`.

The positions for protein coding genes are saved to `./ReferenceFiles/GenePositionsProteinCoding.bed` (included in this repository). The file is originally downloaded from [Ensembl's Biomart](http://www.ensembl.org/biomart/martview/), with the following settings:
![Ensembl Biomart](/ReferenceFiles/GenePositionsProteinCoding.png?raw=true)

## Computing and software
The analyses are run using [R](https://www.r-project.org) on a UNIX-like system (FreeBSD). Large files (~40gb) are read into memory, so a powerful computer/server will be required (we have run analyses on a server with 32 cores and 256gb ram).

* All R packages are loaded, and installed if necessary, in the scripts.
* [bedtools](https://bedtools.readthedocs.io/en/latest/) must be installed so that it can be invoked by typing `bedtools` in the command line/terminal. If this is not possible, the `MethylationPCA.R` script must be updated so that bedtools can be run on your system.

## Analyses
All code is saved to `./Code/`. The entire pipeline is implemented in the file `./Code/EntirePipeline.R`. Briefly, the analyses should be run in the following order:

```r
source("./Code/SampleOverview.R")
source("./Code/FormatExpression.R")
source("./Code/Calculate_GII.R")
source("./Code/MethylationPCA.R")
source("./Code/MethylationPCA_Metrics.R")
source("./Code/Models.R")
source("./Code/GatherModels.R")
source("./Code/Models-allSamples.R")
source("./Code/PlotModels.R")
```

Note that it may be a good idea to restart the R session after running each script.

## Output
All output is saved to the `./Output/` folder. Figures are saved to `./Output/Figures/`. Intermediate/tempporary files are saved to `./Output/TempFiles/`. 
