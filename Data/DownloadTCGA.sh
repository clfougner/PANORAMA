#!/bin/bash

# Set path to /path/to/GeneExpressionModulation/
# Change directory to /Data/
cd ./Data/

# Gene expression (RNA)
## RNA batch corrected matrix
wget https://api.gdc.cancer.gov/data/9a4679c3-855d-4055-8be9-3577ce10f66e
mv 9a4679c3-855d-4055-8be9-3577ce10f66e EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv

# Copy number
## ISAR-corrected SNP6 whitelisted copy number segments file
wget https://api.gdc.cancer.gov/data/74abec30-668d-45ef-b1e3-ff49e8f132c7
mv 74abec30-668d-45ef-b1e3-ff49e8f132c7 ISAR_corrected.PANCAN_Genome_Wide_SNP_6_whitelisted.seg

## ABSOLUTE purity/ploidy file
wget https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5
mv 4f277128-f793-4354-a13d-30cc7fe9f6b5 TCGA_mastercalls.abs_tables_JSedit.fixed.txt

## ABSOLUTE-annotated seg file
wget https://api.gdc.cancer.gov/data/0f4f5701-7b61-41ae-bda9-2805d1ca9781
mv 0f4f5701-7b61-41ae-bda9-2805d1ca9781 TCGA_mastercalls.abs_segtabs.fixed.txt

## Aneuploidy scores and arm calls file
wget https://api.gdc.cancer.gov/data/4c35f34f-b0f3-4891-8794-4840dd748aad
mv 4c35f34f-b0f3-4891-8794-4840dd748aad PANCAN_ArmCallsAndAneuploidyScore_092817.txt

# Methylation
## DNA methylation 450k only beta value data matrix
wget https://api.gdc.cancer.gov/data/99b0c493-9e94-4d99-af9f-151e46bab989
mv 0f4f5701-7b61-41ae-bda9-2805d1ca9781 jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv

# Mutation
## Public mutation annotation file
wget https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc
mv 1c8cfe5f-e52d-41ba-94da-f15ea1337efc mc3.v0.2.8.PUBLIC.maf.gz

## Mutation load
wget https://api.gdc.cancer.gov/data/ff3f962c-3573-44ae-a8f4-e5ac0aea64b6
mv ff3f962c-3573-44ae-a8f4-e5ac0aea64b6 mutation-load-updated.txt

# miRNA
## miRNA batch corrected matrix
wget https://api.gdc.cancer.gov/data/1c6174d9-8ffb-466e-b5ee-07b204c15cf8
mv 1c6174d9-8ffb-466e-b5ee-07b204c15cf8 pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv

## miRNA sample information
wget https://api.gdc.cancer.gov/data/55d9bf6f-0712-4315-b588-e6f8e295018e
mv 55d9bf6f-0712-4315-b588-e6f8e295018e PanCanAtlas_miRNA_sample_information_list.txt

# Protein
## RPPA batch corrected matrix
wget https://api.gdc.cancer.gov/data/fcbb373e-28d4-4818-92f3-601ede3da5e1
mv fcbb373e-28d4-4818-92f3-601ede3da5e1 TCGA-RPPA-pancan-clean.txt

# Clinical (https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018)
## TCGA-CDR
wget https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
mv 1b5f413e-a8d1-4d10-92eb-7c4ae739ed81 TCGA-CDR-SupplementalTableS1.xlsx
