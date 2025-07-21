# =============================================================================
# IL6 Genetic Instrument Construction
# =============================================================================

# Clean environment and load required libraries
rm(list=ls())
library(data.table)
library(MendelianRandomization)
library(plyr)
library(TwoSampleMR)
library(readxl)
library(dplyr)
library(metafor) 
library(ggplot2)
library(meta)

# =============================================================================
# Load and process CRP GWAS data
# =============================================================================
# Load CRP GWAS meta-analysis data (UKB + CHARGE)
crp <- fread("data/raw/GCST90029070_buildGRCh37.tsv.gz")

# Standardize column names
colnames(crp)[1] <- "unique_id"
colnames(crp)[2] <- "variant_id"
colnames(crp)[3] <- "chromosome"
colnames(crp)[4] <- "bp"
colnames(crp)[5] <- "effect_allele"
colnames(crp)[6] <- "other_allele"
colnames(crp)[7] <- "beta"
colnames(crp)[8] <- "se"
colnames(crp)[9] <- "pvalue"

# =============================================================================
# Define IL6 gene regions and extract variants
# =============================================================================
# Filter for chromosome 7 variants
crp_il6_unf <- subset(crp, chromosome == "7")
crp_il6_unf$pos <- as.numeric(crp_il6_unf$bp)

# IL6 gene location based on GRCh37/hg19 (NCBI): chr7:22,766,819-22,771,617
# Extract variants in three different windows around IL6 gene

# 1. +/- 300kB window (main instrument)
crp_il6_unf <- subset(crp_il6_unf, pos > 22466819 & pos < 23071617) 

# 2. +/- 100kB window (sensitivity analysis)
crp_il6_unf_100 <- subset(crp_il6_unf, pos > 22666819 & pos < 22871617) 

# 3. +/- 10kB window (sensitivity analysis)
crp_il6_unf_10 <- subset(crp_il6_unf, pos > 22756819 & pos < 22781617) 

# =============================================================================
# Prepare data for clumping
# =============================================================================
# Standardize column names for TwoSampleMR compatibility
colnames(crp_il6_unf)[3] <- "chr_name"
colnames(crp_il6_unf)[4] <- "chr_start"
colnames(crp_il6_unf)[2] <- "SNP"
colnames(crp_il6_unf)[9] <- "pval.exposure"

colnames(crp_il6_unf_100)[3] <- "chr_name"
colnames(crp_il6_unf_100)[4] <- "chr_start"
colnames(crp_il6_unf_100)[2] <- "SNP"
colnames(crp_il6_unf_100)[9] <- "pval.exposure"

colnames(crp_il6_unf_10)[3] <- "chr_name"
colnames(crp_il6_unf_10)[4] <- "chr_start"
colnames(crp_il6_unf_10)[2] <- "SNP"
colnames(crp_il6_unf_10)[9] <- "pval.exposure"

# Filter for genome-wide significant variants (p < 5e-8)
crp_il6_unf$pval.exposure <- as.numeric(crp_il6_unf$pval.exposure)
crp_il6_unf <- subset(crp_il6_unf, pval.exposure < 5e-8)

crp_il6_unf_100$pval.exposure <- as.numeric(crp_il6_unf_100$pval.exposure)
crp_il6_unf_100 <- subset(crp_il6_unf_100, pval.exposure < 5e-8)

crp_il6_unf_10$pval.exposure <- as.numeric(crp_il6_unf_10$pval.exposure)
crp_il6_unf_10 <- subset(crp_il6_unf_10, pval.exposure < 5e-8)

# =============================================================================
# Perform LD clumping to select independent variants
# =============================================================================
# Clump variants (r² < 0.1) to ensure independence
crp_il6_unf_c <- clump_data(crp_il6_unf, clump_r2 = 0.1, clump_p1 = 1, clump_p2 = 1)
crp_il6_unf_100_c <- clump_data(crp_il6_unf_100, clump_r2 = 0.1, clump_p1 = 1, clump_p2 = 1)
crp_il6_unf_10_c <- clump_data(crp_il6_unf_10, clump_r2 = 0.1, clump_p1 = 1, clump_p2 = 1)

# Clean up intermediate objects
rm(crp, crp_il6_unf)

# =============================================================================
# Save genetic instruments
# =============================================================================
# Save the three genetic instruments
write.table(crp_il6_unf_c, "data/instruments/il6_instruments_300.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(crp_il6_unf_100_c, "data/instruments/il6_instruments_100.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(crp_il6_unf_10_c, "data/instruments/il6_instruments_10.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")