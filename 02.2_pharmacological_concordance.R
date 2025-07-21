# =============================================================================
# Pharmacological Concordance Analysis
# =============================================================================
# Purpose: Validate IL6 genetic instruments against clinical trial biomarkers
#          and autoimmune outcomes to demonstrate concordance with IL-6 inhibition

# Load required libraries
library(data.table)
library(MendelianRandomization)
library(plyr)
library(TwoSampleMR)
library(readxl)
library("survival")
library(dplyr)
library(SUMnlmr)
library(metafor) 
library(ggplot2)
library(meta)

# =============================================================================
# Load genetic instruments
# =============================================================================
# Read the three IL6 genetic instruments constructed in step 01
crp_il6_unf_c <- fread("data/instruments/il6_instruments_300.txt")
crp_il6_unf_100_c <- fread("data/instruments/il6_instruments_100.txt")
crp_il6_unf_10_c <- fread("data/instruments/il6_instruments_10.txt")

# =============================================================================
# Biomarker Validation: Lipoprotein(a)
# =============================================================================
# Load lipoprotein(a) GWAS data
lpa <- fread("data/raw/Lipoprotein_A.imp.gz")
colnames(lpa)[1] <- "chr_name"
colnames(lpa)[2] <- "chr_start"

# Merge with genetic instruments
lpa_il6 <- merge(lpa, crp_il6_unf_c, by = c("chr_name", "chr_start"))
lpa_il6_100 <- merge(lpa, crp_il6_unf_100_c, by = c("chr_name", "chr_start"))
lpa_il6_10 <- merge(lpa, crp_il6_unf_10_c, by = c("chr_name", "chr_start"))
rm(lpa)

# Prepare data for MR analysis - 300kb instrument
lpa_il6$gx <- lpa_il6$beta
lpa_il6$gx_se <- lpa_il6$se
lpa_il6$gy <- ifelse(tolower(lpa_il6$effect_allele) == tolower(lpa_il6$ALT), 
                     lpa_il6$Effect, -lpa_il6$Effect)
lpa_il6$gy_se <- lpa_il6$StdErr

# Run MR analysis
mr_allmethods(mr_input(bx = lpa_il6$gx, bxse = lpa_il6$gx_se, 
                       by = lpa_il6$gy, 
                       byse = lpa_il6$gy_se), method = "all")
mr_ivw_fe(lpa_il6$gx, lpa_il6$gy, lpa_il6$gx_se, lpa_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
lpa_il6_100$gx <- lpa_il6_100$beta
lpa_il6_100$gx_se <- lpa_il6_100$se
lpa_il6_100$gy <- ifelse(tolower(lpa_il6_100$effect_allele) == tolower(lpa_il6_100$ALT), 
                         lpa_il6_100$Effect, -lpa_il6_100$Effect)
lpa_il6_100$gy_se <- lpa_il6_100$StdErr

mr_allmethods(mr_input(bx = lpa_il6_100$gx, bxse = lpa_il6_100$gx_se, 
                       by = lpa_il6_100$gy, 
                       byse = lpa_il6_100$gy_se), method = "all")
mr_ivw_fe(lpa_il6_100$gx, lpa_il6_100$gy, lpa_il6_100$gx_se, lpa_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
lpa_il6_10$gx <- lpa_il6_10$beta
lpa_il6_10$gx_se <- lpa_il6_10$se
lpa_il6_10$gy <- ifelse(tolower(lpa_il6_10$effect_allele) == tolower(lpa_il6_10$ALT), 
                        lpa_il6_10$Effect, -lpa_il6_10$Effect)
lpa_il6_10$gy_se <- lpa_il6_10$StdErr

mr_allmethods(mr_input(bx = lpa_il6_10$gx, bxse = lpa_il6_10$gx_se, 
                       by = lpa_il6_10$gy, 
                       byse = lpa_il6_10$gy_se), method = "all")
mr_ivw_fe(lpa_il6_10$gx, lpa_il6_10$gy, lpa_il6_10$gx_se, lpa_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(lpa_il6, "results/lpa_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(lpa_il6_100, "results/lpa_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(lpa_il6_10, "results/lpa_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Biomarker Validation: Fibrinogen
# =============================================================================
# Load fibrinogen GWAS data
fibr <- fread("data/raw/ln_fibrinogen_CHARGEconsortium_2015_filtered.txt")
colnames(fibr)[2] <- "SNP"

# Merge with genetic instruments
fibr_il6 <- merge(fibr, crp_il6_unf_c, by = "SNP")
fibr_il6_100 <- merge(fibr, crp_il6_unf_100_c, by = "SNP")
fibr_il6_10 <- merge(fibr, crp_il6_unf_10_c, by = "SNP")
rm(fibr)

# Prepare data for MR analysis - 300kb instrument
fibr_il6$gx <- fibr_il6$beta
fibr_il6$gx_se <- fibr_il6$se
fibr_il6$gy <- ifelse(tolower(fibr_il6$Allele1) == tolower(fibr_il6$effect_allele), 
                      fibr_il6$Effect, -fibr_il6$Effect)
fibr_il6$gy_se <- fibr_il6$StdErr

# Run MR analysis
mr_allmethods(mr_input(bx = fibr_il6$gx, bxse = fibr_il6$gx_se, 
                       by = fibr_il6$gy, 
                       byse = fibr_il6$gy_se), method = "all")
mr_ivw_fe(fibr_il6$gx, fibr_il6$gy, fibr_il6$gx_se, fibr_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
fibr_il6_100$gx <- fibr_il6_100$beta
fibr_il6_100$gx_se <- fibr_il6_100$se
fibr_il6_100$gy <- ifelse(tolower(fibr_il6_100$Allele1) == tolower(fibr_il6_100$effect_allele), 
                          fibr_il6_100$Effect, -fibr_il6_100$Effect)
fibr_il6_100$gy_se <- fibr_il6_100$StdErr

mr_allmethods(mr_input(bx = fibr_il6_100$gx, bxse = fibr_il6_100$gx_se, 
                       by = fibr_il6_100$gy, 
                       byse = fibr_il6_100$gy_se), method = "all")
mr_ivw_fe(fibr_il6_100$gx, fibr_il6_100$gy, fibr_il6_100$gx_se, fibr_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
fibr_il6_10$gx <- fibr_il6_10$beta
fibr_il6_10$gx_se <- fibr_il6_10$se
fibr_il6_10$gy <- ifelse(tolower(fibr_il6_10$Allele1) == tolower(fibr_il6_10$effect_allele), 
                         fibr_il6_10$Effect, -fibr_il6_10$Effect)
fibr_il6_10$gy_se <- fibr_il6_10$StdErr

mr_allmethods(mr_input(bx = fibr_il6_10$gx, bxse = fibr_il6_10$gx_se, 
                       by = fibr_il6_10$gy, 
                       byse = fibr_il6_10$gy_se), method = "all")
mr_ivw_fe(fibr_il6_10$gx, fibr_il6_10$gy, fibr_il6_10$gx_se, fibr_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(fibr_il6, "results/fibr_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(fibr_il6_100, "results/fibr_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(fibr_il6_10, "results/fibr_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Biomarker Validation: Serum Amyloid A (SAA)
# =============================================================================
# Load SAA GWAS data from Pietzner 2020, Nat Commun
saa <- fread("data/raw/saa_GCST90019394_buildGRCh37.tsv")
colnames(saa)[3] <- "SNP"

# Merge with genetic instruments
saa_il6 <- merge(saa, crp_il6_unf_c, by = "SNP")
saa_il6_100 <- merge(saa, crp_il6_unf_100_c, by = "SNP")
saa_il6_10 <- merge(saa, crp_il6_unf_10_c, by = "SNP")
rm(saa)

# Prepare data for MR analysis - 300kb instrument
saa_il6$gx <- saa_il6$beta.y
saa_il6$gx_se <- saa_il6$se
saa_il6$gy <- ifelse(tolower(saa_il6$effect_allele.x) == tolower(saa_il6$effect_allele.y), 
                     saa_il6$beta.x, -saa_il6$beta.x)
saa_il6$gy_se <- saa_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = saa_il6$gx, bxse = saa_il6$gx_se, 
                       by = saa_il6$gy, 
                       byse = saa_il6$gy_se), method = "all")
mr_ivw_fe(saa_il6$gx, saa_il6$gy, saa_il6$gx_se, saa_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
saa_il6_100$gx <- saa_il6_100$beta.y
saa_il6_100$gx_se <- saa_il6_100$se
saa_il6_100$gy <- ifelse(tolower(saa_il6_100$effect_allele.x) == tolower(saa_il6_100$effect_allele.y), 
                         saa_il6_100$beta.x, -saa_il6_100$beta.x)
saa_il6_100$gy_se <- saa_il6_100$standard_error

mr_allmethods(mr_input(bx = saa_il6_100$gx, bxse = saa_il6_100$gx_se, 
                       by = saa_il6_100$gy, 
                       byse = saa_il6_100$gy_se), method = "all")
mr_ivw_fe(saa_il6_100$gx, saa_il6_100$gy, saa_il6_100$gx_se, saa_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
saa_il6_10$gx <- saa_il6_10$beta.y
saa_il6_10$gx_se <- saa_il6_10$se
saa_il6_10$gy <- ifelse(tolower(saa_il6_10$effect_allele.x) == tolower(saa_il6_10$effect_allele.y), 
                        saa_il6_10$beta.x, -saa_il6_10$beta.x)
saa_il6_10$gy_se <- saa_il6_10$standard_error

mr_allmethods(mr_input(bx = saa_il6_10$gx, bxse = saa_il6_10$gx_se, 
                       by = saa_il6_10$gy, 
                       byse = saa_il6_10$gy_se), method = "all")
mr_ivw_fe(saa_il6_10$gx, saa_il6_10$gy, saa_il6_10$gx_se, saa_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(saa_il6, "results/saa_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(saa_il6_100, "results/saa_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(saa_il6_10, "results/saa_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Biomarker Validation: Haptoglobin
# =============================================================================
# Load haptoglobin GWAS data from Pietzner 2020, Nat Commun
hapto <- fread("data/raw/haptoglobin_GCST90019470_buildGRCh37.tsv")
colnames(hapto)[3] <- "SNP"

# Merge with genetic instruments
hapto_il6 <- merge(hapto, crp_il6_unf_c, by = "SNP")
hapto_il6_100 <- merge(hapto, crp_il6_unf_100_c, by = "SNP")
hapto_il6_10 <- merge(hapto, crp_il6_unf_10_c, by = "SNP")
rm(hapto)

# Prepare data for MR analysis - 300kb instrument
hapto_il6$gx <- hapto_il6$beta.y
hapto_il6$gx_se <- hapto_il6$se
hapto_il6$gy <- ifelse(tolower(hapto_il6$effect_allele.x) == tolower(hapto_il6$effect_allele.y), 
                       hapto_il6$beta.x, -hapto_il6$beta.x)
hapto_il6$gy_se <- hapto_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = hapto_il6$gx, bxse = hapto_il6$gx_se, 
                       by = hapto_il6$gy, 
                       byse = hapto_il6$gy_se), method = "all")
mr_ivw_fe(hapto_il6$gx, hapto_il6$gy, hapto_il6$gx_se, hapto_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
hapto_il6_100$gx <- hapto_il6_100$beta.y
hapto_il6_100$gx_se <- hapto_il6_100$se
hapto_il6_100$gy <- ifelse(tolower(hapto_il6_100$effect_allele.x) == tolower(hapto_il6_100$effect_allele.y), 
                           hapto_il6_100$beta.x, -hapto_il6_100$beta.x)
hapto_il6_100$gy_se <- hapto_il6_100$standard_error

mr_allmethods(mr_input(bx = hapto_il6_100$gx, bxse = hapto_il6_100$gx_se, 
                       by = hapto_il6_100$gy, 
                       byse = hapto_il6_100$gy_se), method = "all")
mr_ivw_fe(hapto_il6_100$gx, hapto_il6_100$gy, hapto_il6_100$gx_se, hapto_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
hapto_il6_10$gx <- hapto_il6_10$beta.y
hapto_il6_10$gx_se <- hapto_il6_10$se
hapto_il6_10$gy <- ifelse(tolower(hapto_il6_10$effect_allele.x) == tolower(hapto_il6_10$effect_allele.y), 
                          hapto_il6_10$beta.x, -hapto_il6_10$beta.x)
hapto_il6_10$gy_se <- hapto_il6_10$standard_error

mr_allmethods(mr_input(bx = hapto_il6_10$gx, bxse = hapto_il6_10$gx_se, 
                       by = hapto_il6_10$gy, 
                       byse = hapto_il6_10$gy_se), method = "all")
mr_ivw_fe(hapto_il6_10$gx, hapto_il6_10$gy, hapto_il6_10$gx_se, hapto_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(hapto_il6, "results/hapto_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(hapto_il6_100, "results/hapto_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(hapto_il6_10, "results/hapto_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# IL-6 Protein Level Analysis
# =============================================================================
# Load IL-6 protein GWAS data from Konieczny et al 2025 Communications Biology
# Same code for IL6 in CSF GWAS data from Western et al 2024 Nat Genet.
il6 <- fread("data/raw/IL_61.TBL")
colnames(il6)[1] <- "SNP"

# Merge with genetic instruments
il6_il6 <- merge(il6, crp_il6_unf_c, by = "SNP")
il6_il6_100 <- merge(il6, crp_il6_unf_100_c, by = "SNP")
il6_il6_10 <- merge(il6, crp_il6_unf_10_c, by = "SNP")
rm(il6)

# Prepare data for MR analysis - 300kb instrument
il6_il6$gx <- il6_il6$beta
il6_il6$gx_se <- il6_il6$se
# Calculate standard error for outcome using Z-score method
il6_il6$gy_se <- 1 / sqrt((2 * il6_il6$maf) * (1 - il6_il6$maf) * (il6_il6$Weight + il6_il6$Zscore^2))
il6_il6$Effect <- il6_il6$gy_se * il6_il6$Zscore
il6_il6$gy <- ifelse(tolower(il6_il6$Allele1) == tolower(il6_il6$effect_allele), 
                     il6_il6$Effect, -il6_il6$Effect)

# Run MR analysis
mr_allmethods(mr_input(bx = il6_il6$gx, bxse = il6_il6$gx_se, 
                       by = il6_il6$gy, 
                       byse = il6_il6$gy_se), method = "all")
mr_ivw_fe(il6_il6$gx, il6_il6$gy, il6_il6$gx_se, il6_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
il6_il6_100$gx <- il6_il6_100$beta
il6_il6_100$gx_se <- il6_il6_100$se
il6_il6_100$gy_se <- 1 / sqrt((2 * il6_il6_100$maf) * (1 - il6_il6_100$maf) * (il6_il6_100$Weight + il6_il6_100$Zscore^2))
il6_il6_100$Effect <- il6_il6_100$gy_se * il6_il6_100$Zscore
il6_il6_100$gy <- ifelse(tolower(il6_il6_100$Allele1) == tolower(il6_il6_100$effect_allele), 
                         il6_il6_100$Effect, -il6_il6_100$Effect)

mr_allmethods(mr_input(bx = il6_il6_100$gx, bxse = il6_il6_100$gx_se, 
                       by = il6_il6_100$gy, 
                       byse = il6_il6_100$gy_se), method = "all")
mr_ivw_fe(il6_il6_100$gx, il6_il6_100$gy, il6_il6_100$gx_se, il6_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
il6_il6_10$gx <- il6_il6_10$beta
il6_il6_10$gx_se <- il6_il6_10$se
il6_il6_10$gy_se <- 1 / sqrt((2 * il6_il6_10$maf) * (1 - il6_il6_10$maf) * (il6_il6_10$Weight + il6_il6_10$Zscore^2))
il6_il6_10$Effect <- il6_il6_10$gy_se * il6_il6_10$Zscore
il6_il6_10$gy <- ifelse(tolower(il6_il6_10$Allele1) == tolower(il6_il6_10$effect_allele), 
                        il6_il6_10$Effect, -il6_il6_10$Effect)

mr_allmethods(mr_input(bx = il6_il6_10$gx, bxse = il6_il6_10$gx_se, 
                       by = il6_il6_10$gy, 
                       byse = il6_il6_10$gy_se), method = "all")
mr_ivw_fe(il6_il6_10$gx, il6_il6_10$gy, il6_il6_10$gx_se, il6_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(il6_il6, "results/il6_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(il6_il6_100, "results/il6_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(il6_il6_10, "results/il6_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Autoimmune Disease Validation: Polymyalgia Rheumatica
# =============================================================================
# Load polymyalgia rheumatica GWAS data
polym_rheum <- fread("data/raw/Meta_PMR.txt")
colnames(polym_rheum)[1] <- "SNP"

# Merge with genetic instruments
polym_rheum_il6 <- merge(polym_rheum, crp_il6_unf_c, by = "SNP")
polym_rheum_il6_100 <- merge(polym_rheum, crp_il6_unf_100_c, by = "SNP")
polym_rheum_il6_10 <- merge(polym_rheum, crp_il6_unf_10_c, by = "SNP")
rm(polym_rheum)

# Prepare data for MR analysis - 300kb instrument
polym_rheum_il6$gx <- polym_rheum_il6$beta
polym_rheum_il6$gx_se <- polym_rheum_il6$se
polym_rheum_il6$gy <- ifelse(tolower(polym_rheum_il6$EA) == tolower(polym_rheum_il6$effect_allele), 
                             polym_rheum_il6$Beta, -polym_rheum_il6$Beta)
polym_rheum_il6$gy_se <- polym_rheum_il6$StdErr

# Run MR analysis
mr_allmethods(mr_input(bx = polym_rheum_il6$gx, bxse = polym_rheum_il6$gx_se, 
                       by = polym_rheum_il6$gy, 
                       byse = polym_rheum_il6$gy_se), method = "all")
mr_ivw_fe(polym_rheum_il6$gx, polym_rheum_il6$gy, polym_rheum_il6$gx_se, polym_rheum_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
polym_rheum_il6_100$gx <- polym_rheum_il6_100$beta
polym_rheum_il6_100$gx_se <- polym_rheum_il6_100$se
polym_rheum_il6_100$gy <- ifelse(tolower(polym_rheum_il6_100$EA) == tolower(polym_rheum_il6_100$effect_allele), 
                                 polym_rheum_il6_100$Beta, -polym_rheum_il6_100$Beta)
polym_rheum_il6_100$gy_se <- polym_rheum_il6_100$StdErr

mr_allmethods(mr_input(bx = polym_rheum_il6_100$gx, bxse = polym_rheum_il6_100$gx_se, 
                       by = polym_rheum_il6_100$gy, 
                       byse = polym_rheum_il6_100$gy_se), method = "all")
mr_ivw_fe(polym_rheum_il6_100$gx, polym_rheum_il6_100$gy, polym_rheum_il6_100$gx_se, polym_rheum_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
polym_rheum_il6_10$gx <- polym_rheum_il6_10$beta
polym_rheum_il6_10$gx_se <- polym_rheum_il6_10$se
polym_rheum_il6_10$gy <- ifelse(tolower(polym_rheum_il6_10$EA) == tolower(polym_rheum_il6_10$effect_allele), 
                                polym_rheum_il6_10$Beta, -polym_rheum_il6_10$Beta)
polym_rheum_il6_10$gy_se <- polym_rheum_il6_10$StdErr

mr_allmethods(mr_input(bx = polym_rheum_il6_10$gx, bxse = polym_rheum_il6_10$gx_se, 
                       by = polym_rheum_il6_10$gy, 
                       byse = polym_rheum_il6_10$gy_se), method = "all")
mr_ivw_fe(polym_rheum_il6_10$gx, polym_rheum_il6_10$gy, polym_rheum_il6_10$gx_se, polym_rheum_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(polym_rheum_il6, "results/polym_rheum_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(polym_rheum_il6_100, "results/polym_rheum_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(polym_rheum_il6_10, "results/polym_rheum_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Autoimmune Disease Validation: Rheumatoid Arthritis
# =============================================================================
# Load rheumatoid arthritis GWAS data
rheum <- fread("data/raw/GCST90132222_buildGRCh37.tsv")
colnames(rheum)[1] <- "SNP"

# Merge with genetic instruments
rheum_il6 <- merge(rheum, crp_il6_unf_c, by = "SNP")
rheum_il6_100 <- merge(rheum, crp_il6_unf_100_c, by = "SNP")
rheum_il6_10 <- merge(rheum, crp_il6_unf_10_c, by = "SNP")
rm(rheum)

# Prepare data for MR analysis - 300kb instrument
rheum_il6$gx <- rheum_il6$beta.y
rheum_il6$gx_se <- rheum_il6$se
rheum_il6$gy <- ifelse(tolower(rheum_il6$effect_allele.x) == tolower(rheum_il6$effect_allele.y), 
                       rheum_il6$beta.x, -rheum_il6$beta.x)
rheum_il6$gy_se <- rheum_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = rheum_il6$gx, bxse = rheum_il6$gx_se, 
                       by = rheum_il6$gy, 
                       byse = rheum_il6$gy_se), method = "all")
mr_ivw_fe(rheum_il6$gx, rheum_il6$gy, rheum_il6$gx_se, rheum_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
rheum_il6_100$gx <- rheum_il6_100$beta.y
rheum_il6_100$gx_se <- rheum_il6_100$se
rheum_il6_100$gy <- ifelse(tolower(rheum_il6_100$effect_allele.x) == tolower(rheum_il6_100$effect_allele.y), 
                           rheum_il6_100$beta.x, -rheum_il6_100$beta.x)
rheum_il6_100$gy_se <- rheum_il6_100$standard_error

mr_allmethods(mr_input(bx = rheum_il6_100$gx, bxse = rheum_il6_100$gx_se, 
                       by = rheum_il6_100$gy, 
                       byse = rheum_il6_100$gy_se), method = "all")
mr_ivw_fe(rheum_il6_100$gx, rheum_il6_100$gy, rheum_il6_100$gx_se, rheum_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
rheum_il6_10$gx <- rheum_il6_10$beta.y
rheum_il6_10$gx_se <- rheum_il6_10$se
rheum_il6_10$gy <- ifelse(tolower(rheum_il6_10$effect_allele.x) == tolower(rheum_il6_10$effect_allele.y), 
                          rheum_il6_10$beta.x, -rheum_il6_10$beta.x)
rheum_il6_10$gy_se <- rheum_il6_10$standard_error

mr_allmethods(mr_input(bx = rheum_il6_10$gx, bxse = rheum_il6_10$gx_se, 
                       by = rheum_il6_10$gy, 
                       byse = rheum_il6_10$gy_se), method = "all")
mr_ivw_fe(rheum_il6_10$gx, rheum_il6_10$gy, rheum_il6_10$gx_se, rheum_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(rheum_il6, "results/rheum_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(rheum_il6_100, "results/rheum_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(rheum_il6_10, "results/rheum_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")