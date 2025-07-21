# =============================================================================
# Cardiovascular Endpoints Analysis
# =============================================================================
# Purpose: Test associations between IL6 genetic instruments and cardiovascular
#          diseases including CAD, stroke, PAD, and subclinical atherosclerosis

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
# Coronary Artery Disease (CAD)
# =============================================================================
# Load CAD GWAS data
cad <- fread("data/raw/GCST90132315.h.tsv.gz")
colnames(cad)[25] <- "SNP"

# Merge with genetic instruments
cad_il6 <- merge(cad, crp_il6_unf_c, by = "SNP")
cad_il6_100 <- merge(cad, crp_il6_unf_100_c, by = "SNP")
cad_il6_10 <- merge(cad, crp_il6_unf_10_c, by = "SNP")
rm(cad)

# Prepare data for MR analysis - 300kb instrument
cad_il6$gx <- cad_il6$beta.y
cad_il6$gx_se <- cad_il6$se
cad_il6$gy <- ifelse(tolower(cad_il6$effect_allele.x) == tolower(cad_il6$effect_allele.y), 
                     cad_il6$beta.x, -cad_il6$beta.x)
cad_il6$gy_se <- cad_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = cad_il6$gx, bxse = cad_il6$gx_se, 
                       by = cad_il6$gy, 
                       byse = cad_il6$gy_se), method = "all")
mr_ivw_fe(cad_il6$gx, cad_il6$gy, cad_il6$gx_se, cad_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
cad_il6_100$gx <- cad_il6_100$beta.y
cad_il6_100$gx_se <- cad_il6_100$se
cad_il6_100$gy <- ifelse(tolower(cad_il6_100$effect_allele.x) == tolower(cad_il6_100$effect_allele.y), 
                         cad_il6_100$beta.x, -cad_il6_100$beta.x)
cad_il6_100$gy_se <- cad_il6_100$standard_error

mr_allmethods(mr_input(bx = cad_il6_100$gx, bxse = cad_il6_100$gx_se, 
                       by = cad_il6_100$gy, 
                       byse = cad_il6_100$gy_se), method = "all")
mr_ivw_fe(cad_il6_100$gx, cad_il6_100$gy, cad_il6_100$gx_se, cad_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
cad_il6_10$gx <- cad_il6_10$beta.y
cad_il6_10$gx_se <- cad_il6_10$se
cad_il6_10$gy <- ifelse(tolower(cad_il6_10$effect_allele.x) == tolower(cad_il6_10$effect_allele.y), 
                        cad_il6_10$beta.x, -cad_il6_10$beta.x)
cad_il6_10$gy_se <- cad_il6_10$standard_error

mr_allmethods(mr_input(bx = cad_il6_10$gx, bxse = cad_il6_10$gx_se, 
                       by = cad_il6_10$gy, 
                       byse = cad_il6_10$gy_se), method = "all")
mr_ivw_fe(cad_il6_10$gx, cad_il6_10$gy, cad_il6_10$gx_se, cad_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(cad_il6, "results/cad_il6_300.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cad_il6_100, "results/cad_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cad_il6_10, "results/cad_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Myocardial Infarction (MI)
# =============================================================================
# Load MI GWAS data from CARDIoGRAMplusC4D
mi <- fread("data/raw/mi.add.030315.website.txt")
colnames(mi)[1] <- "SNP"

# Merge with genetic instruments
mi_il6 <- merge(mi, crp_il6_unf_c, by = "SNP")
mi_il6_100 <- merge(mi, crp_il6_unf_100_c, by = "SNP")
mi_il6_10 <- merge(mi, crp_il6_unf_10_c, by = "SNP")
rm(mi)

# Prepare data for MR analysis - 300kb instrument
mi_il6$gx <- mi_il6$beta.y
mi_il6$gx_se <- mi_il6$se
mi_il6$gy <- ifelse(tolower(mi_il6$effect_allele.x) == tolower(mi_il6$effect_allele.y), 
                    mi_il6$beta.x, -mi_il6$beta.x)
mi_il6$gy_se <- mi_il6$se_dgc

# Run MR analysis
mr_allmethods(mr_input(bx = mi_il6$gx, bxse = mi_il6$gx_se, 
                       by = mi_il6$gy, 
                       byse = mi_il6$gy_se), method = "all")
mr_ivw_fe(mi_il6$gx, mi_il6$gy, mi_il6$gx_se, mi_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
mi_il6_100$gx <- mi_il6_100$beta.y
mi_il6_100$gx_se <- mi_il6_100$se
mi_il6_100$gy <- ifelse(tolower(mi_il6_100$effect_allele.x) == tolower(mi_il6_100$effect_allele.y), 
                        mi_il6_100$beta.x, -mi_il6_100$beta.x)
mi_il6_100$gy_se <- mi_il6_100$se_dgc

mr_allmethods(mr_input(bx = mi_il6_100$gx, bxse = mi_il6_100$gx_se, 
                       by = mi_il6_100$gy, 
                       byse = mi_il6_100$gy_se), method = "all")
mr_ivw_fe(mi_il6_100$gx, mi_il6_100$gy, mi_il6_100$gx_se, mi_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
mi_il6_10$gx <- mi_il6_10$beta.y
mi_il6_10$gx_se <- mi_il6_10$se
mi_il6_10$gy <- ifelse(tolower(mi_il6_10$effect_allele.x) == tolower(mi_il6_10$effect_allele.y), 
                       mi_il6_10$beta.x, -mi_il6_10$beta.x)
mi_il6_10$gy_se <- mi_il6_10$se_dgc

mr_allmethods(mr_input(bx = mi_il6_10$gx, bxse = mi_il6_10$gx_se, 
                       by = mi_il6_10$gy, 
                       byse = mi_il6_10$gy_se), method = "all")
mr_ivw_fe(mi_il6_10$gx, mi_il6_10$gy, mi_il6_10$gx_se, mi_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(mi_il6, "results/mi_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mi_il6_100, "results/mi_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mi_il6_10, "results/mi_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Peripheral Artery Disease (PAD)
# =============================================================================
# Load PAD GWAS data from MVP
pad <- fread("data/raw/MVP.PAD.Transethnic_split.tsv.gz")
colnames(pad)[1] <- "chr_name"
colnames(pad)[2] <- "chr_start"

# Merge with genetic instruments
pad_il6 <- merge(pad, crp_il6_unf_c, by = c("chr_name", "chr_start"))
pad_il6_100 <- merge(pad, crp_il6_unf_100_c, by = c("chr_name", "chr_start"))
pad_il6_10 <- merge(pad, crp_il6_unf_10_c, by = c("chr_name", "chr_start"))
rm(pad)

# Prepare data for MR analysis - 300kb instrument
pad_il6$gx <- pad_il6$beta
pad_il6$gx_se <- pad_il6$se
pad_il6$gy <- ifelse(tolower(pad_il6$Allele1) == tolower(pad_il6$effect_allele), 
                     pad_il6$Effect, -pad_il6$Effect)
pad_il6$gy_se <- pad_il6$StdErr

# Run MR analysis
mr_allmethods(mr_input(bx = pad_il6$gx, bxse = pad_il6$gx_se, 
                       by = pad_il6$gy, 
                       byse = pad_il6$gy_se), method = "all")
mr_ivw_fe(pad_il6$gx, pad_il6$gy, pad_il6$gx_se, pad_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
pad_il6_100$gx <- pad_il6_100$beta
pad_il6_100$gx_se <- pad_il6_100$se
pad_il6_100$gy <- ifelse(tolower(pad_il6_100$Allele1) == tolower(pad_il6_100$effect_allele), 
                         pad_il6_100$Effect, -pad_il6_100$Effect)
pad_il6_100$gy_se <- pad_il6_100$StdErr

mr_allmethods(mr_input(bx = pad_il6_100$gx, bxse = pad_il6_100$gx_se, 
                       by = pad_il6_100$gy, 
                       byse = pad_il6_100$gy_se), method = "all")
mr_ivw_fe(pad_il6_100$gx, pad_il6_100$gy, pad_il6_100$gx_se, pad_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
pad_il6_10$gx <- pad_il6_10$beta
pad_il6_10$gx_se <- pad_il6_10$se
pad_il6_10$gy <- ifelse(tolower(pad_il6_10$Allele1) == tolower(pad_il6_10$effect_allele), 
                        pad_il6_10$Effect, -pad_il6_10$Effect)
pad_il6_10$gy_se <- pad_il6_10$StdErr

mr_allmethods(mr_input(bx = pad_il6_10$gx, bxse = pad_il6_10$gx_se, 
                       by = pad_il6_10$gy, 
                       byse = pad_il6_10$gy_se), method = "all")
mr_ivw_fe(pad_il6_10$gx, pad_il6_10$gy, pad_il6_10$gx_se, pad_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(pad_il6, "results/pad_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(pad_il6_100, "results/pad_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(pad_il6_10, "results/pad_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Any Ischemic Stroke
# =============================================================================
# Load any ischemic stroke GWAS data from GIGASTROKE
anyis <- fread("data/raw/GCST90104535_buildGRCh37_ANYIS.tsv.gz")
colnames(anyis)[1] <- "chr_name"
colnames(anyis)[2] <- "chr_start"

# Merge with genetic instruments
anyis_il6 <- merge(anyis, crp_il6_unf_c, by = c("chr_name", "chr_start"))
anyis_il6_100 <- merge(anyis, crp_il6_unf_100_c, by = c("chr_name", "chr_start"))
anyis_il6_10 <- merge(anyis, crp_il6_unf_10_c, by = c("chr_name", "chr_start"))
rm(anyis)

# Prepare data for MR analysis - 300kb instrument
anyis_il6$gx <- anyis_il6$beta.y
anyis_il6$gx_se <- anyis_il6$se 
anyis_il6$gy <- ifelse(tolower(anyis_il6$effect_allele.x) == tolower(anyis_il6$effect_allele.y), 
                       anyis_il6$beta.x, -anyis_il6$beta.x)
anyis_il6$gy_se <- anyis_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = anyis_il6$gx, bxse = anyis_il6$gx_se, 
                       by = anyis_il6$gy, 
                       byse = anyis_il6$gy_se), method = "all")
mr_ivw_fe(anyis_il6$gx, anyis_il6$gy, anyis_il6$gx_se, anyis_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
anyis_il6_100$gx <- anyis_il6_100$beta.y
anyis_il6_100$gx_se <- anyis_il6_100$se 
anyis_il6_100$gy <- ifelse(tolower(anyis_il6_100$effect_allele.x) == tolower(anyis_il6_100$effect_allele.y), 
                           anyis_il6_100$beta.x, -anyis_il6_100$beta.x)
anyis_il6_100$gy_se <- anyis_il6_100$standard_error

mr_allmethods(mr_input(bx = anyis_il6_100$gx, bxse = anyis_il6_100$gx_se, 
                       by = anyis_il6_100$gy, 
                       byse = anyis_il6_100$gy_se), method = "all")
mr_ivw_fe(anyis_il6_100$gx, anyis_il6_100$gy, anyis_il6_100$gx_se, anyis_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
anyis_il6_10$gx <- anyis_il6_10$beta.y
anyis_il6_10$gx_se <- anyis_il6_10$se 
anyis_il6_10$gy <- ifelse(tolower(anyis_il6_10$effect_allele.x) == tolower(anyis_il6_10$effect_allele.y), 
                          anyis_il6_10$beta.x, -anyis_il6_10$beta.x)
anyis_il6_10$gy_se <- anyis_il6_10$standard_error

mr_allmethods(mr_input(bx = anyis_il6_10$gx, bxse = anyis_il6_10$gx_se, 
                       by = anyis_il6_10$gy, 
                       byse = anyis_il6_10$gy_se), method = "all")
mr_ivw_fe(anyis_il6_10$gx, anyis_il6_10$gy, anyis_il6_10$gx_se, anyis_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(anyis_il6, "results/anyis_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(anyis_il6_100, "results/anyis_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(anyis_il6_10, "results/anyis_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Large Artery Atherosclerotic Stroke
# =============================================================================
# Load large artery stroke GWAS data from GIGASTROKE
las <- fread("data/raw/GCST90104538_buildGRCh37_LAA.tsv.gz")
colnames(las)[1] <- "chr_name"
colnames(las)[2] <- "chr_start"

# Merge with genetic instruments
las_il6 <- merge(las, crp_il6_unf_c, by = c("chr_name", "chr_start"))
las_il6_100 <- merge(las, crp_il6_unf_100_c, by = c("chr_name", "chr_start"))
las_il6_10 <- merge(las, crp_il6_unf_10_c, by = c("chr_name", "chr_start"))
rm(las)

# Prepare data for MR analysis - 300kb instrument
las_il6$gx <- las_il6$beta.y
las_il6$gx_se <- las_il6$se 
las_il6$gy <- ifelse(tolower(las_il6$effect_allele.x) == tolower(las_il6$effect_allele.y), 
                     las_il6$beta.x, -las_il6$beta.x)
las_il6$gy_se <- las_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = las_il6$gx, bxse = las_il6$gx_se, 
                       by = las_il6$gy, 
                       byse = las_il6$gy_se), method = "all")
mr_ivw_fe(las_il6$gx, las_il6$gy, las_il6$gx_se, las_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
las_il6_100$gx <- las_il6_100$beta.y
las_il6_100$gx_se <- las_il6_100$se 
las_il6_100$gy <- ifelse(tolower(las_il6_100$effect_allele.x) == tolower(las_il6_100$effect_allele.y), 
                         las_il6_100$beta.x, -las_il6_100$beta.x)
las_il6_100$gy_se <- las_il6_100$standard_error

mr_allmethods(mr_input(bx = las_il6_100$gx, bxse = las_il6_100$gx_se, 
                       by = las_il6_100$gy, 
                       byse = las_il6_100$gy_se), method = "all")
mr_ivw_fe(las_il6_100$gx, las_il6_100$gy, las_il6_100$gx_se, las_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
las_il6_10$gx <- las_il6_10$beta.y
las_il6_10$gx_se <- las_il6_10$se 
las_il6_10$gy <- ifelse(tolower(las_il6_10$effect_allele.x) == tolower(las_il6_10$effect_allele.y), 
                        las_il6_10$beta.x, -las_il6_10$beta.x)
las_il6_10$gy_se <- las_il6_10$standard_error

mr_allmethods(mr_input(bx = las_il6_10$gx, bxse = las_il6_10$gx_se, 
                       by = las_il6_10$gy, 
                       byse = las_il6_10$gy_se), method = "all")
mr_ivw_fe(las_il6_10$gx, las_il6_10$gy, las_il6_10$gx_se, las_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(las_il6, "results/las_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(las_il6_100, "results/las_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(las_il6_10, "results/las_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Cardioembolic Stroke
# =============================================================================
# Load cardioembolic stroke GWAS data from GIGASTROKE
ces <- fread("data/raw/GCST90104536_buildGRCh37_CE.tsv.gz")
colnames(ces)[1] <- "chr_name"
colnames(ces)[2] <- "chr_start"

# Merge with genetic instruments
ces_il6 <- merge(ces, crp_il6_unf_c, by = c("chr_name", "chr_start"))
ces_il6_100 <- merge(ces, crp_il6_unf_100_c, by = c("chr_name", "chr_start"))
ces_il6_10 <- merge(ces, crp_il6_unf_10_c, by = c("chr_name", "chr_start"))
rm(ces)

# Prepare data for MR analysis - 300kb instrument
ces_il6$gx <- ces_il6$beta.y
ces_il6$gx_se <- ces_il6$se 
ces_il6$gy <- ifelse(tolower(ces_il6$effect_allele.x) == tolower(ces_il6$effect_allele.y), 
                     ces_il6$beta.x, -ces_il6$beta.x)
ces_il6$gy_se <- ces_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = ces_il6$gx, bxse = ces_il6$gx_se, 
                       by = ces_il6$gy, 
                       byse = ces_il6$gy_se), method = "all")
mr_ivw_fe(ces_il6$gx, ces_il6$gy, ces_il6$gx_se, ces_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
ces_il6_100$gx <- ces_il6_100$beta.y
ces_il6_100$gx_se <- ces_il6_100$se 
ces_il6_100$gy <- ifelse(tolower(ces_il6_100$effect_allele.x) == tolower(ces_il6_100$effect_allele.y), 
                         ces_il6_100$beta.x, -ces_il6_100$beta.x)
ces_il6_100$gy_se <- ces_il6_100$standard_error

mr_allmethods(mr_input(bx = ces_il6_100$gx, bxse = ces_il6_100$gx_se, 
                       by = ces_il6_100$gy, 
                       byse = ces_il6_100$gy_se), method = "all")
mr_ivw_fe(ces_il6_100$gx, ces_il6_100$gy, ces_il6_100$gx_se, ces_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
ces_il6_10$gx <- ces_il6_10$beta.y
ces_il6_10$gx_se <- ces_il6_10$se 
ces_il6_10$gy <- ifelse(tolower(ces_il6_10$effect_allele.x) == tolower(ces_il6_10$effect_allele.y), 
                        ces_il6_10$beta.x, -ces_il6_10$beta.x)
ces_il6_10$gy_se <- ces_il6_10$standard_error

mr_allmethods(mr_input(bx = ces_il6_10$gx, bxse = ces_il6_10$gx_se, 
                       by = ces_il6_10$gy, 
                       byse = ces_il6_10$gy_se), method = "all")
mr_ivw_fe(ces_il6_10$gx, ces_il6_10$gy, ces_il6_10$gx_se, ces_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(ces_il6, "results/ces_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(ces_il6_100, "results/ces_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(ces_il6_10, "results/ces_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Small Vessel Stroke
# =============================================================================
# Load small vessel stroke GWAS data from GIGASTROKE
svs <- fread("data/raw/GCST90104537_buildGRCh37_SVS.tsv.gz")
colnames(svs)[1] <- "chr_name"
colnames(svs)[2] <- "chr_start"

# Merge with genetic instruments
svs_il6 <- merge(svs, crp_il6_unf_c, by = c("chr_name", "chr_start"))
svs_il6_100 <- merge(svs, crp_il6_unf_100_c, by = c("chr_name", "chr_start"))
svs_il6_10 <- merge(svs, crp_il6_unf_10_c, by = c("chr_name", "chr_start"))
rm(svs)

# Prepare data for MR analysis - 300kb instrument
svs_il6$gx <- svs_il6$beta.y
svs_il6$gx_se <- svs_il6$se 
svs_il6$gy <- ifelse(tolower(svs_il6$effect_allele.x) == tolower(svs_il6$effect_allele.y), 
                     svs_il6$beta.x, -svs_il6$beta.x)
svs_il6$gy_se <- svs_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = svs_il6$gx, bxse = svs_il6$gx_se, 
                       by = svs_il6$gy, 
                       byse = svs_il6$gy_se), method = "all")
mr_ivw_fe(svs_il6$gx, svs_il6$gy, svs_il6$gx_se, svs_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
svs_il6_100$gx <- svs_il6_100$beta.y
svs_il6_100$gx_se <- svs_il6_100$se 
svs_il6_100$gy <- ifelse(tolower(svs_il6_100$effect_allele.x) == tolower(svs_il6_100$effect_allele.y), 
                         svs_il6_100$beta.x, -svs_il6_100$beta.x)
svs_il6_100$gy_se <- svs_il6_100$standard_error

mr_allmethods(mr_input(bx = svs_il6_100$gx, bxse = svs_il6_100$gx_se, 
                       by = svs_il6_100$gy, 
                       byse = svs_il6_100$gy_se), method = "all")
mr_ivw_fe(svs_il6_100$gx, svs_il6_100$gy, svs_il6_100$gx_se, svs_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
svs_il6_10$gx <- svs_il6_10$beta.y
svs_il6_10$gx_se <- svs_il6_10$se 
svs_il6_10$gy <- ifelse(tolower(svs_il6_10$effect_allele.x) == tolower(svs_il6_10$effect_allele.y), 
                        svs_il6_10$beta.x, -svs_il6_10$beta.x)
svs_il6_10$gy_se <- svs_il6_10$standard_error

mr_allmethods(mr_input(bx = svs_il6_10$gx, bxse = svs_il6_10$gx_se, 
                       by = svs_il6_10$gy, 
                       byse = svs_il6_10$gy_se), method = "all")
mr_ivw_fe(svs_il6_10$gx, svs_il6_10$gy, svs_il6_10$gx_se, svs_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(svs_il6, "results/svs_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(svs_il6_100, "results/svs_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(svs_il6_10, "results/svs_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Cryptogenic Stroke
# =============================================================================
# Load cryptogenic stroke GWAS data from SiGN
crypt <- fread("data/raw/stroke.meta_analysis.toastUNDETER.EUR.all.out.gz")

# Merge with genetic instruments
crypt_il6 <- merge(crypt, crp_il6_unf_c, by = "SNP")
crypt_il6_100 <- merge(crypt, crp_il6_unf_100_c, by = "SNP")
crypt_il6_10 <- merge(crypt, crp_il6_unf_10_c, by = "SNP")
rm(crypt)

# Prepare data for MR analysis - 300kb instrument
crypt_il6$gx <- crypt_il6$beta
crypt_il6$gx_se <- crypt_il6$se
crypt_il6$gy <- ifelse(tolower(crypt_il6$A1) == tolower(crypt_il6$effect_allele), 
                       crypt_il6$BETA_FIXED, -crypt_il6$BETA_FIXED)
crypt_il6$gy_se <- crypt_il6$SE_FIXED

# Run MR analysis
mr_allmethods(mr_input(bx = crypt_il6$gx, bxse = crypt_il6$gx_se, 
                       by = crypt_il6$gy, 
                       byse = crypt_il6$gy_se), method = "all")
mr_ivw_fe(crypt_il6$gx, crypt_il6$gy, crypt_il6$gx_se, crypt_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
crypt_il6_100$gx <- crypt_il6_100$beta
crypt_il6_100$gx_se <- crypt_il6_100$se
crypt_il6_100$gy <- ifelse(tolower(crypt_il6_100$A1) == tolower(crypt_il6_100$effect_allele), 
                           crypt_il6_100$BETA_FIXED, -crypt_il6_100$BETA_FIXED)
crypt_il6_100$gy_se <- crypt_il6_100$SE_FIXED

mr_allmethods(mr_input(bx = crypt_il6_100$gx, bxse = crypt_il6_100$gx_se, 
                       by = crypt_il6_100$gy, 
                       byse = crypt_il6_100$gy_se), method = "all")
mr_ivw_fe(crypt_il6_100$gx, crypt_il6_100$gy, crypt_il6_100$gx_se, crypt_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
crypt_il6_10$gx <- crypt_il6_10$beta
crypt_il6_10$gx_se <- crypt_il6_10$se
crypt_il6_10$gy <- ifelse(tolower(crypt_il6_10$A1) == tolower(crypt_il6_10$effect_allele), 
                          crypt_il6_10$BETA_FIXED, -crypt_il6_10$BETA_FIXED)
crypt_il6_10$gy_se <- crypt_il6_10$SE_FIXED

mr_allmethods(mr_input(bx = crypt_il6_10$gx, bxse = crypt_il6_10$gx_se, 
                       by = crypt_il6_10$gy, 
                       byse = crypt_il6_10$gy_se), method = "all")
mr_ivw_fe(crypt_il6_10$gx, crypt_il6_10$gy, crypt_il6_10$gx_se, crypt_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(crypt_il6, "results/crypt_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(crypt_il6_100, "results/crypt_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(crypt_il6_10, "results/crypt_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Subclinical Atherosclerosis: Carotid Plaque
# =============================================================================
# Load carotid plaque GWAS data from CHARGE
plaque <- fread("data/raw/Plaque_meta_032218.csv")
colnames(plaque)[17] <- "SNP"

# Merge with genetic instruments
plaque_il6 <- merge(plaque, crp_il6_unf_c, by = "SNP")
plaque_il6_100 <- merge(plaque, crp_il6_unf_100_c, by = "SNP")
plaque_il6_10 <- merge(plaque, crp_il6_unf_10_c, by = "SNP")
rm(plaque)

# Prepare data for MR analysis - 300kb instrument
plaque_il6$gx <- plaque_il6$beta
plaque_il6$gx_se <- plaque_il6$se
plaque_il6$gy <- ifelse(tolower(plaque_il6$Allele1) == tolower(plaque_il6$effect_allele), 
                        plaque_il6$Effect, -plaque_il6$Effect)
plaque_il6$gy_se <- plaque_il6$StdErr

# Run MR analysis
mr_allmethods(mr_input(bx = plaque_il6$gx, bxse = plaque_il6$gx_se, 
                       by = plaque_il6$gy, 
                       byse = plaque_il6$gy_se), method = "all")
mr_ivw_fe(plaque_il6$gx, plaque_il6$gy, plaque_il6$gx_se, plaque_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
plaque_il6_100$gx <- plaque_il6_100$beta
plaque_il6_100$gx_se <- plaque_il6_100$se
plaque_il6_100$gy <- ifelse(tolower(plaque_il6_100$Allele1) == tolower(plaque_il6_100$effect_allele), 
                            plaque_il6_100$Effect, -plaque_il6_100$Effect)
plaque_il6_100$gy_se <- plaque_il6_100$StdErr

mr_allmethods(mr_input(bx = plaque_il6_100$gx, bxse = plaque_il6_100$gx_se, 
                       by = plaque_il6_100$gy, 
                       byse = plaque_il6_100$gy_se), method = "all")
mr_ivw_fe(plaque_il6_100$gx, plaque_il6_100$gy, plaque_il6_100$gx_se, plaque_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
plaque_il6_10$gx <- plaque_il6_10$beta
plaque_il6_10$gx_se <- plaque_il6_10$se
plaque_il6_10$gy <- ifelse(tolower(plaque_il6_10$Allele1) == tolower(plaque_il6_10$effect_allele), 
                           plaque_il6_10$Effect, -plaque_il6_10$Effect)
plaque_il6_10$gy_se <- plaque_il6_10$StdErr

mr_allmethods(mr_input(bx = plaque_il6_10$gx, bxse = plaque_il6_10$gx_se, 
                       by = plaque_il6_10$gy, 
                       byse = plaque_il6_10$gy_se), method = "all")
mr_ivw_fe(plaque_il6_10$gx, plaque_il6_10$gy, plaque_il6_10$gx_se, plaque_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(plaque_il6, "results/plaque_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(plaque_il6_100, "results/plaque_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(plaque_il6_10, "results/plaque_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# =============================================================================
# Coronary Artery Calcification (CACS)
# =============================================================================
# Load coronary artery calcification GWAS data
cacs <- fread("data/raw/GCST90278455.tsv")
colnames(cacs)[1] <- "chr_name"
colnames(cacs)[2] <- "chr_start"

# Merge with genetic instruments
cacs_il6 <- merge(cacs, crp_il6_unf_c, by = c("chr_name", "chr_start"))
cacs_il6_100 <- merge(cacs, crp_il6_unf_100_c, by = c("chr_name", "chr_start"))
cacs_il6_10 <- merge(cacs, crp_il6_unf_10_c, by = c("chr_name", "chr_start"))
rm(cacs)

# Prepare data for MR analysis - 300kb instrument
cacs_il6$gx <- cacs_il6$beta.y
cacs_il6$gx_se <- cacs_il6$se
cacs_il6$gy <- ifelse(tolower(cacs_il6$effect_allele.x) == tolower(cacs_il6$effect_allele.y), 
                      cacs_il6$beta.x, -cacs_il6$beta.x)
cacs_il6$gy_se <- cacs_il6$standard_error

# Run MR analysis
mr_allmethods(mr_input(bx = cacs_il6$gx, bxse = cacs_il6$gx_se, 
                       by = cacs_il6$gy, 
                       byse = cacs_il6$gy_se), method = "all")
mr_ivw_fe(cacs_il6$gx, cacs_il6$gy, cacs_il6$gx_se, cacs_il6$gy_se, parameters = default_parameters())

# Repeat for 100kb instrument
cacs_il6_100$gx <- cacs_il6_100$beta.y
cacs_il6_100$gx_se <- cacs_il6_100$se
cacs_il6_100$gy <- ifelse(tolower(cacs_il6_100$effect_allele.x) == tolower(cacs_il6_100$effect_allele.y), 
                          cacs_il6_100$beta.x, -cacs_il6_100$beta.x)
cacs_il6_100$gy_se <- cacs_il6_100$standard_error

mr_allmethods(mr_input(bx = cacs_il6_100$gx, bxse = cacs_il6_100$gx_se, 
                       by = cacs_il6_100$gy, 
                       byse = cacs_il6_100$gy_se), method = "all")
mr_ivw_fe(cacs_il6_100$gx, cacs_il6_100$gy, cacs_il6_100$gx_se, cacs_il6_100$gy_se, parameters = default_parameters())

# Repeat for 10kb instrument
cacs_il6_10$gx <- cacs_il6_10$beta.y
cacs_il6_10$gx_se <- cacs_il6_10$se
cacs_il6_10$gy <- ifelse(tolower(cacs_il6_10$effect_allele.x) == tolower(cacs_il6_10$effect_allele.y), 
                         cacs_il6_10$beta.x, -cacs_il6_10$beta.x)
cacs_il6_10$gy_se <- cacs_il6_10$standard_error

mr_allmethods(mr_input(bx = cacs_il6_10$gx, bxse = cacs_il6_10$gx_se, 
                       by = cacs_il6_10$gy, 
                       byse = cacs_il6_10$gy_se), method = "all")
mr_ivw_fe(cacs_il6_10$gx, cacs_il6_10$gy, cacs_il6_10$gx_se, cacs_il6_10$gy_se, parameters = default_parameters())

# Save results
write.table(cacs_il6, "results/cacs_il6.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cacs_il6_100, "results/cacs_il6_100.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cacs_il6_10, "results/cacs_il6_10.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")