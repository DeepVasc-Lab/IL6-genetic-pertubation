# =============================================================================
# Metabolic Outcomes Analysis (Simplified)
# =============================================================================
# Purpose: Analyze metabolic traits using IL6 genetic instruments
# =============================================================================

# NOTE: This script provides representative examples for each major outcome category.
# The complete analysis includes all traits listed below, using identical
# formatting parameters for traits from the same data source:
#
# INFECTIOUS DISEASES:
# - Any infection, Pneumonia (UK Biobank + FinnGen)
# - Urinary tract infection, Skin/soft tissue infections (UK Biobank)  
# - Candida infection (UK Biobank + FinnGen), Influenza (MVP)
# - Sepsis outcomes: admission, <75y, critical sepsis (UK Biobank + FinnGen)
# - COVID-19: vs population, hospitalized, severe respiratory (COVID19-hg)
#
# ALLERGIC ENDPOINTS:
# - Asthma (UK Biobank)
# - Atopic dermatitis (multi-cohort meta-analysis)
#
# DIABETES & METABOLIC TRAITS:
# - Type 2 diabetes (DIAGRAM/DIAMANTE/T2DGGI)
# - Glycemic traits: random glucose, HbA1c, fasting insulin/glucose, 
#   2h-glucose, insulin sensitivity (MAGIC)
# - Anthropometric: BMI (GIANT), WHR, VAT, ASAT, GFAT (UK Biobank)
#
# BLOOD CELL COUNTS (Blood Cell Consortium):
# - RBC traits: count, hemoglobin, hematocrit, MCH, MCV, MCHC, RDW
# - WBC traits: total count, neutrophils, lymphocytes, monocytes, basophils, eosinophils  
# - Platelet traits: count, mean volume
#
# For code efficiency, representative examples are shown for each data source.
# All related traits use identical formatting parameters and analytical approach.
# =============================================================================

# Load required libraries
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MendelianRandomization)

# =============================================================================
# Load and prepare genetic instruments
# =============================================================================
# Load IL6 instruments
IL6_12 <- fread("data/instruments/il6_instruments_300.txt", data.table = FALSE)
IL6_exp <- format_data(
  IL6_12,
  effect_allele_col = 'effect_allele',
  other_allele_col  = 'other_allele',
  type              = 'exposure',
  se_col            = 'se',
  snp_col           = 'SNP',
  eaf_col           = 'maf',
  beta_col          = 'beta',
  pval_col          = 'pval.exposure'
) %>% mutate(exposure = 'IL6')

# Load IL6R instruments for comparison from previous Geogarkis et al 2022 BMC Med 
IL6R_26 <- fread("data/instruments/IL6R_26SNP_instrument.csv", data.table = FALSE)
IL6R_exp <- format_data(
  IL6R_26,
  effect_allele_col = 'effect_allele',
  other_allele_col  = 'other_allele',
  type              = 'exposure',
  se_col            = 'stderr',
  snp_col           = 'SNP',
  eaf_col           = 'eaf',
  beta_col          = 'beta',
  pval_col          = 'pval'
) %>% mutate(exposure = 'IL6R')

# =============================================================================
# Helper function for MR analysis
# =============================================================================
run_mr_analysis <- function(exposure_data, outcome_file, format_params, analysis_name, binary_outcome = FALSE) {
  
  # Load outcome data
  outcome_data <- fread(outcome_file, data.table = FALSE)
  
  # Format outcome data
  outcome_format <- do.call(format_data, c(list(outcome_data, type = 'outcome'), format_params))
  
  # Harmonize data
  mydata <- harmonise_data(exposure_data, outcome_format, action = 2)
  
  # Run MR analysis
  res <- mr(mydata, method_list = c("mr_ivw_fe", "mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  
  # Scale results to 24% CRP reduction (ln(0.76) ≈ -0.2744)
  k <- log(0.76)
  
  if (binary_outcome) {
    # For binary outcomes - calculate OR and confidence intervals
    res2 <- res %>%
      mutate(
        beta_new = b * k,
        se_new   = abs(se * k),
        OR       = exp(beta_new),
        lci      = exp(beta_new - 1.96 * se_new),
        hci      = exp(beta_new + 1.96 * se_new)
      ) %>%
      select(method, nsnp, beta_new, OR, lci, hci, pval)
  } else {
    # For continuous outcomes - calculate beta and confidence intervals
    res2 <- res %>%
      mutate(
        beta_new = b * k,
        se_new   = abs(se * k),
        lci      = beta_new - 1.96 * se_new,
        hci      = beta_new + 1.96 * se_new
      ) %>%
      select(method, nsnp, beta_new, lci, hci, pval)
  }
  
  # Perform heterogeneity and pleiotropy tests
  heterogeneity_test <- mr_heterogeneity(mydata, method_list = "mr_ivw")
  pleiotropy_test <- mr_pleiotropy_test(mydata)
  
  # Save results
  result_name <- paste0("results/", analysis_name, "_", exposure_data$exposure[1])
  write.table(res2, paste0(result_name, "_scaled.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(heterogeneity_test, paste0(result_name, "_heterogeneity.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(pleiotropy_test, paste0(result_name, "_pleiotropy.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  # Print results
  cat("\n=== Results for", analysis_name, "- Exposure:", exposure_data$exposure[1], "===\n")
  print(res2, digits = 4)
  cat("\nHeterogeneity test:\n")
  print(heterogeneity_test)
  cat("\nPleiotropy test:\n")
  print(pleiotropy_test)
  
  return(list(results = res2, heterogeneity = heterogeneity_test, pleiotropy = pleiotropy_test))
}

# =============================================================================
# Type 2 Diabetes
# =============================================================================
t2d_params <- list(
  effect_allele_col = 'EffectAllele',
  other_allele_col  = 'NonEffectAllele',
  se_col            = 'SE',
  snp_col           = 'variant_id',
  eaf_col           = 'EAF',
  beta_col          = 'Beta',
  pval_col          = 'Pval'
)

# Run analysis for both IL6 and IL6R (binary outcome)
run_mr_analysis(IL6_exp, "data/raw/All_Metal_LDSC-CORR_Neff.v2.txt", t2d_params, "T2D", binary_outcome = TRUE)
run_mr_analysis(IL6R_exp, "data/raw/All_Metal_LDSC-CORR_Neff.v2.txt", t2d_params, "T2D", binary_outcome = TRUE)

# =============================================================================
# HbA1c
# =============================================================================
hba1c_params <- list(
  effect_allele_col = 'effect_allele',
  other_allele_col  = 'other_allele',
  se_col            = 'standard_error',
  snp_col           = 'variant_id',
  eaf_col           = 'effect_allele_frequency',
  beta_col          = 'beta',
  pval_col          = 'p_value'
)

# Run analysis for both IL6 and IL6R (continuous outcome)
run_mr_analysis(IL6_exp, "data/raw/MAGIC_RG_MA_transEthnic_2020Oct08.txt.gz", hba1c_params, "HbA1c", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/MAGIC_RG_MA_transEthnic_2020Oct08.txt.gz", hba1c_params, "HbA1c", binary_outcome = FALSE)

# =============================================================================
# Waist-to-Hip Ratio
# =============================================================================
whr_params <- list(
  effect_allele_col = 'Tested_Allele',
  other_allele_col  = 'Other_Allele',
  se_col            = 'SE',
  snp_col           = 'SNP',
  eaf_col           = 'Freq_Tested_Allele',
  beta_col          = 'BETA',
  pval_col          = 'P'
)

# Load WHR data and process SNP column
whr_data <- fread("data/raw/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz", data.table = FALSE)
whr_data <- whr_data %>% mutate(SNP = sub(":.*", "", SNP))  # Remove everything after first colon

# Save processed data temporarily
write.table(whr_data, "temp_whr.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

run_mr_analysis(IL6_exp, "temp_whr.txt", whr_params, "WHR", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "temp_whr.txt", whr_params, "WHR", binary_outcome = FALSE)

# Clean up temp file
file.remove("temp_whr.txt")

# =============================================================================
# Visceral Adipose Tissue
# =============================================================================
vat_params <- list(
  effect_allele_col = 'ALLELE1',
  other_allele_col  = 'ALLELE0',
  se_col            = 'SE',
  snp_col           = 'SNP',
  eaf_col           = 'A1FREQ',
  beta_col          = 'BETA',
  pval_col          = 'P_BOLT_LMM_INF'
)

run_mr_analysis(IL6_exp, "data/raw/0321_vat_bgen_stats.gz", vat_params, "VAT", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/0321_vat_bgen_stats.gz", vat_params, "VAT", binary_outcome = FALSE)

# =============================================================================
# Complete Blood Count Traits (Hematological Analysis)
# =============================================================================
# Blood Cell Consortium (BCX) data - all traits use same formatting parameters
cbc_params <- list(
  effect_allele_col = 'reference_allele',
  other_allele_col  = 'other_allele',
  se_col            = 'se',
  snp_col           = 'SNP',
  eaf_col           = 'eaf',
  beta_col          = 'beta',
  pval_col          = 'p-value'
)

# Red Blood Cell (RBC) traits - continuous outcomes
run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_RBC_Trans_GWAMA.out.gz.csv", cbc_params, "RBC_count", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_RBC_Trans_GWAMA.out.gz.csv", cbc_params, "RBC_count", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_HGB_Trans_GWAMA.out.gz.csv", cbc_params, "Hemoglobin", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_HGB_Trans_GWAMA.out.gz.csv", cbc_params, "Hemoglobin", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_HCT_Trans_GWAMA.out.gz.csv", cbc_params, "Hematocrit", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_HCT_Trans_GWAMA.out.gz.csv", cbc_params, "Hematocrit", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_MCH_Trans_GWAMA.out.gz.csv", cbc_params, "MCH", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_MCH_Trans_GWAMA.out.gz.csv", cbc_params, "MCH", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_MCV_Trans_GWAMA.out.gz.csv", cbc_params, "MCV", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_MCV_Trans_GWAMA.out.gz.csv", cbc_params, "MCV", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_MCHC_Trans_GWAMA.out.gz.csv", cbc_params, "MCHC", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_MCHC_Trans_GWAMA.out.gz.csv", cbc_params, "MCHC", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_RDW_Trans_GWAMA.out.gz.csv", cbc_params, "RDW", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_RDW_Trans_GWAMA.out.gz.csv", cbc_params, "RDW", binary_outcome = FALSE)

# White Blood Cell (WBC) traits - continuous outcomes
run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_WBC_Trans_GWAMA.out.gz.csv", cbc_params, "WBC_count", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_WBC_Trans_GWAMA.out.gz.csv", cbc_params, "WBC_count", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_NEU_Trans_GWAMA.out.gz.csv", cbc_params, "Neutrophils", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_NEU_Trans_GWAMA.out.gz.csv", cbc_params, "Neutrophils", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_LYM_Trans_GWAMA.out.gz.csv", cbc_params, "Lymphocytes", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_LYM_Trans_GWAMA.out.gz.csv", cbc_params, "Lymphocytes", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_MON_Trans_GWAMA.out.gz.csv", cbc_params, "Monocytes", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_MON_Trans_GWAMA.out.gz.csv", cbc_params, "Monocytes", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_BAS_Trans_GWAMA.out.gz.csv", cbc_params, "Basophils", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_BAS_Trans_GWAMA.out.gz.csv", cbc_params, "Basophils", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_EOS_Trans_GWAMA.out.gz.csv", cbc_params, "Eosinophils", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_EOS_Trans_GWAMA.out.gz.csv", cbc_params, "Eosinophils", binary_outcome = FALSE)

# Platelet traits - continuous outcomes
run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_PLT_Trans_GWAMA.out.gz.csv", cbc_params, "Platelet_count", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_PLT_Trans_GWAMA.out.gz.csv", cbc_params, "Platelet_count", binary_outcome = FALSE)

run_mr_analysis(IL6_exp, "data/raw/processed_BCX2_MPV_Trans_GWAMA.out.gz.csv", cbc_params, "MPV", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/processed_BCX2_MPV_Trans_GWAMA.out.gz.csv", cbc_params, "MPV", binary_outcome = FALSE)

# =============================================================================
# COVID-19 Hospitalization
# =============================================================================
covid_params <- list(
  effect_allele_col = 'ALT',
  other_allele_col  = 'REF',
  se_col            = 'all_inv_var_meta_sebeta',
  snp_col           = 'rsid',
  eaf_col           = 'all_meta_AF',
  beta_col          = 'all_inv_var_meta_beta',
  pval_col          = 'all_inv_var_meta_p'
)

# Infectious diseases - binary outcomes
run_mr_analysis(IL6_exp, "data/raw/COVID19_HGI_A2_ALL_leave_23andme_20220403_GRCh37.tsv.gz", covid_params, "COVID19", binary_outcome = TRUE)
run_mr_analysis(IL6R_exp, "data/raw/COVID19_HGI_A2_ALL_leave_23andme_20220403_GRCh37.tsv.gz", covid_params, "COVID19", binary_outcome = TRUE)

# =============================================================================
# Sepsis
# =============================================================================
sepsis_params <- list(
  effect_allele_col = 'effect_allele',
  other_allele_col  = 'other_allele',
  se_col            = 'standard_error',
  snp_col           = 'rs_id',
  eaf_col           = 'effect_allele_frequency',
  beta_col          = 'beta',
  pval_col          = 'p_value'
)

# Sepsis - binary outcome
run_mr_analysis(IL6_exp, "data/raw/GCST90281174.tsv.gz", sepsis_params, "Sepsis", binary_outcome = TRUE)
run_mr_analysis(IL6R_exp, "data/raw/GCST90281174.tsv.gz", sepsis_params, "Sepsis", binary_outcome = TRUE)

# =============================================================================
# Carotid Intima-Media Thickness (cIMT)
# =============================================================================
cimt_params <- list(
  effect_allele_col = 'effect_allele',
  other_allele_col  = 'other_allele',
  se_col            = 'standard_error',
  snp_col           = 'marker',
  eaf_col           = 'effect_allele_frequency',
  beta_col          = 'beta',
  pval_col          = 'p_value'
)

# Subclinical atherosclerosis - continuous outcome
run_mr_analysis(IL6_exp, "data/raw/GCST90100575_buildGRCh37.tsv", cimt_params, "cIMT", binary_outcome = FALSE)
run_mr_analysis(IL6R_exp, "data/raw/GCST90100575_buildGRCh37.tsv", cimt_params, "cIMT", binary_outcome = FALSE)


# =============================================================================
# Summary Results and Multiple Testing Correction
# =============================================================================
# Purpose: Combine all MR results and apply multiple testing correction
#          for IVW fixed effects analysis across all outcomes

library(data.table)
library(dplyr)

# Read data and process
data1 <- read.table("res/binary_oucomes.txt", header=T, sep="\t")
data2 <- read.table("res/continuous_oucomes.txt", header=T, sep="\t")

# Extract IVW fixed-effects only
ivw_fe1 <- data1[data1$Method == "Fixed-effects IVW",]
ivw_fe2 <- data2[data2$Method == "Fixed-effects IVW",]

# FDR correction by panel
ivw_fe1 = ivw_fe1 %>% group_by(Panel) %>% mutate(fdr = p.adjust(p.value, "BH"))
ivw_fe1 = ivw_fe1 %>% group_by(Panel) %>% mutate(fdr = p.adjust(P.value, "BH"))

# Save results
write.csv(ivw_fe1, "res/binary_results.csv")
write.csv(ivw_fe1, "res/continuous_results.csv")

# =============================================================================
# Metabolites 
# =============================================================================
# Load required packages
library(TwoSampleMR)   # For Mendelian Randomization analysis
library(data.table)    # Data manipulation
library(readr)         # Read CSV files
library(dplyr)         # Data processing
library(pacman)        # Package management
library(yulab.utils)   # Utility tools
library(ieugwasr)      # GWAS data access
library(tidyr)         # Data tidying

# Main MR function
MR_function <- function(exposure.data, outcome.data, outcome_name, adjustment_factor) {
  # Format outcome data
  outcome_data <- outcome.data %>%
    rename(
      SNP = "SNP",                               
      beta.outcome = "beta",                     
      se.outcome = "standard_error",             
      effect_allele.outcome = "effect_allele",   
      other_allele.outcome = "other_allele",     
      eaf.outcome = "effect_allele_frequency",   
      pval.outcome = "P_value"                   
    )
  
  # Add required columns for TwoSampleMR
  outcome_data$id.outcome <- outcome_name
  outcome_data$outcome <- outcome_name
  
  # Format data for TwoSampleMR
  outcome_data <- format_data(
    outcome_data,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    eaf_col = "eaf.outcome",
    pval_col = "pval.outcome"
  )
  
  # Harmonize exposure and outcome data
  final.data <- harmonise_data(exposure.data, outcome_data, action = 1)
  final.data <- final.data[which(!duplicated(final.data)),]  # Remove duplicates
  
  # Perform MR analysis (IVW methods only)
  results <- mr(final.data, method_list = c("mr_ivw", "mr_ivw_fe"))
  
  # Adjust beta and se for continuous variables
  k <- log(adjustment_factor)  # Calculate adjustment coefficient
  results <- results %>%
    mutate(
      beta_adjusted = b * k,           # Adjusted beta
      se_adjusted = abs(se * k),       # Adjusted se (keep positive)
      LCI = beta_adjusted - 1.96 * se_adjusted,  # Lower confidence interval
      UCI = beta_adjusted + 1.96 * se_adjusted   # Upper confidence interval
    )
  
  # Add SNP count
  results$SNPs <- nrow(final.data)
  
  # Heterogeneity test
  heterogeneity <- mr_heterogeneity(final.data)
  
  # Pleiotropy test
  pleiotropy <- mr_pleiotropy_test(final.data)
  
  # Add heterogeneity and pleiotropy results
  results$IVW.Qpval <- heterogeneity$Q_pval[heterogeneity$method == "Inverse variance weighted"]
  results$pleiotropy.pval <- pleiotropy$pval
  
  return(results)
}

# Set input parameters (modify these paths for your data)
data_path <- "data/metabolites/"  # Path to metabolite files
exposure_data <- IL6_exp  # Or IL6R_exp for IL-6R analysis

# Get all metabolite file paths
outcome_files <- list.files(
  path = data_path,
  pattern = "*.csv",
  full.names = TRUE
)

# Create empty list to store results
all_results <- list()

# Set adjustment factor
adjustment_factor <- log(0.76)

# Loop through each metabolite file and perform MR analysis
for (file in outcome_files) {
  tryCatch({
    outcome_data <- read_csv(file, show_col_types = FALSE)  # Read outcome data
    outcome_name <- tools::file_path_sans_ext(basename(file))  # Get metabolite name
    
    # Run MR analysis
    result <- MR_function(exposure_data, outcome_data, outcome_name, adjustment_factor)
    
    # Mark result source
    result$outcome <- outcome_name
    all_results[[outcome_name]] <- result
    
    cat("Processed:", outcome_name, "\n")
  }, error = function(e) {
    cat("Error processing", basename(file), ":", e$message, "\n")
  })
}

# Combine all analysis results
final_results <- do.call(rbind, all_results)

# Select appropriate method for each metabolite
selected_results <- final_results %>%
  group_by(outcome) %>%
  mutate(
    # Select method based on heterogeneity
    selected_method = ifelse(IVW.Qpval < 0.05, "Inverse variance weighted", "Inverse variance weighted (fixed effects)")
  ) %>%
  # Keep only selected method
  filter(method == selected_method) %>%
  ungroup()

# Perform FDR correction
selected_results$pval_fdr <- p.adjust(selected_results$pval, method = "fdr")

# Organize final results
final_output <- selected_results %>%
  select(
    exposure, outcome, method, nsnp, 
    beta_original = b, se_original = se,
    beta_adjusted, se_adjusted, LCI, UCI, 
    pval, pval_fdr, 
    SNPs, IVW.Qpval, pleiotropy.pval
  ) 

# View results
print(final_output)

# Save results
output_file <- paste0("results/", deparse(substitute(exposure_data)), "_metabolite_MR_results.csv")
write_csv(final_output, output_file)
