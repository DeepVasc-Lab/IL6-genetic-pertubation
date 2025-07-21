# PheWAS Mendelian Randomization Analysis - Based on Original Code
# Load required packages
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(RColorBrewer)
library(export)

# ============================================================================
# 1. MR Analysis Function
# ============================================================================
MR_function <- function(exposure.data, outcome.data) {
  # Harmonize data
  final.data <- harmonise_data(exposure.data, outcome.data, action = 1)
  final.data <- final.data[!duplicated(final.data),]
  
  # Perform MR analysis
  results <- mr(final.data, method_list = c("mr_ivw_fe", "mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  # Heterogeneity and pleiotropy tests
  heterogeneity <- mr_heterogeneity(final.data)
  pleiotropy <- mr_pleiotropy_test(final.data)
  
  # Add test results
  ivw_qpval <- heterogeneity$Q_pval[heterogeneity$method == "Inverse variance weighted"]
  results$IVW.Qpval <- ivw_qpval
  results$pleiotropy.pval <- pleiotropy$pval
  results$pleiotropy.intercept <- pleiotropy$egger_intercept
  results$pleiotropy.se <- pleiotropy$se
  
  return(results)
}

# ============================================================================
# 2. Load Exposure Data
# ============================================================================
# Load IL6 exposure
IL6 <- fread("IL6_12SNP_instrument.txt", data.table = FALSE) 

IL6_exp <- format_data(IL6, 
                       effect_allele_col = 'effect_allele', 
                       other_allele_col = 'other_allele',
                       type = 'exposure', 
                       se_col = 'se', 
                       snp_col = 'SNP', 
                       eaf_col = "eaf",
                       beta_col = 'beta', 
                       pval_col = 'pval.exposure') %>% 
  mutate(exposure = 'IL6')

# ============================================================================
# 3. Load Outcome Data
# ============================================================================
# Load manifest
finngen_r12 <- fread("finngen_R12_manifest.tsv", data.table = FALSE)

# Note: Merge phenotype-specific GWAS files that were pre-filtered on the server 
# to include only SNPs overlapping with IL6/IL6R instruments.
# This script assumes the filtered files are already saved locally as .txt/.csv.

merge_files_correctly <- function(file_list, mapping_df, source_type) {
  combined_df <- data.frame()
  
  for (file_path in file_list) {
    data <- read.csv(file_path, stringsAsFactors = FALSE, sep = "\t")
    file_name <- basename(file_path)
    file_base <- tools::file_path_sans_ext(file_name)
    phenocode <- sub("^(IL6?_finngen_R12_)(.*)$", "\\2", file_base)
    
    mapping_info <- mapping_df %>% filter(phenocode == !!phenocode)
    
    if (nrow(mapping_info) == 0) {
      warning(paste("Phenocode", phenocode, "not found"))
      next
    }
    
    data$phenocode <- mapping_info$phenocode
    data$phenotype <- mapping_info$phenotype
    data$category <- mapping_info$category
    data$source <- source_type
    
    combined_df <- bind_rows(combined_df, data)
  }
  
  return(combined_df)
}

# Get file lists
il6_files <- list.files("IL6_results/", pattern = "^IL6_finngen_R12_.*\\.(csv|txt)$", full.names = TRUE)

# Merge data
il6_combined <- merge_files_correctly(il6_files, finngen_r12, "IL6")

# ============================================================================
# 4. Format outcome data and run MR
# ============================================================================
# Format outcome data
outcome_il6 <- format_data(
  il6_combined,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  pval_col = "pval",
  effect_allele_col = "alt",
  eaf_col = "af_alt",
  other_allele_col = 'ref'
)
# Harmonize data
har_il6 <- harmonise_data(IL6_exp, outcome_il6, action = 1)

# Perform MR analysis
MR_il6 <- mr(har_il6, method_list = c("mr_ivw_fe", "mr_ivw", "mr_egger_regression", "mr_weighted_median"))

# OR MR_function()
MR_il6 = MR_function(IL6_exp, outcome_il6)

# ============================================================================
# 4. Select appropriate method and apply FDR correction
# ============================================================================
# Filter for appropriate IVW method based on heterogeneity
il6_filtered_data <- MR_il6 %>%
  group_by(exposure, outcome) %>%
  filter(if (any(method == "Inverse variance weighted" & IVW.Qpval < 0.05)) {
    method == "Inverse variance weighted"
  } else {
    method == "Inverse variance weighted (fixed effects)"
  }) %>%
  ungroup()

# Select key columns and apply FDR correction
MR_il6_ivwfe <- il6_filtered_data[, c("outcome", "exposure", "method", 
                                      "nsnp", "b", "se", "pval")] 
MR_il6_ivwfe$pval_adjust = p.adjust(MR_il6_ivwfe$pval, method = "BH")

# ============================================================================
# 5. Add disease categories and prepare for plotting
# ============================================================================
# Add categories for IL6
match_indices <- match(MR_il6_ivwfe$outcome, il6_combined$phenotype)
MR_il6_ivwfe$group_narrow <- il6_combined$category[match_indices]
MR_il6_ivwfe$neg_log10_qpvalue <- -log10(MR_il6_ivwfe$pval_adjust)

# ============================================================================
# 6. Create broad disease groups 
# ============================================================================
# Create simplified names mapping
simplified_names <- c(
  "Alcohol related diseases" = "Alcohol-related",
  "Asthma and related endpoints" = "Asthma-related",
  "Cardiometabolic endpoints" = "Cardiometabolic",
  "Comorbidities of Asthma" = "Asthma comorbidities",
  "Comorbidities of COPD" = "COPD comorbidities",
  "Comorbidities of Diabetes" = "Diabetes comorbidities",
  "Comorbidities of Gastrointestinal endpoints" = "GI comorbidities",
  "Comorbidities of Interstitial lung disease endpoints" = "Interstitial lung comorbidities",
  "Comorbidities of Neurological endpoints" = "Neurological comorbidities",
  "COPD and related endpoints" = "COPD-related",
  "Diabetes endpoints" = "Diabetes",
  "Diseases marked as autimmune origin" = "Autoimmune diseases",
  "Drug purchase endpoints" = "Drug-related",
  "Gastrointestinal endpoints" = "Gastrointestinal",
  "I Certain infectious and parasitic diseases (AB1_)" = "Infectious diseases",
  "II Neoplasms from hospital discharges (CD2_)" = "Neoplasms (hospital)",
  "II Neoplasms, from cancer register (ICD-O-3)" = "Neoplasms (register)",
  "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)" = "Blood/immune",
  "Interstitial lung disease endpoints" = "Interstitial lung",
  "IV Endocrine, nutritional and metabolic diseases (E4_)" = "Endocrine/metabolic",
  "IX Diseases of the circulatory system (I9_)" = "Circulatory system",
  "Miscellaneous, not yet classified endpoints" = "Miscellaneous",
  "Neurological endpoints" = "Neurological",
  "Other, not yet classified endpoints (same as #MISC)" = "Other classified",
  "Psychiatric endpoints from Katri Räikkönen" = "Psychiatric",
  "Quantitative endpoints" = "Quantitative",
  "Rheuma endpoints" = "Rheuma",
  "V Mental and behavioural disorders (F5_)" = "Mental/behavioral",
  "VI Diseases of the nervous system (G6_)" = "Nervous system",
  "VII Diseases of the eye and adnexa (H7_)" = "Eye disorders",
  "VIII Diseases of the ear and mastoid process (H8_)" = "Ear disorders",
  "X Diseases of the respiratory system (J10_)" = "Respiratory system",
  "XI Diseases of the digestive system (K11_)" = "Digestive system",
  "XII Diseases of the skin and subcutaneous tissue (L12_)" = "Skin disorders",
  "XIII Diseases of the musculoskeletal system and connective tissue (M13_)" = "Musculoskeletal system",
  "XIV Diseases of the genitourinary system (N14_)" = "Genitourinary",
  "XIX Injury, poisoning and certain other consequences of external causes (ST19_)" = "Injuries/poisoning",
  "XV Pregnancy, childbirth and the puerperium (O15_)" = "Pregnancy-related",
  "XVI Certain conditions originating in the perinatal period (P16_)" = "Perinatal conditions",
  "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)" = "Congenital anomalies",
  "XVIII Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R18_)" = "Symptoms/signs",
  "XXI Factors influencing health status and contact with health services (Z21_)" = "Health factors",
  "XXII Codes for special purposes (U22_)" = "Special codes"
)

# Create broad groups mapping
broad_groups <- c(
  # Respiratory diseases
  "Asthma-related" = "Respiratory diseases",
  "COPD-related" = "Respiratory diseases",
  "COPD comorbidities" = "Respiratory diseases",
  "Asthma comorbidities" = "Respiratory diseases",
  "Respiratory system" = "Respiratory diseases",
  "Interstitial lung" = "Respiratory diseases",
  "Interstitial lung comorbidities" = "Respiratory diseases",
  
  # Metabolic-related diseases
  "Diabetes" = "Endocrine/Metabolic-related diseases",
  "Diabetes comorbidities" = "Endocrine/Metabolic-related diseases",
  "Endocrine/metabolic" = "Endocrine/Metabolic-related diseases",
  
  # Mental and neurological disorders
  "Mental/behavioral" = "Mental and neurological disorders",
  "Nervous system" = "Mental and neurological disorders",
  "Psychiatric" = "Mental and neurological disorders",
  "Neurological comorbidities" = "Mental and neurological disorders",
  "Neurological" = "Mental and neurological disorders",
  
  # Cardiovascular diseases
  "Circulatory system" = "Cardiovascular diseases",
  "Cardiometabolic" = "Cardiovascular diseases",
  
  # Autoimmune and inflammatory diseases
  "Autoimmune diseases" = "Autoimmune and inflammatory diseases",
  "Rheuma" = "Autoimmune and inflammatory diseases",
  
  # Cancer and neoplasms
  "Neoplasms (hospital)" = "Cancer and neoplasms",
  "Neoplasms (register)" = "Cancer and neoplasms",
  
  # Digestive diseases
  "Gastrointestinal" = "Digestive diseases",
  "GI comorbidities" = "Digestive diseases",
  "Digestive system" = "Digestive diseases",
  
  # Genitourinary diseases
  "Genitourinary" = "Genitourinary diseases",
  
  # Dermatologic diseases
  "Skin disorders" = "Dermatologic diseases",
  
  # Injuries, Symptoms & Signs, Musculoskeletal
  "Injuries/poisoning" = "Injuries, Symptoms & Signs, Musculoskeletal",
  "Symptoms/signs" = "Injuries, Symptoms & Signs, Musculoskeletal",
  "Musculoskeletal system" = "Injuries, Symptoms & Signs, Musculoskeletal",
  
  # Infectious diseases
  "Infectious diseases" = "Infectious diseases",
  
  # Blood/immune diseases
  "Blood/immune" = "Blood/immune diseases",
  
  # Other diseases
  "Alcohol-related" = "Miscellaneous",
  "Drug-related" = "Miscellaneous",
  "Health factors" = "Miscellaneous",
  "Eye disorders" = "Eye/ear disorders",
  "Ear disorders" = "Eye/ear disorders",
  "Pregnancy-related" = "Pregnancy/Perinatal-related diseases",
  "Perinatal conditions" = "Pregnancy/Perinatal-related diseases",
  "Congenital anomalies" = "Developmental disorders",
  "Special codes" = "Miscellaneous",
  "Quantitative" = "Miscellaneous",
  "Miscellaneous" = "Miscellaneous",
  "Other classified" = "Miscellaneous"
)

# Apply mappings to IL6 data
MR_il6_ivwfe <- MR_il6_ivwfe %>%
  mutate(
    group_narrow_simplified = simplified_names[group_narrow],
    broad_group = broad_groups[group_narrow_simplified]
  )


# Check for unmapped groups
check_mapping <- function(data, data_name) {
  unmapped <- unique(data$group_narrow_simplified[is.na(data$broad_group)])
  if (length(unmapped) > 0) {
    cat(paste("Unmapped groups in", data_name, ":\n"))
    print(unmapped)
  } else {
    cat(paste("All groups mapped successfully in", data_name, "\n"))
  }
  
  # Print broad group counts
  broad_group_counts <- data %>%
    group_by(broad_group) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))
  
  cat(paste("Broad group counts for", data_name, ":\n"))
  print(broad_group_counts)
  cat("\n")
}

check_mapping(MR_il6_ivwfe, "IL6")


# ============================================================================
# 7. Manhattan Plot 
# ============================================================================
# Define disease group order 
desired_order <- c(
  "Infectious diseases", 
  "Blood/immune diseases", 
  "Digestive diseases", 
  "Mental and neurological disorders", 
  "Endocrine/Metabolic-related diseases", 
  "Eye/ear disorders", 
  "Injuries, Symptoms & Signs, Musculoskeletal", 
  "Respiratory diseases", 
  "Cardiovascular diseases", 
  "Cancer and neoplasms", 
  "Dermatologic diseases", 
  "Genitourinary diseases", 
  "Pregnancy/Perinatal-related diseases", 
  "Autoimmune and inflammatory diseases", 
  "Developmental disorders", 
  "Miscellaneous"
)

# Use your data (replace with MR_il6_ivwfe or MR_il6r_ivwfe as needed)
data_processed <- MR_il6_ivwfe %>% 
  mutate(broad_group = factor(broad_group, levels = desired_order)) %>%
  arrange(broad_group, neg_log10_qpvalue) %>%
  mutate(Chromosome = row_number())

# Create colors
broad_group_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(data_processed$broad_group)))
names(broad_group_colors) <- levels(data_processed$broad_group)

# Manhattan plot
manhattan_plot <- ggplot(data_processed, aes(
  x = Chromosome,
  y = neg_log10_qpvalue,
  fill = broad_group,
  shape = ifelse(b > 0, "up", "down")  # Using b instead of scaled_beta
)) +
  geom_point(size = 5.5, aes(color = broad_group), alpha = 1) +
  geom_hline(yintercept = -log10(0.05), color = "gray", linetype = "dashed", size = 0.5) +
  annotate("text", x = max(data_processed$Chromosome) * 0.95, 
           y = -log10(0.05), label = "FDR < 0.05", color = "black", size = 5.5, hjust = 0, vjust = -0.5) +
  geom_text_repel(
    data = data_processed %>%
      filter(pval_adjust <= 0.001),  # Label significant results, or add your custom label data
    aes(label = outcome),
    size = 4.5,
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = broad_group_colors) +
  scale_fill_manual(values = broad_group_colors) +
  scale_shape_manual(values = c("up" = 24, "down" = 25)) +
  scale_x_continuous(
    breaks = data_processed %>%
      group_by(broad_group) %>%
      summarise(center = mean(Chromosome)) %>%
      pull(center),
    labels = desired_order, 
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, max(data_processed$neg_log10_qpvalue, na.rm = TRUE) + 1),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "-log10(q-value)") +
  annotate(
    "text", x = max(data_processed$Chromosome) * 0.05, 
    y = max(data_processed$neg_log10_qpvalue, na.rm = TRUE) - 1, 
    label = "\n▲ Beta estimate > 0\n▼ Beta estimate < 0", 
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

print(manhattan_plot)
graph2ppt(x = manhattan_plot, file = "manhattan_plotil6.pptx", width = 21, height = 18)

# ============================================================================
# 8. Save results
# ============================================================================
write.xlsx(MR_il6_ivwfe, "MR_res_il6_ivwfe_finngen.xlsx")
write.xlsx(MR_il6r_ivwfe, "MR_res_il6r_ivwfe_finngen.xlsx")
write.xlsx(il6_combined, "il6_combined.xlsx")
write.xlsx(il6r_combined, "il6r_combined.xlsx")

# MR analysis code for External validation same as previous analyses, not provided here
# GWAS data sources see Supplementary Table 16
