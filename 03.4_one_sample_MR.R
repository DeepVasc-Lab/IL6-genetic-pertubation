# One-Sample Mendelian Randomization Analysis using Cox Regression
# Analysis: IL-6 GRS → logCRP → Cardiovascular Disease Risk
# Age as time-scale with left truncation

library(dplyr)
library(data.table)
library(openxlsx)
library(survival)
# plink command for GRS calculation:

#plink2  --out il6_12var_grs --pfile /cluster/UKB_IMP/IMP_v3/chr7 --score snps_12var.txt 1 4 6 cols=+scoresums no-mean-imputation ignore-dup-ids

#1st column – rsID
#4th column – effect allele
#6th column - beta estimate
# ============================================================================
# Configuration
# ============================================================================
outcomes_file <- "data/UKB_CAD_PAD_IS_survival_ready.xlsx"
grs_file <- "data/il6_12var_grs.sscore"
output_file <- "results/Cox_MR_Results_OneSample.xlsx"
scaling_factor <- log(0.76)  # 24% CRP reduction

# ============================================================================
# 1. Data Loading and Processing
# ============================================================================
# Load survival outcomes
ukb_outcomes <- read.xlsx(outcomes_file)

# Load IL-6 genetic risk score
il6_grs <- fread(grs_file, data.table = FALSE)

# Map Munich study IDs to UK Biobank participant IDs
analysis <- il6_grs %>%
  rename(eid_munich = IID) %>%
  inner_join(bridge_new, by = "eid_munich") %>%
  rename(eid = eid_151281)

# Standardize GRS
analysis$il6_grs_std <- scale(analysis$SCORE1_SUM)
il6_grs <- analysis[, c("eid", "il6_grs_std")]

# Process covariates
pcs_data <- PCA_10 %>%
  setNames(c("eid", paste0("PC", 1:10)))

kinship_processed <- kinship %>%
  mutate(
    kinship_status = factor(p22021),
    chip_batch = factor(p22000)
  ) %>%
  select(eid, kinship_status, chip_batch)

risk_factors_processed <- risk_factors_ukb151281 %>%
  mutate(
    logCRP = scale(log(CRP)),
    blage = Age,
    sex = Sex
  ) %>%
  select(eid, blage, sex, logCRP, CRP)

# ============================================================================
# 2. Prepare Survival Data
# ============================================================================

# Merge all data
mr_data_base <- ukb_outcomes_clean %>%
  inner_join(il6_grs, by = "eid") %>%
  inner_join(risk_factors_processed, by = "eid") %>%
  inner_join(kinship_processed, by = "eid") %>%
  inner_join(pcs_data, by = "eid")

# Complete case analysis
pc_cols <- paste0("PC", 1:10)
key_vars <- c("il6_grs_std", "logCRP", "sex", "kinship_status", "chip_batch", pc_cols)

mr_data_complete <- mr_data_base %>%
  filter(complete.cases(select(., all_of(key_vars))))

# Calculate exit ages for age-scale analysis
mr_data_age <- mr_data_complete %>%
  mutate(
    age_exit_cad = ifelse(CAD == 1, CAD_age_at_event, blage + CAD_followup_years),
    age_exit_pad = ifelse(PAD == 1, PAD_age_at_event, blage + PAD_followup_years),
    age_exit_is = ifelse(IS == 1, IS_age_at_event, blage + IS_followup_years)
  )

# ============================================================================
# 3. First-Stage: GRS → logCRP
# ============================================================================
formula_linear <- as.formula(paste(
  "logCRP ~ il6_grs_std + sex +",
  paste(pc_cols, collapse = " + "),
  "+ kinship_status + chip_batch"))

linear_model <- lm(formula_linear, data = mr_data_complete)
beta_crp <- coef(linear_model)["il6_grs_std"]
se_crp <- summary(linear_model)$coefficients["il6_grs_std", "Std. Error"]

# ============================================================================
# 4. Second-Stage: Cox Regression (Age as Time-Scale)
# ============================================================================
cox_cov <- paste("il6_grs_std + sex +",
                 paste(pc_cols, collapse = " + "),
                 "+ kinship_status + chip_batch")

# Cox models with age as time-scale (left truncation)
cox_cad_age <- coxph(as.formula(paste("Surv(blage, age_exit_cad, CAD) ~", cox_cov)),
                     data = mr_data_age)
cox_pad_age <- coxph(as.formula(paste("Surv(blage, age_exit_pad, PAD) ~", cox_cov)),
                     data = mr_data_age)
cox_is_age <- coxph(as.formula(paste("Surv(blage, age_exit_is, IS) ~", cox_cov)),
                    data = mr_data_age)

# ============================================================================
# 5. Extract coefficients
# ============================================================================
# CAD
beta_cad_age <- coef(cox_cad_age)["il6_grs_std"]
se_cad_age <- summary(cox_cad_age)$coefficients["il6_grs_std", "se(coef)"]
p_cad_age <- summary(cox_cad_age)$coefficients["il6_grs_std", "Pr(>|z|)"]

# PAD
beta_pad_age <- coef(cox_pad_age)["il6_grs_std"]
se_pad_age <- summary(cox_pad_age)$coefficients["il6_grs_std", "se(coef)"]
p_pad_age <- summary(cox_pad_age)$coefficients["il6_grs_std", "Pr(>|z|)"]

# IS
beta_is_age <- coef(cox_is_age)["il6_grs_std"]
se_is_age <- summary(cox_is_age)$coefficients["il6_grs_std", "se(coef)"]
p_is_age <- summary(cox_is_age)$coefficients["il6_grs_std", "Pr(>|z|)"]

# ============================================================================
# 6. Calculate MR Estimates (Age Scale Approach)
# ============================================================================
# MR estimates: ratio of coefficients (GRS→Disease / GRS→CRP)
mr_cad_age <- beta_cad_age / beta_crp
mr_pad_age <- beta_pad_age / beta_crp
mr_is_age <- beta_is_age / beta_crp

# Calculate standard errors using Delta method
mr_se_cad_age <- abs(mr_cad_age) * sqrt((se_cad_age/abs(beta_cad_age))^2 + (se_crp/abs(beta_crp))^2)
mr_se_pad_age <- abs(mr_pad_age) * sqrt((se_pad_age/abs(beta_pad_age))^2 + (se_crp/abs(beta_crp))^2)
mr_se_is_age <- abs(mr_is_age) * sqrt((se_is_age/abs(beta_is_age))^2 + (se_crp/abs(beta_crp))^2)

# ============================================================================
# 7. Scale by ln(0.76) for Clinical Interpretation
# ============================================================================
scaling_factor <- log(0.76)

# Scaled MR estimates for 24% CRP reduction
mr_cad_age_scaled <- mr_cad_age * scaling_factor
mr_pad_age_scaled <- mr_pad_age * scaling_factor
mr_is_age_scaled <- mr_is_age * scaling_factor

mr_se_cad_age_scaled <- mr_se_cad_age * abs(scaling_factor)
mr_se_pad_age_scaled <- mr_se_pad_age * abs(scaling_factor)
mr_se_is_age_scaled <- mr_se_is_age * abs(scaling_factor)

# ============================================================================
# 8. Results Table
# ============================================================================
results <- data.frame(
  Outcome = c("CAD", "PAD", "IS"),
  Cases = c(sum(mr_data_age$CAD), sum(mr_data_age$PAD), sum(mr_data_age$IS)),
  N_total = rep(nrow(mr_data_age), 3),
  
  # Hazard ratios for 24% CRP reduction
  MR_HR_24percent = round(c(exp(mr_cad_age_scaled), exp(mr_pad_age_scaled), exp(mr_is_age_scaled)), 4),
  
  # Lower confidence intervals
  MR_HR_LCI = round(c(exp(mr_cad_age_scaled - 1.96 * mr_se_cad_age_scaled),
                      exp(mr_pad_age_scaled - 1.96 * mr_se_pad_age_scaled),
                      exp(mr_is_age_scaled - 1.96 * mr_se_is_age_scaled)), 4),
  
  # Upper confidence intervals  
  MR_HR_UCI = round(c(exp(mr_cad_age_scaled + 1.96 * mr_se_cad_age_scaled),
                      exp(mr_pad_age_scaled + 1.96 * mr_se_pad_age_scaled),
                      exp(mr_is_age_scaled + 1.96 * mr_se_is_age_scaled)), 4)
)

# Create formatted CI column
results$CI_95 <- paste0("(", results$MR_HR_LCI, "-", results$MR_HR_UCI, ")")

print(results)
write.xlsx(results, output_file, rowNames = FALSE)