# Instrument Strength Analysis - F-statistics Calculation
# Calculate R² and F-statistics for instrumental variables

# 1. Data preparation
instrument_bbj_clean <- il6_crp_instru %>% 
  filter(!is.na(SNP)) %>% 
  rename(
    SE   = Standard,
    MAF  = `Minor Allele Frequency`,
    Pval = `P-value`
  ) %>% 
  mutate(across(c(Beta, SE, MAF, Pval), as.numeric))

# 2. Sample size parameters
N_total <- 575531   # Total sample size
K       <- 1        # Number of instruments in regression

# 3. Calculate R² and F-statistics for each SNP
instrument_bbj_metrics <- instrument_bbj_clean %>% 
  rowwise() %>% 
  mutate(
    # Calculate R² using standard formula
    R2 = 2 * Beta^2 * MAF * (1 - MAF) /
      (2 * Beta^2 * MAF * (1 - MAF) + 2 * SE^2 * N_total * MAF * (1 - MAF)),
    
    # Calculate F-statistic
    F_stat = (R2 / (1 - R2)) * (N_total - K - 1) / K
  ) %>% 
  ungroup()

# 4. View results
instrument_bbj_metrics %>% 
  select(SNP, Beta, SE, MAF, R2, F_stat) %>% 
  print(n = Inf)
