# IL6 genetic perturbation mimicking IL-6 inhibition is associated with lower cardiometabolic risk 
This repository contains the analysis code for our study on genetic perturbation of IL6 mimicking IL-6 inhibition and its association with cardiometabolic risk. We developed a robust genetic proxy for IL-6 signaling downregulation and used it to predict the effects of IL-6 inhibition on cardiovascular and safety endpoints.

## Analysis Workflow

### 1. Construction of IL6 Genetic Instrument
Build genetic instruments using variants in the IL6 locus associated with CRP levels

### 2. Concordance with Pharmacological IL-6 Inhibition
Validate instruments against clinical trial biomarkers and autoimmune outcomes

### 3. Effects on Cardiovascular Endpoints
Test associations with cardiovascular diseases and perform colocalization analysis

### 4. Effects on Metabolic Outcomes and Traits
Analyze type 2 diabetes, lipid profiles, and metabolic markers

### 5. Effects on Key Safety Endpoints
Assess infection risk, hematological traits, and safety profiles

### 6. Phenome-wide Association Study
Comprehensive analysis across 2,469 clinical outcomes in FinnGen

## Requirements
- R ≥ 4.4.3; Python(3.9.1);PLINK (v2.00a3.3LM) 
- Required packages: TwoSampleMR, MendelianRandomization, coloc, data.table, dplyr, metafor

## Data Sources
GWAS summary statistics from public consortia. See Supplementary Table 1 in the manuscript for detailed information on all datasets used.

## Contact
For questions about the analysis or data, please contact:
- Lanyue Zhang - Lanyue.Zhang@med.uni-muenchen.de / zhanglanyue1996@gmail.com
