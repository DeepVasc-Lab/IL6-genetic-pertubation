library(dplyr)
library(coloc)
library(data.table)
library(TwoSampleMR)
library(readxl)


#### CRP data preparation #### 
crp = fread('/Users/a1/projects/IL6_pertrubation_coloc/summ_stats/35459240-GCST90029070-EFO_0004458-Build37.f.tsv.gz')
crp = crp[crp$chromosome ==7,]
crp = crp[(crp$base_pair_location>=  22466819 )& (crp$base_pair_location<=23071617),]

crp = crp[(!is.na(crp$variant_id)),]
crp = crp[!duplicated(crp$variant_id), ]

crp$chr_pos = paste0(crp$chromosome,'_',crp$base_pair_location)

freqs = fread('./crp_freq.frq')
colnames(freqs)[2] = 'variant_id'
colnames(freqs)[5] = 'maf' # for SdY 

freqs = freqs[freqs$maf!=0,]
freqs = freqs[!duplicated(freqs$variant_id), ]

crp =merge(crp,freqs[,c('variant_id','maf')],by = 'variant_id')

il6_form = format_data(as.data.frame(crp), effect_allele_col ='effect_allele', other_allele_col = 'other_allele',
                       type = 'exposure',se_col = 'standard_error',snp_col = 'variant_id',pos_col = 'base_pair_location',
                       beta_col = 'beta', pval_col = 'p_value',eaf_col = 'maf')


#### Function for colocalisation ####
perform_coloc <- function(harmonised, N) {
  DatasetCrp <- list(
    beta = harmonised$beta.exposure,
    varbeta = harmonised$se.exposure^2,
    snp = harmonised$SNP,
    position = harmonised$pos.exposure,
    type = "quant",
    N = as.integer(575531),
    MAF = harmonised$eaf.exposure
  )
  
  DatasetOutcome <- list(
    beta = harmonised$beta.outcome,
    varbeta = harmonised$se.outcome^2,
    snp = harmonised$SNP,
    position = harmonised$pos.outcome,
    type = "cc",
    N = as.integer(N),
    MAF = harmonised$eaf.outcome
  )
  
  c_sum <- coloc::coloc.abf(DatasetCrp, DatasetOutcome)
  
  return(as.data.frame(as.list(c_sum$summary)))
}


### Coronary Artery Disease #### 

cad = fread('../summ_stats/CAD_chr7_GCST90132315.tsv')
cad$chr_pos = paste0(cad$chromosome,'_',cad$base_pair_location)


cad = merge(crp[,c('variant_id','chr_pos')],cad,by = 'chr_pos' )
cad = cad[!duplicated(cad$variant_id), ]


cad_form = format_data(as.data.frame(cad), effect_allele_col ='effect_allele', other_allele_col = 'other_allele',
                       type = 'outcome',se_col = 'standard_error',snp_col = 'variant_id',pos_col = 'base_pair_location',
                       beta_col = 'beta', pval_col = 'p_value',eaf_col = 'minfreq')


for_coloc_cad <- harmonise_data(
  exposure_dat = il6_form, 
  outcome_dat = cad_form
)

cad_res = perform_coloc(for_coloc_cad, 1378170)



#### Ischemic Stroke ####
isstroke = fread('../summ_stats/IS_GCST90104535_buildGRCh37.tsv.gz')
head(isstroke)
isstroke = isstroke[isstroke$chromosome==7,]
isstroke$chr_pos = paste0(isstroke$chromosome,'_',isstroke$base_pair_location)

rm(isstroke)

isstroke = merge(crp[,c('variant_id','chr_pos')],isstroke,by = 'chr_pos' )
isstroke = isstroke[!duplicated(isstroke$variant_id), ]

isstroke_form = format_data(as.data.frame(isstroke), effect_allele_col ='effect_allele', other_allele_col = 'other_allele',
                       type = 'outcome',se_col = 'standard_error',snp_col = 'variant_id',pos_col = 'base_pair_location',
                       beta_col = 'beta', pval_col = 'p_value',eaf_col = 'effect_allele_frequency')


for_coloc_isstroke <- harmonise_data(
  exposure_dat = il6_form, 
  outcome_dat = isstroke_form
)


iss_res = perform_coloc(for_coloc_isstroke, 1590566)


#### Large Artery Stroke ####

las = fread('../summ_stats/LAS.tsv.gz')

las =las[las$chromosome==7,]

las$chr_pos = paste0(las$chromosome,'_',las$base_pair_location)
las = merge(crp[,c('variant_id','chr_pos')],las,by = 'chr_pos' )
las = las[!duplicated(las$variant_id), ]


las_form = format_data(as.data.frame(las), effect_allele_col ='effect_allele', other_allele_col = 'other_allele',
                            type = 'outcome',se_col = 'standard_error',snp_col = 'variant_id',pos_col = 'base_pair_location',
                            beta_col = 'beta', pval_col = 'p_value',eaf_col = 'effect_allele_frequency')


for_coloc_las <- harmonise_data(
  exposure_dat = il6_form, 
  outcome_dat = las_form
)


las_res = perform_coloc(for_coloc_las,1503898+9219)

#### Carotid Plaque ####
plaque = fread('../summ_stats/plaque.tsv')

plaque = plaque[!duplicated(plaque$variant_id), ]

plaque = plaque[plaque$rs_id %in% crp$variant_id,]
plaque_form = format_data(as.data.frame(plaque), effect_allele_col ='effect_allele', other_allele_col = 'other_allele',
                       type = 'outcome',se_col = 'standard_error',snp_col = 'rs_id',pos_col = 'base_pair_location',
                       beta_col = 'beta', pval_col = 'p_value',eaf_col = 'effect_allele_frequency')


for_coloc_plaque <- harmonise_data(
  exposure_dat = il6_form, 
  outcome_dat = plaque_form
)



pl_res = perform_coloc(for_coloc_plaque, 48434)



#### Rheumatoid Arthritis ####

ra = fread('../summ_stats/GCST90132222_buildGRCh37.tsv.gz')

ra = ra[ra$chromosome==7,]


ra = ra[!duplicated(ra$variant_id), ]
ra = ra[!is.na(ra$variant_id), ]
ra=ra[ra$variant_id %in% crp$variant_id,]
ra_form = format_data(as.data.frame(ra), effect_allele_col ='effect_allele', other_allele_col = 'other_allele',
                          type = 'outcome',se_col = 'standard_error',snp_col = 'variant_id',pos_col = 'base_pair_location',
                          beta_col = 'beta', pval_col = 'p_value',eaf_col = 'effect_allele_frequency')
head(ra)

for_coloc_form <- harmonise_data(
  exposure_dat = il6_form, 
  outcome_dat = ra_form
)
for_coloc_form=for_coloc_form[for_coloc_form$mr_keep==T,]
for_coloc_form$eaf.outcome = for_coloc_form$eaf.exposure


ra_res = perform_coloc(for_coloc_form, 35871 + 240149)


#### Polymyalgia Rheumatica ####
pR = fread('../summ_stats/Meta_PMR.tsv')
pR = pR[pR$chr==7,]

head(pR)
pR_form = format_data(as.data.frame(pR), effect_allele_col ='EA', other_allele_col = 'NEA',
                          type = 'outcome',se_col = 'StdErr',snp_col = 'SNP',pos_col = 'position',
                          beta_col = 'Beta', pval_col = 'Pvalue',eaf_col = 'MAF')
rm(pR)

for_coloc_pR <- harmonise_data(
  exposure_dat = il6_form, 
  outcome_dat = pR_form
)
pr_res = perform_coloc(for_coloc_pR, 8156+416495)


#### Type 2 Diabetes ####

t2d = fread('../summ_stats/All_Metal_LDSC-CORR_Neff.v2.txt')
t2d=t2d[t2d$Chromsome==7,]
t2d$chr_pos = paste0(t2d$Chromsome,'_',t2d$Position)
t2d = merge(crp[,c('variant_id','chr_pos')],t2d,by = 'chr_pos' )

t2d_form = format_data(as.data.frame(t2d), effect_allele_col ='EffectAllele', other_allele_col = 'NonEffectAllele',
                       type = 'outcome',se_col = 'SE',snp_col = 'variant_id',pos_col = 'Position',
                       beta_col = 'Beta', pval_col = 'PvalAssociation',eaf_col = 'EAF')

for_coloc_t2d<- harmonise_data(
  exposure_dat = il6_form, 
  outcome_dat = t2d_form
)
t2d_res = perform_coloc(for_coloc_t2d, 428452 + 2107149)



#### Peripheral Artery Disease ####
PAD = fread('../summ_stats/PAD.csv')
pad_form = format_data(as.data.frame(PAD), effect_allele_col ='effect_allele.outcome', other_allele_col = 'other_allele.outcome',
                       type = 'outcome',se_col = 'se.outcome',snp_col = 'SNP',pos_col = 'pos.outcome',
                       beta_col = 'beta.outcome', pval_col = 'pval.outcome',eaf_col = 'eaf.outcome')

for_coloc_pad<- harmonise_data(
  exposure_dat = il6_form, 
  
  outcome_dat = pad_form
)

pad_res  = perform_coloc(for_coloc_pad, 31307 +211753)


#### Plot datasets #### 

DatasetCrp <- list(
  beta = il6_form$beta.exposure,
  varbeta = il6_form$se.exposure^2,
  snp = il6_form$SNP,
  position = il6_form$pos.exposure,
  type = "quant",
  N = as.integer(575531),
  MAF = il6_form$eaf.exposure
)

DatasetCad =  list(beta=for_coloc_cad$beta.outcome,varbeta = for_coloc_cad$se.outcome**2,
                   snp = for_coloc_cad$SNP, N = 1378170,
                   position=for_coloc_cad$pos.outcome,type ='cc')
  
  
  
DatasetIscStroke= list(beta=for_coloc_isstroke$beta.outcome,varbeta = for_coloc_isstroke$se.outcome**2,
                 snp = for_coloc_isstroke$SNP, N = 1590566,
                 position=for_coloc_isstroke$pos.outcome,type ='cc')


DatasetLas = list(beta=for_coloc_las$beta.outcome,varbeta = for_coloc_las$se.outcome**2,
                  snp = for_coloc_las$SNP, N = 1590566,
                  position=for_coloc_las$pos.outcome,type ='cc')

DatasetPAD =  list(beta=for_coloc_pad$beta.outcome,varbeta = for_coloc_pad$se.outcome**2,
                   snp = for_coloc_pad$SNP, N = 31307 +211753,
                   position=for_coloc_pad$pos.outcome,type ='cc')

DatasetPlaque  =  list(beta=for_coloc_plaque$beta.outcome,varbeta = for_coloc_plaque$se.outcome**2,
                       snp = for_coloc_plaque$SNP, N = 48434,
                       position=for_coloc_plaque$pos.outcome,type ='cc')


DatasetT2d= list(beta=for_coloc_t2d$beta.outcome,varbeta = for_coloc_t2d$se.outcome**2,
                 snp = for_coloc_t2d$SNP, N = 428452 + 2107149 ,
                 position=for_coloc_t2d$pos.outcome,type ='cc')

DatasetRa= list(beta=for_coloc_form$beta.outcome,varbeta = for_coloc_form$se.outcome**2,
                 snp = for_coloc_form$SNP, N = 35871 + 240149 ,
                 position=for_coloc_form$pos.outcome,type ='cc')


DatasetpR = list(beta=for_coloc_pR$beta.outcome,varbeta = for_coloc_pR$se.outcome**2,
                 snp = for_coloc_pR$SNP, N = 8156+416495,
                 position=for_coloc_pR$pos.outcome,type ='cc')

  
  
  
datasets <- list(
  "Coronary artery disease" = DatasetCad,
  "Ischemic stroke" = DatasetIscStroke,
  "Large artery stroke" = DatasetLas,
  "Peripheral artery disease" = DatasetPAD,
  "Carotid plaque" = DatasetPlaque,
  "Type 2 diabetes" = DatasetT2d,
  "Rheumatoid arthritis" = DatasetRa,
  "Polymyalgia rheumatica" = DatasetpR
)

## Extract the variants from IL6 instrument
var12= read_excel('~/il6_instruments.xlsx',sheet = '12var')


pdf('figure_coloc_locus_plots.pdf',width = 12,height = 8)

par(mfrow = c(3, 4))  

plot_dataset(DatasetCrp,highlight_list = var12$snp,color = rep("dodgerblue2",12),show_legend = F)
for (i in 1:3) {
  plot.new()  
}

title(main = "CRP")


for (dataset_name in names(datasets)) {
  plot_dataset(datasets[[dataset_name]],highlight_list = var12$snp,color = rep("dodgerblue2",12),show_legend = F)  # Call the function
  title(main = dataset_name)  # Add the dataset name as the title
}


dev.off()



