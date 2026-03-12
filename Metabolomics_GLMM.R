# Environment -------------------------------------------------------------
library(tidyverse)
library(lmerTest)
library(here)
library(magrittr)
setwd('~/Documents/Projects/Short_Term_Broccoli/')
here::i_am('./human_analysis/Untargeted_MultiOmics/Biomarker_Stats.R')

# Data Prep ---------------------------------------------------------------

#Loading in the data
data_neg <- read_csv(here('./Data/Untargeted_Metabolomics/Clean/Metabolomoics_Neg_PQN_Metabolanalyst.csv'))
data_pos <- read_csv(here('./Data/Untargeted_Metabolomics/Clean/Metabolomics_Pos_withNIST_PQN_R.csv'))

#Load in metadata
meta_sprout <- read_csv(here('./Data/Metadata/sprout_meta.csv'))
meta_treatment <- read_csv(here('./Data/Metadata/meta_treatment.csv'))


## Negative Mode -----------------------------------------------------------

#Clean data 
clean_data_neg <- data_neg %>%
  #Load in PQN normalized data from metaboanalyst
  column_to_rownames("Feature") %>%
  #Transpose the matrix to make it how i want
  t() %>%
  #Conver to dataframe
  as.data.frame()%>%
  #Remove rownames so it's just a matrix of features and intensities
  rownames_to_column("sample") %>%
  #Breakup sample names to make metadata
  mutate(np = str_split(sample, '_')) %>%
  mutate(subject_id = map_chr(np, function(x) x[1])) %>%
  mutate(treatment = map_chr(np, function(x) x[2])) %>%
  mutate(time = map_chr(np, function(x) x[3])) %>%
  mutate(time = gsub('h', '', time)) %>%
  dplyr::select(-np) %>%
  #Filter to just human samples
  filter(str_detect(sample, 'BSS')) %>%
  #Agglomerate treatments
  mutate(veg = ifelse(treatment %in% c('AL', 'AU'), 'alf', 'broc')) %>%
  #Factorize variables and set levels
  modify_at('time', factor, levels = c(0,3,6,24,48,72)) %>%
  modify_at('treatment', factor, levels = c('AU', 'AL', 'BU', 'BL')) %>%
  modify_at('veg', factor, levels = c('alf', 'broc'))

#Clean up the data for lmer
lmer_data_neg <- clean_data_neg %>%
  #Make the data tidy
  pivot_longer(cols = c(ends_with('m/z'), ends_with('n')), names_to = 'feature', values_to = 'abd') %>%
  #Join in metadata
  left_join(meta_treatment) %>%
  left_join(meta_sprout) %>%
  #Factorize
  modify_at('treatment', factor, levels = c('AU', 'AL', 'BU', 'BL')) %>%
  group_by(feature) %>%
  #Nest
  nest()


## Positive Mode -----------------------------------------------------------

clean_data_pos <- data_pos %>%
  column_to_rownames("Feature") %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("sample") %>%
  mutate(np = str_split(sample, '_')) %>%
  mutate(subject_id = map_chr(np, function(x) x[1])) %>%
  mutate(treatment = map_chr(np, function(x) x[2])) %>%
  mutate(time = map_chr(np, function(x) x[3])) %>%
  mutate(time = gsub('h', '', time)) %>%
  dplyr::select(-np) %>%
  filter(str_detect(sample, 'BSS')) %>%
  mutate(veg = ifelse(treatment %in% c('AL', 'AU'), 'alf', 'broc')) %>%
  modify_at('time', factor, levels = c(0,3,6,24,48,72)) %>%
  modify_at('treatment', factor, levels = c('AU', 'AL', 'BU', 'BL')) %>%
  modify_at('veg', factor, levels = c('alf', 'broc'))

lmer_data_pos <- clean_data_pos %>%
  pivot_longer(cols = c(ends_with('m/z'), ends_with('n')), names_to = 'feature', values_to = 'abd') %>%
  left_join(meta_treatment) %>%
  left_join(meta_sprout) %>%
  modify_at('treatment', factor, levels = c('AU', 'AL', 'BU', 'BL')) %>%
  group_by(feature) %>%
  nest()

# GLMM --------------------------------------------------------------------

## Setup Contrasts ---------------------------------------------------------

BvA_0h <- c(0,1,0,0,0,0,0,0,0,0,0,0,0)
BvA_3h <- c(0,1,0,0,0,0,0,0,1,0,0,0,0)
BvA_6h <- c(0,1,0,0,0,0,0,0,0,1,0,0,0)
BvA_24h <- c(0,1,0,0,0,0,0,0,0,0,1,0,0)
BvA_48h <- c(0,1,0,0,0,0,0,0,0,0,0,1,0)
BvA_72h <- c(0,1,0,0,0,0,0,0,0,0,0,0,1)


## Negative Mode -----------------------------------------------------------

#Run the model
lmer_mod_neg <- lmer_data_neg %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ veg*time + grams_fed + (1|subject_id), data = x)))

#Run contrasts
lmer_contra_neg <- lmer_mod_neg %>%
  mutate(h0 = map_chr(model, function(x) lmerTest::contest(x, BvA_0h)$`Pr(>F)`)) %>%
  mutate(h3 = map_chr(model, function(x) lmerTest::contest(x, BvA_3h)$`Pr(>F)`)) %>%
  mutate(h6 = map_chr(model, function(x) lmerTest::contest(x, BvA_6h)$`Pr(>F)`)) %>%
  mutate(h24 = map_chr(model, function(x) lmerTest::contest(x, BvA_24h)$`Pr(>F)`)) %>%
  mutate(h48 = map_chr(model, function(x) lmerTest::contest(x, BvA_48h)$`Pr(>F)`)) %>%
  mutate(h72 = map_chr(model, function(x) lmerTest::contest(x, BvA_72h)$`Pr(>F)`)) %>%
  ungroup() 

#Adjust p-values
checkmat_neg <- lmer_contra_neg %>%
  mutate(h0_adj = p.adjust(h0, 'BH')) %>%
  mutate(h3_adj = p.adjust(h3, 'BH')) %>%
  mutate(h6_adj = p.adjust(h6, 'BH')) %>%
  mutate(h24_adj = p.adjust(h24, 'BH')) %>%
  mutate(h48_adj = p.adjust(h48, 'BH')) %>%
  mutate(h72_adj = p.adjust(h72, 'BH')) %>%
  column_to_rownames('feature') 

#Check overall
apply(checkmat_neg %>% dplyr::select(contains('adj')), 2, function(x) sum(x <= 0.05))


## Positive Mode -----------------------------------------------------------

#Run the model
lmer_mod_pos <- lmer_data_pos %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ veg*time + grams_fed + (1|subject_id), data = x)))

#Run contrasts
lmer_contra_pos <- lmer_mod_pos %>%
  mutate(h0 = map_chr(model, function(x) lmerTest::contest(x, BvA_0h)$`Pr(>F)`)) %>%
  mutate(h3 = map_chr(model, function(x) lmerTest::contest(x, BvA_3h)$`Pr(>F)`)) %>%
  mutate(h6 = map_chr(model, function(x) lmerTest::contest(x, BvA_6h)$`Pr(>F)`)) %>%
  mutate(h24 = map_chr(model, function(x) lmerTest::contest(x, BvA_24h)$`Pr(>F)`)) %>%
  mutate(h48 = map_chr(model, function(x) lmerTest::contest(x, BvA_48h)$`Pr(>F)`)) %>%
  mutate(h72 = map_chr(model, function(x) lmerTest::contest(x, BvA_72h)$`Pr(>F)`)) %>%
  ungroup() 

#Adjust p-values
checkmat_pos <- lmer_contra_pos %>%
  mutate(h0_adj = p.adjust(h0, 'BH')) %>%
  mutate(h3_adj = p.adjust(h3, 'BH')) %>%
  mutate(h6_adj = p.adjust(h6, 'BH')) %>%
  mutate(h24_adj = p.adjust(h24, 'BH')) %>%
  mutate(h48_adj = p.adjust(h48, 'BH')) %>%
  mutate(h72_adj = p.adjust(h72, 'BH')) %>%
  column_to_rownames('feature') 

#Check overall
apply(checkmat_pos %>% dplyr::select(contains('adj')), 2, function(x) sum(x <= 0.05))


# Generate Stats Tables ---------------------------------------------------

## Helper Functions --------------------------------------------------------

ReplaceZeros <- function(x){
  #Drop all constant rows (metabolites) from the dataframe:
  xnz <- as.data.frame(x[which(apply(x, 1, sd, na.rm = TRUE) != 0), ])
  rn <- rownames(xnz)
  #Find the minimum value of each column (sample) and divide by 5
  minval <- apply(xnz, 2, function(x) min(abs(x[x != 0]))/5)
  #Replace the 0s with the minimum value:
  rdf <- map2_df(xnz, minval, function(w,y){
   z <- w
   z[z == 0] <- y
   return(z)
  })
  rdf <- as.matrix(rdf)
  rownames(rdf) <- rn
  return(rdf)
}

ProbNorm <- function(x, ref.smpl){
  x/median(as.numeric(x/ref.smpl), na.rm=T)
}

PQN <- function(data, ref){
  out <- ReplaceZeros(data)
  ref.smpl <- out[ ,ref, drop=FALSE];
  data <- apply(out, 2, ProbNorm, ref.smpl);
  return(data)
}

calc_FC <- function(x, condition1, condition2, meta, value){
  mu1 <- mean(x[x[[meta]] == condition1, colnames(x) == value][[1]])
  mu2 <- mean(x[x[[meta]] == condition2, colnames(x) == value][[1]])
  fc <- mu2/mu1
  return(fc)
}


## Negative Mode -----------------------------------------------------------

### Data Prep ---------------------------------------------------------------

data_pqn_neg <- read_csv(here('./Data/Untargeted_Metabolomics/Raw/NegQCFiltered_PQN_RFormat.csv')) 

data_working_neg <- data_pqn_neg %>%
  column_to_rownames('Compound') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  filter(str_detect(sample, 'BSS')) %>%
  mutate(nf = str_split(sample, '_')) %>%
  mutate(treatment = map_chr(nf, function(x) x[2])) %>%
  mutate(veg = ifelse(treatment %in% c('BU', 'BL'), 'broc', 'alf')) %>%
  mutate(time = map_chr(nf, function(x) x[3])) %>%
  dplyr::select(-nf)


### 3h ----------------------------------------------------------------------

sig_3h_neg <- checkmat_neg %>%
  filter(h3_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_3h_neg <- data_working_neg %>%
  filter(time == '3h') %>%
  #dplyr::select(veg, all_of(sig_3h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_neg %>% dplyr::select(h3, h3_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h3_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h3_adj < 0.05, 'sig', 'non-sig')) 

table_3h_neg %>%
  dplyr::select(FT, FC, log2_FC, h3, h3_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h3', 'padj' = 'h3_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/3h_neg.csv'))

table_3h_neg %>%
  dplyr::select(FT, FC, log2_FC, h3, h3_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h3', 'padj' = 'h3_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/3h_neg_all.csv'))

### 6h ----------------------------------------------------------------------

sig_6h_neg <- checkmat_neg %>%
  filter(h6_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_6h_neg <- data_working_neg %>%
  filter(time == '6h') %>%
  #dplyr::select(veg, all_of(sig_6h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_neg %>% dplyr::select(h6, h6_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h6_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h6_adj < 0.05, 'sig', 'non-sig')) 

table_6h_neg %>%
  dplyr::select(FT, FC, log2_FC, h6, h6_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h6', 'padj' = 'h6_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/6h_neg.csv'))

table_6h_neg %>%
  dplyr::select(FT, FC, log2_FC, h6, h6_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h6', 'padj' = 'h6_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/6h_neg_all.csv'))

### 24h ----------------------------------------------------------------------

sig_24h_neg <- checkmat_neg %>%
  filter(h24_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_24h_neg <- data_working_neg %>%
  filter(time == '24h') %>%
  #dplyr::select(veg, all_of(sig_24h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_neg %>% dplyr::select(h24, h24_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h24_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h24_adj < 0.05, 'sig', 'non-sig')) 

table_24h_neg %>%
  dplyr::select(FT, FC, log2_FC, h24, h24_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h24', 'padj' = 'h24_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/24h_neg.csv'))

table_24h_neg %>%
  dplyr::select(FT, FC, log2_FC, h24, h24_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h24', 'padj' = 'h24_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/24h_neg_all.csv'))


### 48h ----------------------------------------------------------------------

sig_48h_neg <- checkmat_neg %>%
  filter(h48_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_48h_neg <- data_working_neg %>%
  filter(time == '48h') %>%
  #dplyr::select(veg, all_of(sig_48h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_neg %>% dplyr::select(h48, h48_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h48_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h48_adj < 0.05, 'sig', 'non-sig')) 

table_48h_neg %>%
  dplyr::select(FT, FC, log2_FC, h48, h48_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h48', 'padj' = 'h48_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/48h_neg.csv'))

table_48h_neg %>%
  dplyr::select(FT, FC, log2_FC, h48, h48_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h48', 'padj' = 'h48_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/48h_neg_all.csv'))

### 72h ----------------------------------------------------------------------

sig_72h_neg <- checkmat_neg %>%
  filter(h72_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_72h_neg <- data_working_neg %>%
  filter(time == '72h') %>%
  #dplyr::select(veg, all_of(sig_72h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_neg %>% dplyr::select(h72, h72_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h72_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h72_adj < 0.05, 'sig', 'non-sig')) 

table_72h_neg %>%
  dplyr::select(FT, FC, log2_FC, h72, h72_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h72', 'padj' = 'h72_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/72h_neg.csv'))

table_72h_neg %>%
  dplyr::select(FT, FC, log2_FC, h72, h72_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h72', 'padj' = 'h72_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/72h_neg_all.csv'))


## Positive Mode -----------------------------------------------------------

dataraw_pos <- read_csv(here('./Data/Untargeted_Metabolomics/Raw/PosQCFiltered_withNIST_RFormat.csv')) %>%
  column_to_rownames('Compound') 
data_pqn_pos <- PQN(dataraw_pos, ref = 'QC_urine_20') 

data_working_pos <- data_pqn_pos %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  filter(str_detect(sample, 'BSS')) %>%
  mutate(nf = str_split(sample, '_')) %>%
  mutate(treatment = map_chr(nf, function(x) x[2])) %>%
  mutate(veg = ifelse(treatment %in% c('BU', 'BL'), 'broc', 'alf')) %>%
  mutate(time = map_chr(nf, function(x) x[3])) %>%
  dplyr::select(-nf)

### 3h ----------------------------------------------------------------------

sig_3h_pos <- checkmat_pos %>%
  filter(h3_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_3h_pos <- data_working_pos %>%
  filter(time == '3h') %>%
  #dplyr::select(veg, all_of(sig_3h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_pos %>% dplyr::select(h3, h3_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h3_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h3_adj < 0.05, 'sig', 'non-sig')) 

table_3h_pos %>%
  dplyr::select(FT, FC, log2_FC, h3, h3_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h3', 'padj' = 'h3_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/3h_pos.csv'))

table_3h_pos %>%
  dplyr::select(FT, FC, log2_FC, h3, h3_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h3', 'padj' = 'h3_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/3h_pos_all.csv'))

### 6h ----------------------------------------------------------------------

sig_6h_pos <- checkmat_pos %>%
  filter(h6_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_6h_pos <- data_working_pos %>%
  filter(time == '6h') %>%
  #dplyr::select(veg, all_of(sig_6h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_pos %>% dplyr::select(h6, h6_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h6_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h6_adj < 0.05, 'sig', 'non-sig')) 

table_6h_pos %>%
  dplyr::select(FT, FC, log2_FC, h6, h6_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h6', 'padj' = 'h6_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/6h_pos.csv'))

table_6h_pos %>%
  dplyr::select(FT, FC, log2_FC, h6, h6_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h6', 'padj' = 'h6_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/6h_pos_all.csv'))

### 24h ----------------------------------------------------------------------

sig_24h_pos <- checkmat_pos %>%
  filter(h24_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_24h_pos <- data_working_pos %>%
  filter(time == '24h') %>%
  #dplyr::select(veg, all_of(sig_24h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_pos %>% dplyr::select(h24, h24_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h24_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h24_adj < 0.05, 'sig', 'non-sig')) 

table_24h_pos %>%
  dplyr::select(FT, FC, log2_FC, h24, h24_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h24', 'padj' = 'h24_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/24h_pos.csv'))

table_24h_pos %>%
  dplyr::select(FT, FC, log2_FC, h24, h24_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h24', 'padj' = 'h24_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/24h_pos_all.csv'))


### 48h ----------------------------------------------------------------------

sig_48h_pos <- checkmat_pos %>%
  filter(h48_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_48h_pos <- data_working_pos %>%
  filter(time == '48h') %>%
  #dplyr::select(veg, all_of(sig_48h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_pos %>% dplyr::select(h48, h48_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h48_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h48_adj < 0.05, 'sig', 'non-sig')) 

table_48h_pos %>%
  dplyr::select(FT, FC, log2_FC, h48, h48_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h48', 'padj' = 'h48_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/48h_pos.csv'))

table_48h_pos %>%
  dplyr::select(FT, FC, log2_FC, h48, h48_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h48', 'padj' = 'h48_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/48h_pos_all.csv'))

### 72h ----------------------------------------------------------------------

sig_72h_pos <- checkmat_pos %>%
  filter(h72_adj <= 0.05) %>%
  rownames_to_column('feature') %>%
  pull(feature)

table_72h_pos <- data_working_pos %>%
  filter(time == '72h') %>%
  #dplyr::select(veg, all_of(sig_72h_pos)) %>%
  pivot_longer(ends_with(c('m/z', 'n')), names_to = 'FT', values_to = 'int') %>%
  group_by(FT) %>%
  nest() %>%
  #If alfalfa is higher there will be a negative FC, if broccoli is higher there will be a positive FC
  mutate(FC = map_dbl(data, function(x) calc_FC(x, condition1 = 'alf', condition2 = 'broc', meta = 'veg', value = 'int'))) %>%
  mutate(log2_FC = log2(FC)) %>%
  left_join(., (checkmat_pos %>% dplyr::select(h72, h72_adj) %>% rownames_to_column('FT'))) %>%
  mutate(log_p = -log10(h72_adj + 1E-6)) %>%
  mutate(highest_mean = ifelse(sign(log2_FC) == -1, 'alf', 'broc')) %>%
  mutate(sig = ifelse(h72_adj < 0.05, 'sig', 'non-sig')) 

table_72h_pos %>%
  dplyr::select(FT, FC, log2_FC, h72, h72_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h72', 'padj' = 'h72_adj') %>%
  filter(padj <= 0.05) %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/Sig/72h_pos.csv'))

table_72h_pos %>%
  dplyr::select(FT, FC, log2_FC, h72, h72_adj, highest_mean) %>%
  rename('Compound' = 'FT', 'pvalue' = 'h72', 'padj' = 'h72_adj') %>%
  write_csv(here('./human_analysis/Untargeted_MultiOmics/stats_tables/All/72h_pos_all.csv'))


