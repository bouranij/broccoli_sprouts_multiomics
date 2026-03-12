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

#At 0h, L vs UL
h0_ALvAU <- c(0,1,rep(0, 23))
h0_BLvBU <- c(0,0,-1,1,rep(0, 21))

#At 3h, L vs UL
h3_ALvAU <- c(0,1,0,0,0,0,0,0,0,0,1,rep(0,14))
h3_BLvBU <- c(0,0,-1,1,0,0,0,0,0,0,-1,1,rep(0,13))

#At 6h, L vs UL
h6_ALvAU <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,1,rep(0,11))
h6_BLvBU <- c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0,-1,1,rep(0,9))

#At 24h, L vs UL
h24_ALvAU <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,rep(0,8))
h24_BLvBU <- c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,rep(0,6))

#At 48h, L vs UL
h48_ALvAU <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,rep(0,5))
h48_BLvBU <- c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,rep(0,3))

#At 72h, L vs UL
h72_ALvAU <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
h72_BLvBU <- c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1)

## Negative Mode -----------------------------------------------------------

#Run the model for ever feature
lmer_label_neg <- lmer_data_neg %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ treatment*time + grams_fed + (1|subject_id), data = x)))

label_contra_neg <- lmer_label_neg %>%
  mutate(ALvsAU_0h = map_chr(model, function(x) lmerTest::contest(x, h0_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_0h = map_chr(model, function(x) lmerTest::contest(x, h0_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_3h = map_chr(model, function(x) lmerTest::contest(x, h3_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_3h = map_chr(model, function(x) lmerTest::contest(x, h3_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_6h = map_chr(model, function(x) lmerTest::contest(x, h6_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_6h = map_chr(model, function(x) lmerTest::contest(x, h6_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_24h = map_chr(model, function(x) lmerTest::contest(x, h24_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_24h = map_chr(model, function(x) lmerTest::contest(x, h24_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_48h = map_chr(model, function(x) lmerTest::contest(x, h48_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_48h = map_chr(model, function(x) lmerTest::contest(x, h48_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_72h = map_chr(model, function(x) lmerTest::contest(x, h72_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_72h = map_chr(model, function(x) lmerTest::contest(x, h72_BLvBU)$`Pr(>F)`)) %>%
  ungroup() 

#Adjust P-Values
checkmat_label_neg <- label_contra_neg %>%
  mutate(A0h_adj = p.adjust(ALvsAU_0h, 'BH')) %>%
  mutate(B0h_adj = p.adjust(BLvsBU_0h, 'BH')) %>%
  mutate(A3h_adj = p.adjust(ALvsAU_3h, 'BH')) %>%
  mutate(B3h_adj = p.adjust(BLvsBU_3h, 'BH')) %>%
  mutate(A6h_adj = p.adjust(ALvsAU_6h, 'BH')) %>%
  mutate(B6h_adj = p.adjust(BLvsBU_6h, 'BH')) %>%
  mutate(A24h_adj = p.adjust(ALvsAU_24h, 'BH')) %>%
  mutate(B24h_adj = p.adjust(BLvsBU_24h, 'BH')) %>%
  mutate(A48h_adj = p.adjust(ALvsAU_48h, 'BH')) %>%
  mutate(B48h_adj = p.adjust(BLvsBU_48h, 'BH')) %>%
  mutate(A72h_adj = p.adjust(ALvsAU_72h, 'BH')) %>%
  mutate(B72h_adj = p.adjust(BLvsBU_72h, 'BH')) %>%
  column_to_rownames('feature') %>%
  dplyr::select(contains('adj'))
#How many sig differences? Cool - we can combine L and UL
apply(checkmat_label_neg, 2, function(x) sum(x <= 0.05))


# Positive Mode -----------------------------------------------------------


#Run the model for ever feature
lmer_label_pos <- lmer_data_pos %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ treatment*time + grams_fed + (1|subject_id), data = x)))

label_contra_pos <- lmer_label_pos %>%
  mutate(ALvsAU_0h = map_chr(model, function(x) lmerTest::contest(x, h0_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_0h = map_chr(model, function(x) lmerTest::contest(x, h0_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_3h = map_chr(model, function(x) lmerTest::contest(x, h3_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_3h = map_chr(model, function(x) lmerTest::contest(x, h3_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_6h = map_chr(model, function(x) lmerTest::contest(x, h6_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_6h = map_chr(model, function(x) lmerTest::contest(x, h6_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_24h = map_chr(model, function(x) lmerTest::contest(x, h24_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_24h = map_chr(model, function(x) lmerTest::contest(x, h24_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_48h = map_chr(model, function(x) lmerTest::contest(x, h48_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_48h = map_chr(model, function(x) lmerTest::contest(x, h48_BLvBU)$`Pr(>F)`)) %>%
  mutate(ALvsAU_72h = map_chr(model, function(x) lmerTest::contest(x, h72_ALvAU)$`Pr(>F)`)) %>%
  mutate(BLvsBU_72h = map_chr(model, function(x) lmerTest::contest(x, h72_BLvBU)$`Pr(>F)`)) %>%
  ungroup() 

#Adjust P-Values
checkmat_label_pos <- label_contra_pos %>%
  mutate(A0h_adj = p.adjust(ALvsAU_0h, 'BH')) %>%
  mutate(B0h_adj = p.adjust(BLvsBU_0h, 'BH')) %>%
  mutate(A3h_adj = p.adjust(ALvsAU_3h, 'BH')) %>%
  mutate(B3h_adj = p.adjust(BLvsBU_3h, 'BH')) %>%
  mutate(A6h_adj = p.adjust(ALvsAU_6h, 'BH')) %>%
  mutate(B6h_adj = p.adjust(BLvsBU_6h, 'BH')) %>%
  mutate(A24h_adj = p.adjust(ALvsAU_24h, 'BH')) %>%
  mutate(B24h_adj = p.adjust(BLvsBU_24h, 'BH')) %>%
  mutate(A48h_adj = p.adjust(ALvsAU_48h, 'BH')) %>%
  mutate(B48h_adj = p.adjust(BLvsBU_48h, 'BH')) %>%
  mutate(A72h_adj = p.adjust(ALvsAU_72h, 'BH')) %>%
  mutate(B72h_adj = p.adjust(BLvsBU_72h, 'BH')) %>%
  column_to_rownames('feature') %>%
  dplyr::select(contains('adj'))

#How many sig differences? Cool - we can combine L and UL
apply(checkmat_label_pos, 2, function(x) sum(x <= 0.05))




