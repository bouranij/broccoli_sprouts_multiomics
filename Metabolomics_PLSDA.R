#Code to run PLS-DAs and PCAs in R 


# Environment -------------------------------------------------------------

library(tidyverse)
library(here)
library(magrittr)
setwd('~/Documents/Projects/Short_Term_Broccoli/')
here::i_am('./human_analysis/Untargeted_MultiOmics/Biomarker_Stats.R')

# Data Prep ---------------------------------------------------------------

#Loading in the data
data <- read_csv(here('./Data/Untargeted_Metabolomics/Clean/Metabolomoics_Neg_PQN_Metabolanalyst.csv'))
data_pos_nist <- read_csv(here('./Data/Untargeted_Metabolomics/Clean/Metabolomics_Pos_withNIST_PQN_R.csv'))

#Load in metadata
meta_sprout <- read_csv(here('./Data/Metadata/sprout_meta.csv'))
meta_treatment <- read_csv(here('./Data/Metadata/meta_treatment.csv'))


## Negative Mode -----------------------------------------------------------

#Clean data 
clean_data <- data %>%
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
lmer_data <- clean_data %>%
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

clean_data_pos_nist <- data_pos_nist %>%
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

lmer_data_pos_nist <- clean_data_pos_nist %>%
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
lmer_mod2 <- lmer_data %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ veg*time + grams_fed + (1|subject_id), data = x)))

#Run contrasts
lmer_contra2 <- lmer_mod2 %>%
  mutate(h0 = map_chr(model, function(x) lmerTest::contest(x, BvA_0h)$`Pr(>F)`)) %>%
  mutate(h3 = map_chr(model, function(x) lmerTest::contest(x, BvA_3h)$`Pr(>F)`)) %>%
  mutate(h6 = map_chr(model, function(x) lmerTest::contest(x, BvA_6h)$`Pr(>F)`)) %>%
  mutate(h24 = map_chr(model, function(x) lmerTest::contest(x, BvA_24h)$`Pr(>F)`)) %>%
  mutate(h48 = map_chr(model, function(x) lmerTest::contest(x, BvA_48h)$`Pr(>F)`)) %>%
  mutate(h72 = map_chr(model, function(x) lmerTest::contest(x, BvA_72h)$`Pr(>F)`)) %>%
  ungroup() 

#Adjust p-values
checkmat2 <- lmer_contra2 %>%
  mutate(h0_adj = p.adjust(h0, 'BH')) %>%
  mutate(h3_adj = p.adjust(h3, 'BH')) %>%
  mutate(h6_adj = p.adjust(h6, 'BH')) %>%
  mutate(h24_adj = p.adjust(h24, 'BH')) %>%
  mutate(h48_adj = p.adjust(h48, 'BH')) %>%
  mutate(h72_adj = p.adjust(h72, 'BH')) %>%
  column_to_rownames('feature') %>%
  dplyr::select(contains('adj'))

#Check overall
apply(checkmat2, 2, function(x) sum(x <= 0.05))


## Positive Mode -----------------------------------------------------------

#Run the model
lmer_mod2_pos_nist <- lmer_data_pos_nist %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ veg*time + grams_fed + (1|subject_id), data = x)))

#Run contrasts
lmer_contra2_pos_nist <- lmer_mod2_pos_nist %>%
  mutate(h0 = map_chr(model, function(x) lmerTest::contest(x, BvA_0h)$`Pr(>F)`)) %>%
  mutate(h3 = map_chr(model, function(x) lmerTest::contest(x, BvA_3h)$`Pr(>F)`)) %>%
  mutate(h6 = map_chr(model, function(x) lmerTest::contest(x, BvA_6h)$`Pr(>F)`)) %>%
  mutate(h24 = map_chr(model, function(x) lmerTest::contest(x, BvA_24h)$`Pr(>F)`)) %>%
  mutate(h48 = map_chr(model, function(x) lmerTest::contest(x, BvA_48h)$`Pr(>F)`)) %>%
  mutate(h72 = map_chr(model, function(x) lmerTest::contest(x, BvA_72h)$`Pr(>F)`)) %>%
  ungroup() 

#Adjust p-values
checkmat2_pos_nist <- lmer_contra2_pos_nist %>%
  mutate(h0_adj = p.adjust(h0, 'BH')) %>%
  mutate(h3_adj = p.adjust(h3, 'BH')) %>%
  mutate(h6_adj = p.adjust(h6, 'BH')) %>%
  mutate(h24_adj = p.adjust(h24, 'BH')) %>%
  mutate(h48_adj = p.adjust(h48, 'BH')) %>%
  mutate(h72_adj = p.adjust(h72, 'BH')) %>%
  column_to_rownames('feature') %>%
  dplyr::select(contains('adj'))

#Check overall
apply(checkmat2_pos_nist, 2, function(x) sum(x <= 0.05))

# PLS-DA Data Prep --------------------------------------------------------

#Alter feature names so we can tell Pos and Neg mode apart
data_neg <- data
feature_map_neg <- data_neg %>%
  dplyr::select(Feature) %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_neg', Feature))) 

feature_map_pos <- data_pos_nist %>%
  dplyr::select(Feature) %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_pos', Feature))) 

data_clean_neg <- data_neg %>%
  filter(!Feature %in% ft_remove) %>%
  right_join(feature_map_neg, .) %>%
  dplyr::select(-Feature) %>%
  column_to_rownames("fixed_name") %>%
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

data_clean_pos <- data_pos_nist %>%
  filter(!Feature %in% ft_remove) %>%
  right_join(feature_map_pos, .) %>%
  dplyr::select(-Feature) %>%
  column_to_rownames("fixed_name") %>%
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

#Extract only siginificant features:
pos_sig <- checkmat2_pos_nist %>%
  rownames_to_column('feature') %>%
  pivot_longer(cols = ends_with('adj'), names_to = 'test', values_to = 'pval') %>%
  group_by(feature) %>%
  filter(any(pval <= 0.05)) %>%
  pivot_wider(names_from = 'test', values_from = 'pval') %>%
  pull(feature) %>%
  paste0('FT_', .) %>%
  gsub('(m/z|n)', '_pos', .)

pos_block <- data_clean_pos %>%
  dplyr::select(subject_id, treatment, time, veg, all_of(pos_sig)) 

neg_sig <- checkmat2 %>%
  rownames_to_column('feature') %>%
  pivot_longer(cols = ends_with('adj'), names_to = 'test', values_to = 'pval') %>%
  group_by(feature) %>%
  filter(any(pval <= 0.05)) %>%
  pivot_wider(names_from = 'test', values_from = 'pval') %>%
  pull(feature) %>%
  paste0('FT_', .) %>%
  gsub('(m/z|n)', '_neg', .)

neg_block <- data_clean_neg %>%
  dplyr::select(subject_id, treatment, time, veg, all_of(neg_sig)) 

#Combine into one big block
full_block <- left_join(pos_block, neg_block)

#Extract out individual timepoints for the PLS-DA
metab_3h <- full_block %>%
  filter(time == 3) %>%
  dplyr::select(starts_with('FT')) %>%
  as.matrix()

metab_3h <- full_block %>%
  filter(time == 3) %>%
  dplyr::select(starts_with('FT')) %>%
  as.matrix()

metab_6h <- full_block %>%
  filter(time == 6) %>%
  dplyr::select(starts_with('FT')) %>%
  as.matrix()

metab_24h <- full_block %>%
  filter(time == 24) %>%
  dplyr::select(starts_with('FT')) %>%
  as.matrix()

metab_48h <- full_block %>%
  filter(time == 48) %>%
  dplyr::select(starts_with('FT')) %>%
  as.matrix()

metab_72h <- full_block %>%
  filter(time == 72) %>%
  dplyr::select(starts_with('FT')) %>%
  as.matrix()

#Y-block:
treatment_block <- full_block %>%
  filter(time == 3) %>%
  dplyr::pull(veg) 

# PLS-DA ------------------------------------------------------------------

library(mixOmics)
library(pROC)
#Set up aesthetics:
pca_key <- data.frame(treatment = treatment_block) %>%
  mutate(color = ifelse(treatment == 'alf', 'orchid4', 'darkgreen'),
         pch = ifelse(treatment == 'alf', 15, 16))

#Pull out class information:
classmat <- as.data.frame(treatment_block) %>%
  mutate(alf = ifelse(treatment_block == 'alf', 1, 0)) %>%
  mutate(broc = ifelse(treatment_block == 'broc', 1, 0)) %>%
  dplyr::select(-treatment_block)

## 3h ----------------------------------------------------------------------

#PCA - Run and Plot
metab_3h_pca <- pca(X = metab_3h)
plotIndiv(metab_3h_pca, col = unique(pca_key$color), pch = unique(pca_key$pch), ellipse = T, 
          group = pca_key$treatment, title = '3h PCA', legend = T)

#PLS-DA - Run and plot:
metab_3h_plsda <- plsda(X = metab_3h, Y = treatment_block)
#Cross validate using LOO-CV
metab_3h_perf <- perf(metab_3h_plsda, validation = 'loo', auc = T)

#Pull out LOO results
nc1_3h <- as.data.frame(metab_3h_perf$predict$comp1) %>%
  set_colnames(c('alf', 'broc'))
#Generate ROC Curve
roc3h <- roc(classmat$broc, nc1_3h$broc)
#Get the AUC
auc(roc3h)
#Plot the ROC curve
ggroc(roc3h, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc3h)))

#Cross validate using Mfold-CV
metab_3h_perf_m <- perf(metab_3h_plsda, validation = 'Mfold', folds = 3, auc = T, nrepeat = 10)
out_3h <- array_branch(metab_3h_perf_m$predict$comp1, 3) %>%
  purrr::map_dfc(., function(x) x[,2]) %>%
  apply(., 1, mean) %>%
  as.data.frame() %>%
  set_colnames('broc')
roc3h_m <- roc(classmat$broc, out_3h$broc)
#Get the AUC
auc(roc3h_m)
#Plot the ROC curve
ggroc(roc3h_m, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc3h_m)))

## 6h ----------------------------------------------------------------------

#PCA - Run and Plot
metab_6h_pca <- pca(X = metab_6h)
plotIndiv(metab_6h_pca, col = unique(pca_key$color), pch = unique(pca_key$pch), ellipse = T, 
          group = pca_key$treatment, title = '6h PCA', legend = T)

#PLS-DA - Run and plot:
metab_6h_plsda <- plsda(X = metab_6h, Y = treatment_block)
#Cross validate using LOO-CV
metab_6h_perf <- perf(metab_6h_plsda, validation = 'loo', auc = T)

#Pull out LOO results
nc1_6h <- as.data.frame(metab_6h_perf$predict$comp1) %>%
  set_colnames(c('alf', 'broc'))
#Generate ROC Curve
roc6h <- roc(classmat$broc, nc1_6h$broc)
#Get the AUC
auc(roc6h)
#Plot the ROC curve
ggroc(roc6h, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc6h)))

#Cross validate using Mfold-CV
metab_6h_perf_m <- perf(metab_6h_plsda, validation = 'Mfold', folds = 3, auc = T, nrepeat = 10)
out_6h <- array_branch(metab_6h_perf_m$predict$comp1, 3) %>%
  purrr::map_dfc(., function(x) x[,2]) %>%
  apply(., 1, mean) %>%
  as.data.frame() %>%
  set_colnames('broc')
roc6h_m <- roc(classmat$broc, out_6h$broc)
#Get the AUC
auc(roc6h_m)
#Plot the ROC curve
ggroc(roc6h_m, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc6h_m)))

## 24h ----------------------------------------------------------------------

#PCA - Run and Plot
metab_24h_pca <- pca(X = metab_24h)
plotIndiv(metab_24h_pca, col = unique(pca_key$color), pch = unique(pca_key$pch), ellipse = T, 
          group = pca_key$treatment, title = '24h PCA', legend = T)

#PLS-DA - Run and plot:
metab_24h_plsda <- plsda(X = metab_24h, Y = treatment_block)
#Cross validate using LOO-CV
metab_24h_perf <- perf(metab_24h_plsda, validation = 'loo', auc = T)

#Pull out LOO results
nc1_24h <- as.data.frame(metab_24h_perf$predict$comp1) %>%
  set_colnames(c('alf', 'broc'))
#Generate ROC Curve
roc24h <- roc(classmat$broc, nc1_24h$broc)
#Get the AUC
auc(roc24h)
#Plot the ROC curve
ggroc(roc24h, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc24h)))

#Cross validate using Mfold-CV
metab_24h_perf_m <- perf(metab_24h_plsda, validation = 'Mfold', folds = 3, auc = T, nrepeat = 10)
out_24h <- array_branch(metab_24h_perf_m$predict$comp1, 3) %>%
  purrr::map_dfc(., function(x) x[,2]) %>%
  apply(., 1, mean) %>%
  as.data.frame() %>%
  set_colnames('broc')
roc24h_m <- roc(classmat$broc, out_24h$broc)
#Get the AUC
auc(roc24h_m)
#Plot the ROC curve
ggroc(roc24h_m, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc24h_m)))

## 48h ----------------------------------------------------------------------

#PCA - Run and Plot
metab_48h_pca <- pca(X = metab_48h)
plotIndiv(metab_48h_pca, col = unique(pca_key$color), pch = unique(pca_key$pch), ellipse = T, 
          group = pca_key$treatment, title = '48h PCA', legend = T)

#PLS-DA - Run and plot:
metab_48h_plsda <- plsda(X = metab_48h, Y = treatment_block)
#Cross validate using LOO-CV
metab_48h_perf <- perf(metab_48h_plsda, validation = 'loo', auc = T)

#Pull out LOO results
nc1_48h <- as.data.frame(metab_48h_perf$predict$comp1) %>%
  set_colnames(c('alf', 'broc'))
#Generate ROC Curve
roc48h <- roc(classmat$broc, nc1_48h$broc)
#Get the AUC
auc(roc48h)
#Plot the ROC curve
ggroc(roc48h, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc48h)))

#Cross validate using Mfold-CV
metab_48h_perf_m <- perf(metab_48h_plsda, validation = 'Mfold', folds = 3, auc = T, nrepeat = 10)
out_48h <- array_branch(metab_48h_perf_m$predict$comp1, 3) %>%
  purrr::map_dfc(., function(x) x[,2]) %>%
  apply(., 1, mean) %>%
  as.data.frame() %>%
  set_colnames('broc')
roc48h_m <- roc(classmat$broc, out_48h$broc)
#Get the AUC
auc(roc48h_m)
#Plot the ROC curve
ggroc(roc48h_m, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc48h_m)))

## 72h ----------------------------------------------------------------------

#PCA - Run and Plot
metab_72h_pca <- pca(X = metab_72h)
plotIndiv(metab_72h_pca, col = unique(pca_key$color), pch = unique(pca_key$pch), ellipse = T, 
          group = pca_key$treatment, title = '72h PCA', legend = T)

#PLS-DA - Run and plot:
metab_72h_plsda <- plsda(X = metab_72h, Y = treatment_block)
#Cross validate using LOO-CV
metab_72h_perf <- perf(metab_72h_plsda, validation = 'loo', auc = T)

#Pull out LOO results
nc1_72h <- as.data.frame(metab_72h_perf$predict$comp1) %>%
  set_colnames(c('alf', 'broc'))
#Generate ROC Curve
roc72h <- roc(classmat$broc, nc1_72h$broc)
#Get the AUC
auc(roc72h)
#Plot the ROC curve
ggroc(roc72h, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc72h)))

#Cross validate using Mfold-CV
metab_72h_perf_m <- perf(metab_72h_plsda, validation = 'Mfold', folds = 3, auc = T, nrepeat = 10)
out_72h <- array_branch(metab_72h_perf_m$predict$comp1, 3) %>%
  purrr::map_dfc(., function(x) x[,2]) %>%
  apply(., 1, mean) %>%
  as.data.frame() %>%
  set_colnames('broc')
roc72h_m <- roc(classmat$broc, out_72h$broc)
#Get the AUC
auc(roc72h_m)
#Plot the ROC curve
ggroc(roc72h_m, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', auc(roc72h_m)))


# Generate Result Tables --------------------------------------------------

## 3h ----------------------------------------------------------------------

#Pull VIP Scores
vip_3h <- as.data.frame(vip(metab_3h_plsda)) %>%
  rownames_to_column('feature') %>%
  arrange(desc(comp1))

#Pull group contribution information
loadings_3h <- plotLoadings(metab_3h_plsda, plot = F, contrib = 'max') %>%
  rownames_to_column('feature') %>%
  dplyr::select(feature, GroupContrib) %>%
  rename('Highest_Group' = 'GroupContrib')

#Combine
final_var_3h <- left_join(vip_3h, loadings_3h)

#Write to CSV
write_csv(final_var_3h, file = './PLSDA_Output/3h_metab.csv')

## 6h ----------------------------------------------------------------------

#Pull VIP Scores
vip_6h <- as.data.frame(vip(metab_6h_plsda)) %>%
  rownames_to_column('feature') %>%
  arrange(desc(comp1))

#Pull group contribution information
loadings_6h <- plotLoadings(metab_6h_plsda, plot = F, contrib = 'max') %>%
  rownames_to_column('feature') %>%
  dplyr::select(feature, GroupContrib) %>%
  rename('Highest_Group' = 'GroupContrib')

#Combine
final_var_6h <- left_join(vip_6h, loadings_6h)

#Write to CSV
write_csv(final_var_6h, file = './PLSDA_Output/6h_metab.csv')

## 24h ----------------------------------------------------------------------

#Pull VIP Scores
vip_24h <- as.data.frame(vip(metab_24h_plsda)) %>%
  rownames_to_column('feature') %>%
  arrange(desc(comp1))

#Pull group contribution information
loadings_24h <- plotLoadings(metab_24h_plsda, plot = F, contrib = 'max') %>%
  rownames_to_column('feature') %>%
  dplyr::select(feature, GroupContrib) %>%
  rename('Highest_Group' = 'GroupContrib')

#Combine
final_var_24h <- left_join(vip_24h, loadings_24h)

#Write to CSV
write_csv(final_var_24h, file = './PLSDA_Output/24h_metab.csv')


## 48h ----------------------------------------------------------------------

#Pull VIP Scores
vip_48h <- as.data.frame(vip(metab_48h_plsda)) %>%
  rownames_to_column('feature') %>%
  arrange(desc(comp1))

#Pull group contribution information
loadings_48h <- plotLoadings(metab_48h_plsda, plot = F, contrib = 'max') %>%
  rownames_to_column('feature') %>%
  dplyr::select(feature, GroupContrib) %>%
  rename('Highest_Group' = 'GroupContrib')

#Combine
final_var_48h <- left_join(vip_48h, loadings_48h)

#Write to CSV
write_csv(final_var_48h, file = './PLSDA_Output/48h_metab.csv')

## 72h ----------------------------------------------------------------------

#Pull VIP Scores
vip_72h <- as.data.frame(vip(metab_72h_plsda)) %>%
  rownames_to_column('feature') %>%
  arrange(desc(comp1))

#Pull group contribution information
loadings_72h <- plotLoadings(metab_72h_plsda, plot = F, contrib = 'max') %>%
  rownames_to_column('feature') %>%
  dplyr::select(feature, GroupContrib) %>%
  rename('Highest_Group' = 'GroupContrib')

#Combine
final_var_72h <- left_join(vip_72h, loadings_72h)

#Write to CSV
write_csv(final_var_72h, file = './PLSDA_Output/72h_metab.csv')


