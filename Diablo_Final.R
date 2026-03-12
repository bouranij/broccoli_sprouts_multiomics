# Environment -------------------------------------------------------------
library(tidyverse)
library(here)
library(magrittr)
setwd('~/Documents/Projects/PhD/Short_Term_Broccoli/')
here::i_am('./human_analysis/Untargeted_MultiOmics/Diablo_Final.R')


# Metabolomcis Data Prep --------------------------------------------------

#Loading in the data
data_neg <- read_csv(here('./Data/Untargeted_Metabolomics/Clean/Metabolomoics_Neg_PQN_Metabolanalyst.csv'))
data_pos <- read_csv(here('./Data/Untargeted_Metabolomics/Clean/Metabolomics_Pos_withNIST_PQN_R.csv'))

#Load in metadata
meta_sprout <- read_csv(here('./Data/Metadata/sprout_meta.csv'))
meta_treatment <- read_csv(here('./Data/Metadata/meta_treatment.csv'))



## Data Cleaning -----------------------------------------------------------

### Negative Mode -----------------------------------------------------------

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

### Positive Mode -----------------------------------------------------------

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

## GLMM --------------------------------------------------------------------

### Setup Contrasts ---------------------------------------------------------

BvA_0h <- c(0,1,0,0,0,0,0,0,0,0,0,0,0)
BvA_3h <- c(0,1,0,0,0,0,0,0,1,0,0,0,0)
BvA_6h <- c(0,1,0,0,0,0,0,0,0,1,0,0,0)
BvA_24h <- c(0,1,0,0,0,0,0,0,0,0,1,0,0)
BvA_48h <- c(0,1,0,0,0,0,0,0,0,0,0,1,0)
BvA_72h <- c(0,1,0,0,0,0,0,0,0,0,0,0,1)


### Negative Mode -----------------------------------------------------------

#Run the model
lmer_mod_neg <- lmer_data_neg %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ veg*time + grams_fed + (1|subject_id), data = x)))

#Run contrasts
lmer_contra_neg <- lmer_mod_neg %>%
  mutate(h0 = map_dbl(model, function(x) lmerTest::contest(x, BvA_0h)$`Pr(>F)`)) %>%
  mutate(h3 = map_dbl(model, function(x) lmerTest::contest(x, BvA_3h)$`Pr(>F)`)) %>%
  mutate(h6 = map_dbl(model, function(x) lmerTest::contest(x, BvA_6h)$`Pr(>F)`)) %>%
  mutate(h24 = map_dbl(model, function(x) lmerTest::contest(x, BvA_24h)$`Pr(>F)`)) %>%
  mutate(h48 = map_dbl(model, function(x) lmerTest::contest(x, BvA_48h)$`Pr(>F)`)) %>%
  mutate(h72 = map_dbl(model, function(x) lmerTest::contest(x, BvA_72h)$`Pr(>F)`)) %>%
  ungroup() 


#Adjust p-values
checkmat_neg <- lmer_contra_neg %>%
  mutate(h0_adj = p.adjust(h0, 'BH')) %>%
  mutate(h3_adj = p.adjust(h3, 'BH')) %>%
  mutate(h6_adj = p.adjust(h6, 'BH')) %>%
  mutate(h24_adj = p.adjust(h24, 'BH')) %>%
  mutate(h48_adj = p.adjust(h48, 'BH')) %>%
  mutate(h72_adj = p.adjust(h72, 'BH')) %>%
  column_to_rownames('feature') %>%
  dplyr::select(contains('adj'))

#Check overall
apply(checkmat_neg, 2, function(x) sum(x <= 0.05))

#Adjust p-values
checkmat_neg <- lmer_contra_neg %>%
  mutate(h0_adj = p.adjust(h0, 'BH')) %>%
  mutate(h3_adj = p.adjust(h3, 'BH')) %>%
  mutate(h6_adj = p.adjust(h6, 'BH')) %>%
  mutate(h24_adj = p.adjust(h24, 'BH')) %>%
  mutate(h48_adj = p.adjust(h48, 'BH')) %>%
  mutate(h72_adj = p.adjust(h72, 'BH')) %>%
  column_to_rownames('feature') %>%
  dplyr::select(contains('adj'))

apply(checkmat_neg, 2, function(x) sum(x <= 0.05))

### Positive Mode -----------------------------------------------------------

#Run the model
lmer_mod_pos <- lmer_data_pos %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(abd ~ veg*time + grams_fed + (1|subject_id), data = x)))

#Run contrasts
lmer_contra_pos <- lmer_mod_pos %>%
  mutate(h0 = map_dbl(model, function(x) lmerTest::contest(x, BvA_0h)$`Pr(>F)`)) %>%
  mutate(h3 = map_dbl(model, function(x) lmerTest::contest(x, BvA_3h)$`Pr(>F)`)) %>%
  mutate(h6 = map_dbl(model, function(x) lmerTest::contest(x, BvA_6h)$`Pr(>F)`)) %>%
  mutate(h24 = map_dbl(model, function(x) lmerTest::contest(x, BvA_24h)$`Pr(>F)`)) %>%
  mutate(h48 = map_dbl(model, function(x) lmerTest::contest(x, BvA_48h)$`Pr(>F)`)) %>%
  mutate(h72 = map_dbl(model, function(x) lmerTest::contest(x, BvA_72h)$`Pr(>F)`)) %>%
  ungroup() 

#Adjust p-values
checkmat_pos <- lmer_contra_pos %>%
  mutate(h0_adj = p.adjust(h0, 'BH')) %>%
  mutate(h3_adj = p.adjust(h3, 'BH')) %>%
  mutate(h6_adj = p.adjust(h6, 'BH')) %>%
  mutate(h24_adj = p.adjust(h24, 'BH')) %>%
  mutate(h48_adj = p.adjust(h48, 'BH')) %>%
  mutate(h72_adj = p.adjust(h72, 'BH')) %>%
  column_to_rownames('feature') %>%
  dplyr::select(contains('adj'))

#Check overall
apply(checkmat_pos, 2, function(x) sum(x <= 0.05))




## Final Data Wrangling ----------------------------------------------------

#Import a list of our redundant features
ft_remove <- readxl::read_excel('./Data/Untargeted_Metabolomics/BSS_Adduct_Data.xlsx') %>%
  filter(`Keep/Remove`== 'Remove') %>%
  pull(Compound)

#Alter feature names so we can tell Pos and Neg mode apart
feature_map_neg <- data_neg %>%
  dplyr::select(Feature) %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_neg', Feature))) 

feature_map_pos <- data_pos %>%
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

data_clean_pos <- data_pos %>%
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
pos_sig <- checkmat_pos %>%
  rownames_to_column('feature') %>%
  filter(!feature %in% ft_remove) %>%
  pivot_longer(cols = ends_with('adj'), names_to = 'test', values_to = 'pval') %>%
  group_by(feature) %>%
  filter(any(pval <= 0.05)) %>%
  pivot_wider(names_from = 'test', values_from = 'pval') %>%
  pull(feature) %>%
  paste0('FT_', .) %>%
  gsub('(m/z|n)', '_pos', .)

pos_block <- data_clean_pos %>%
  dplyr::select(subject_id, treatment, time, veg, all_of(pos_sig)) 

neg_sig <- checkmat_neg %>%
  rownames_to_column('feature') %>%
  filter(!feature %in% ft_remove) %>%
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
metab_full <- left_join(pos_block, neg_block)

#Extract out individual timepoints for the PLS-DA
metab_3h <- metab_full %>%
  filter(time == 3) %>%
  dplyr::select(subject_id, veg, starts_with('FT'))

metab_6h <- metab_full %>%
  filter(time == 6) %>%
  dplyr::select(subject_id, veg, starts_with('FT'))

metab_24h <- metab_full %>%
  filter(time == 24) %>%
  dplyr::select(subject_id, veg, starts_with('FT'))

metab_48h <- metab_full %>%
  filter(time == 48) %>%
  dplyr::select(subject_id, veg, starts_with('FT'))

metab_72h <- metab_full %>%
  filter(time == 72) %>%
  dplyr::select(subject_id, veg, starts_with('FT'))

# Microbiome Prep ---------------------------------------------------------
library(phyloseq)
library(magrittr)

#Prep data for loading into phyloseq
asvtab <- readRDS('~/Documents/Projects/PhD/Short_Term_Broccoli/Data/BSS_Microbiome/seqtab_nochim.rds')
taxtab <- readRDS('~/Documents/Projects/PhD/Short_Term_Broccoli/Data/BSS_Microbiome/tax.rds')
metadata <- as.data.frame(read_csv('~/Documents/Projects/PhD/Short_Term_Broccoli/Data/BSS_Microbiome/metadata_microbiome.csv')) 
metadata$time %<>% factor(levels = c(0, 3, 24, 48, 72))
metadata$cohort %<>% factor()
#Load in information on how much sprouts each group ate
brocinfo <- read_csv(here('./Data/Metadata/sprout_meta.csv')) %>%
  modify_at('cohort', as.factor)
rownames(metadata) <- metadata$sample
rownames(asvtab) <- metadata$sample
ps_raw <- phyloseq(otu_table(asvtab, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxtab))


## Microbiome Cleaning -----------------------------------------------------

#Give arbitrary names to the taxa as opposed to keeping as just DNA-sequences which identify them
taxa_names(ps_raw) <- paste0("ASV", seq(ntaxa(ps_raw)))

#Fix the metadata
mdata_veg <- data.frame(sample_data(ps_raw)) %>%
  mutate(veg = ifelse(treatment %in% c('BL', 'BU'), 'broc', 'alf')) 
sample_data(ps_raw) <- mdata_veg

#Fill in missing genus names:
renames <- rownames(tax_table(ps_raw)[is.na(tax_table(ps_raw)[, 'Genus'])])
taxdf <- tax_table(ps_raw)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
tax_table(ps_raw)[renames, 'Genus'] <- renamed_genus

#Fill in missing species names
renames_sp <- rownames(tax_table(ps_raw)[is.na(tax_table(ps_raw)[, 'Species'])])
taxdf_sp <- tax_table(ps_raw)[renames_sp,]
renamed_species <- unname(sapply(taxa_names(taxdf_sp), function(x) paste0('g_', taxdf_sp[x, 'Genus'], '_', x)))
tax_table(ps_raw)[renames_sp, 'Species'] <- renamed_species

#ps_raw has 4060 ASVs
#Remove taxa that do not appear at least 3 times in at least 20% of all samples
ps_counts <- ps_raw %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE) #317
#Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) ) #317
#Filter out low abundance (>1e-5) taxa
ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE) #317

ps_c <- ps_counts 

subjects <- unique(metab_full$subject_id)

meta_micro <- sample_data(ps_c) %>%
  as.data.frame()

micro_clr <- ps_c %>%
  otu_table() %>%
  vegan::decostand('clr', pseudocount = 1) %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  left_join(meta_micro) %>%
  filter(!sample %in% c('BSS022_72h-1_S68', 'BSS032_3h_S103', 'BSS067_48h-2_S226', 'BSS083_72h-1_S283')) %>%
  filter(subject_id %in% subjects)

micro_0h <- micro_clr %>%
  filter(time == 0) %>%
  dplyr::select(subject_id, veg, starts_with('ASV'))

micro_24h <- micro_clr %>%
  filter(time == 24) %>%
  dplyr::select(subject_id, veg, starts_with('ASV'))

micro_48h <- micro_clr %>%
  filter(time == 48) %>%
  dplyr::select(subject_id, veg, starts_with('ASV'))

micro_72h <- micro_clr %>%
  filter(time == 72) %>%
  dplyr::select(subject_id, veg, starts_with('ASV'))

#Prep the dietary data
parse_diet_logs <- function(path){
  #Grab the subject ID
  ID <- suppressMessages(read_tsv(path, col_names = F)[[1,1]]) %>%
   gsub('Spreadsheet: ', '', .) %>%
   gsub('[ -]', '', .)
  #Read in the log
  testlog <- suppressMessages(read_tsv(path, skip = 2, na = '--'))
  #Pull out just the summary info
  summary_info <- testlog[str_detect(testlog$`Item Name`, 'Day'), ]
  summary_clean <- summary_info[,-c(2,3)] %>%
   mutate(rdsplit = str_split(string = `Item Name`, pattern = '\\(')) %>%
   mutate(relative_day = map_chr(rdsplit, function(x) gsub(pattern = ' ', replacement = '', x = x[1], ))) %>%
   dplyr::select(-`Item Name`, -rdsplit, -`Wgt (g)`) %>%
   pivot_longer(ends_with(')'), names_to = 'nutrient') %>%
   cbind(ID, .)
  return(summary_clean)
}

read_in_logs <- function(folder){
  map_df(folder, parse_diet_logs)
}

#Dietary Data
logs <- list.files(here('./Data/Metadata/Diet_logs/', full.names = T))
all_logs <- read_in_logs(list.files(here('./Data/Metadata/Diet_logs/'), full.names = T)) 

#Diet Data
diet_summary <- all_logs %>%
  #Keep only days leading up to Day 0
  filter(relative_day %in% paste0('Day', 1:7)) %>%
  group_by(ID, nutrient) %>%
  #Take the mean over the 7 days
  summarise(mean_prestudy = mean(value, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'nutrient', values_from = 'mean_prestudy') %>%
  rename_with(~gsub('\\)', '', gsub('\\(', '', gsub(' ', '_', .x)))) %>%
  #Add a prefix to each variable name for easier handling downstream
  rename_with(~paste0('diet_', .x)) %>%
  rename('subject_id' = 'diet_ID') %>%
  #Normalize diets to kCals consumed
  mutate(across(starts_with('diet_'), ~ .x/diet_Cals_kcal)) %>%
  dplyr::select(-diet_Cals_kcal) 
  
#Remove redundant information
diet_small <- diet_summary %>%
  dplyr::select(-diet_Fib16_g, -diet_FatCals_kcal, -diet_SolFib16_g, -diet_SatCals_kcal, -diet_Carb_g, -diet_Fat_g, -diet_Folate_mcg)



# Prep Data for Integration -----------------------------------------------
#Helpful Functions
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}

#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

#Auto Scaling:
AS_helper <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(mtb){
  mtb_scaled <- apply(mtb, 2, AS_helper) 
  return(mtb_scaled)
}

#Transformation Functions:
#Log Scaling:
log_transform <- function(mtb){
  mtb_nz <- mtb[ ,which(apply(mtb, 2, sum) != 0)]
  min.val <- min(abs(mtb_nz[mtb_nz!=0]))/10
  mtb_log_trans <- apply(mtb_nz, 2, log_helper, min.val)
  return(mtb_log_trans)
}


#Join to make sure samples are in the same order
full_3h <- left_join(metab_3h, micro_0h) %>%
  left_join(diet_small)
  
#Extract out the individual data blocks of interest
h3_sig_pos <- checkmat_pos %>% 
  rownames_to_column('Feature') %>%
  filter(!Feature %in% ft_remove) %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_pos', Feature))) %>%
  filter(h3_adj < 0.05) %>%
  dplyr::pull(fixed_name)

h3_sig_neg <- checkmat_neg %>% 
  rownames_to_column('Feature') %>%
  filter(!Feature %in% ft_remove) %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_neg', Feature))) %>%
  filter(h3_adj < 0.05) %>%
  dplyr::pull(fixed_name)

X_metab_3h <- full_3h %>%
  dplyr::select(all_of(c(h3_sig_pos, h3_sig_neg)))

X_micro_3h <- full_3h %>%
  dplyr::select(starts_with('ASV'))

X_diet_3h <- full_3h %>%
  dplyr::select(starts_with('diet_')) %>%
  log_transform() %>%
  auto_scale()

rownames(X_diet_3h) <- rownames(X_metab_3h)

Y_3h <- full_3h$veg

#Combine into a final data block
blocks_3h <- list(metab = X_metab_3h, micro = X_micro_3h, treatment = Y_3h) 

#Combine into a final data block
blocks_3h_diet <- list(metab = X_metab_3h, micro = X_micro_3h, diet = X_diet_3h, treatment = Y_3h)

# 3h Diablo ---------------------------------------------------------------
library(mixOmics)

#Metabolome and Microbiome
pls_mim <- spls(blocks_3h_diet$metab, blocks_3h_diet$micro, keepX = c(50, 50), keepY = c(50, 50), ncomp = 2, mode = 'canonical')
plotVar(pls_mim, cutoff = 0.5, title = 'Metabolome vs Microbiome', legend = c('Metabolome', 'Microbiome'), 
        var.names = F, pch = c(15,16), col = c('darkorchid', 'brown')) #Both positive and negative correlations
cor(pls_mim$variates$X, pls_mim$variates$Y) #0.80, 0.84

#Metabolome and Diet
pls_md <- spls(blocks_3h_diet$metab, blocks_3h_diet$diet, keepX = c(50, 50), keepY = c(25, 25), ncomp = 2, mode = 'canonical')
plotVar(pls_md, cutoff = 0.5, title = 'Metabolome vs Diet', legend = c('Metabolome', 'Diet'), 
        var.names = F, pch = c(16,17), col = c('darkorchid', 'green')) #Both positive and negative correlations
cor(pls_md$variates$X, pls_md$variates$Y) #0.54, 0.66

#Microbiome and Diet
pls_mid <- spls(blocks_3h_diet$micro, blocks_3h_diet$diet, keepX = c(50, 50), keepY = c(25, 25), ncomp = 2, mode = 'canonical')
plotVar(pls_mid, cutoff = 0.5, title = 'Microbiome vs Diet', legend = c('Metabolome', 'Microbiome'), 
        var.names = F, pch = c(15,17), col = c('brown', 'green')) #Both positive and negative correlations
cor(pls_mid$variates$X, pls_mid$variates$Y) #0.86, 0.86

#Metabolome and Microbiome
plsda_m <- splsda(blocks_3h_diet$metab, blocks_3h_diet$treatment, keepX = c(50, 50), ncomp = 2)
pmp <- perf(plsda_m, auc = T, method = 'loo')
pmp$error.rate #0.03125
pmp$auc #1

plsda_mi <- splsda(blocks_3h_diet$micro, blocks_3h_diet$treatment, keepX = c(50, 50), ncomp = 2)
pmip <- perf(plsda_mi, auc = T, method = 'loo')
pmip$error.rate #0.6875
pmip$auc #0.28

plsda_d <- splsda(blocks_3h_diet$diet, blocks_3h_diet$treatment, keepX = c(25, 25), ncomp = 2)
pdp <- perf(plsda_d, auc = T, method = 'loo')
pdp$error.rate #0.5
pdp$auc #0.48

## Model Tuning ------------------------------------------------------------

library(mixOmics)
library(pROC)

#Set up the model design
design <- matrix(1, ncol = 3, nrow = 3)
diag(design) <- 0

#Tune the model looking at correlations between components
#Run multiple models with different number of features kept to evaluate the sparsity
metab_test <- c(seq(5,175, 5))
micro_test <- c(seq(5,300, 5))
finalstats <- matrix(ncol = 4)
colnames(finalstats) <- c('comp1', 'comp2', 'param_m', 'param_mi')
for(dt in metab_test){
  for(mt in micro_test){
    keepX_temp <- list(metab = rep(dt, 2),
                       micro = rep(mt, 2))
    tempmod <- block.splsda(X = blocks_3h, indY =  3, ncomp = 2, keepX = keepX_temp, design = design)
  outstats <- cbind(t(diag(cor(tempmod$variates$metab, tempmod$variates$micro))),
      data.frame(param_m = dt,
                 param_mi = mt))
  finalstats <- rbind(finalstats, outstats)
  }
}

finalclean_3h <- finalstats %>%
  drop_na() %>%
  pivot_longer(cols = c('comp1', 'comp2'), names_to = 'comp', values_to = 'cor')

outplot_3h <- ggplot(finalclean_3h, aes(x = param_m, y = param_mi, size = cor, color = cor)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_3h) #MvD wants high but MvS wants low - Split the difference

tuneX <- list(metab =  metab_test, micro = micro_test)

tuner_3h <- tune.block.splsda(X = blocks_3h, indY = 3, ncomp = 2, test.keepX = tuneX, design = design,
                           progressBar = T, validation = 'Mfold')

cer_3h <- as.data.frame(tuner_3h$error.rate) %>%
  rownames_to_column('param') %>% 
  mutate(nms = str_split(param, '_')) %>%
  mutate(param_m = map_chr(nms, function(x) x[1])) %>%
  mutate(param_mi = map_chr(nms, function(x) x[2])) %>%
  dplyr::select(param_m, param_mi, comp1, comp2) %>%
  pivot_longer(c('comp1', 'comp2'), names_to = 'comp', values_to = 'err') %>%
  modify_at(c('param_m', 'param_mi', 'err'), as.numeric)

outplot_cer_3h <- ggplot(cer_3h, aes(x = param_m, y = param_mi, size = err, color = err)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_cer_3h) #120 by 120 wins again

#Run final model
diablo_3h_final <- block.splsda(blocks_3h, indY = 3, ncomp = 2, 
                                keepX = list(metab = c(120, 105),
                                             micro = c(230, 235)),
                                design = design)


# Generate ROC Curves -----------------------------------------------------
#Perf Diablo 
perf_diablo_3h <- perf(diablo_3h_final, validation = 'loo', auc = T)

#Generate Class Info
classmat <- as.data.frame(Y_3h) %>%
  mutate(alf = ifelse(Y_3h == 'alf', 1, 0)) %>%
  mutate(broc = ifelse(Y_3h == 'broc', 1, 0)) %>%
  dplyr::select(-Y_3h)

#Microbiome curves
#Comp 1
micro_3h_comp1 <- as.data.frame(perf_diablo_3h$predict$nrep1$micro$comp1)
#Generate ROC Curve
micro_3h_roc_comp1 <- roc(classmat$broc, micro_3h_comp1$broc)
#Get the AUC
auc(micro_3h_roc_comp1)
ggroc(micro_3h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_3h_roc_comp1), 3))) +
  ggtitle('Microbiome Comp 1')

#Comp 2
micro_3h_comp2 <- as.data.frame(perf_diablo_3h$predict$nrep1$micro$comp2)
#Generate ROC Curve
micro_3h_roc_comp2 <- roc(classmat$broc, micro_3h_comp2$broc)
#Get the AUC
auc(micro_3h_roc_comp2)
ggroc(micro_3h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_3h_roc_comp2), 3))) +
  ggtitle('Microbiome Comp 2')
  
#Metabolome curves
#Comp 1
metab_3h_comp1 <- as.data.frame(perf_diablo_3h$predict$nrep1$metab$comp1)
#Generate ROC Curve
metab_3h_roc_comp1 <- roc(classmat$broc, metab_3h_comp1$broc)
#Get the AUC
auc(metab_3h_roc_comp1)
ggroc(metab_3h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_3h_roc_comp1), 3))) +
  ggtitle('Metabolome Comp 1')

#Comp 2
metab_3h_comp2 <- as.data.frame(perf_diablo_3h$predict$nrep1$metab$comp2)
#Generate ROC Curve
metab_3h_roc_comp2 <- roc(classmat$broc, metab_3h_comp2$broc)
#Get the AUC
auc(metab_3h_roc_comp2)
ggroc(metab_3h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_3h_roc_comp2), 3))) +
  ggtitle('Metabolome Comp 2')

vegPal <- c('orchid4', 'darkgreen')
names(vegPal) <- c('alf', 'broc')


selectVar(diablo_3h_final)

plotIndiv(diablo_3h_final, legend = T)
plotLoadings(diablo_3h_final)
circosPlot(diablo_3h_final, cutoff = 0.6, showIntraLinks = F, line = T, size.variables = 0.5,
           color.Y = vegPal)


metabLoadings1_3h <- plotLoadings(diablo_3h_final, block = 'metab', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

metabLoadings2_3h <- plotLoadings(diablo_3h_final, block = 'metab', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

microLoadings1_3h <- plotLoadings(diablo_3h_final, block = 'micro', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

microLoadings2_3h <- plotLoadings(diablo_3h_final, block = 'micro', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

diablo_out_3h <- rbind(metabLoadings1_3h, metabLoadings2_3h) %>%
  rbind(., microLoadings1_3h) %>%
  rbind(., microLoadings2_3h)

dcir_3h <- mixOmics::plotVar(diablo_3h_final, style = 'ggplot2', var.names = T, plot = F) %>%
  filter(names %in% diablo_out_3h$feature)

#Circle function to make plotting easier
circleFun <- function(center = c(-1,1),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#Make the circle
circle <- circleFun(c(0,0), 2, npoints = 100)

#Put it all together
dcir2_3h <- dcir_3h %>%
  dplyr::select(x, y, names) %>%
  left_join(diablo_out_3h, by = c(c('names' = 'feature')))



#Plot it
corcir_3h <- ggplot(dcir2_3h, aes(x, y, color = GroupContrib, shape = Block)) + 
  geom_point(size = 4, aes(text = names)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1) +
  xlab('Component 1') +
  ylab('Component 2') +
  ggtitle('Correlation Circle')  +
  #Make the circle independent of other aesthetics
  geom_path(aes(x,y), data = circle, inherit.aes = F)  +
  scale_color_manual(values = vegPal) +
  scale_y_continuous(labels = signs::signs_format()) +
  scale_x_continuous(labels = signs::signs_format()) +
  theme(legend.position = 'none')

plotly::ggplotly(corcir_3h, tooltip = 'text')

alltax <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column('ASV')

#Heatmap
hData_3h <- full_3h %>%
  filter(veg == 'broc') %>%
  dplyr::select(all_of(diablo_out_3h$feature)) %>%
  pivot_longer(starts_with('FT_'), names_to = 'feature', values_to = 'int') %>%
  pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'clr') %>% 
  group_by(feature, ASV) %>%
  nest() %>%
  mutate(rho = map_dbl(data, function(x) cor(x$int, x$clr))) %>%
  dplyr::select(-data) %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(test = map_lgl(data, function(x) any(abs(x$rho) > 0.6))) %>%
  filter(test) %>%
  dplyr::select(-test) %>%
  unnest('data') %>%
  group_by(feature) %>%
  nest() %>%
  mutate(test = map_lgl(data, function(x) any(abs(x$rho) > 0.6))) %>%
  filter(test) %>%
  dplyr::select(-test) %>%
  unnest('data') %>%
  left_join(alltax) %>%
  mutate(gsa = ifelse(str_detect(Species, 'g_'), paste0(Genus, '_', Species),
                      paste0(Genus, '_', Species, '_', ASV))) %>%
  ungroup() %>%
  dplyr::select(feature, gsa, rho) %>%
  pivot_wider(names_from = 'gsa', values_from = 'rho') %>%
  column_to_rownames('feature') 

pheatmap::pheatmap(hData_3h)



# 6h Diablo ---------------------------------------------------------------

#Join to make sure samples are in the same order
full_6h <- left_join(metab_6h, micro_0h)

#Extract out the individual data blocks of interest
sig_pos_6h <- checkmat_pos %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_pos', Feature))) %>%
  filter(h6_adj < 0.05) %>%
  dplyr::pull(fixed_name)

sig_neg_6h <- checkmat_neg %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_neg', Feature))) %>%
  filter(h6_adj < 0.05) %>%
  dplyr::pull(fixed_name)

X_metab_6h <- full_6h %>%
  dplyr::select(all_of(c(sig_pos_6h, sig_neg_6h)))

X_micro_6h <- full_6h %>%
  dplyr::select(starts_with('ASV'))

Y_6h <- full_6h$veg

#Combine into a final data block
blocks_6h <- list(metab = X_metab_6h, micro = X_micro_6h, treatment = Y_6h)

## Model Tuning ------------------------------------------------------------
#Tune the model looking at correlations between components
#Run multiple models with different number of features kept to evaluate the sparsity
metab_test <- c(seq(5,248, 5))
micro_test <- c(seq(5,300, 5))
finalstats <- matrix(ncol = 4)
colnames(finalstats) <- c('comp1', 'comp2', 'param_m', 'param_mi')
for(dt in metab_test){
  for(mt in micro_test){
    keepX_temp <- list(metab = rep(dt, 2),
                       micro = rep(mt, 2))
    tempmod <- block.splsda(X = blocks_6h, indY =  3, ncomp = 2, keepX = keepX_temp, design = design)
  outstats <- cbind(t(diag(cor(tempmod$variates$metab, tempmod$variates$micro))),
      data.frame(param_m = dt,
                 param_mi = mt))
  finalstats <- rbind(finalstats, outstats)
  }
}

finalclean_6h <- finalstats %>%
  drop_na() %>%
  pivot_longer(cols = c('comp1', 'comp2'), names_to = 'comp', values_to = 'cor')

outplot_6h <- ggplot(finalclean_6h, aes(x = param_m, y = param_mi, size = cor, color = cor)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_6h) #MvD wants high but MvS wants low - Split the difference

tuneX <- list(metab =  metab_test, micro = micro_test)

tuner_6h <- tune.block.splsda(X = blocks_6h, indY = 3, ncomp = 2, test.keepX = tuneX, design = design,
                           progressBar = T, validation = 'Mfold')

cer_6h <- as.data.frame(tuner_6h$error.rate) %>%
  rownames_to_column('param') %>% 
  mutate(nms = str_split(param, '_')) %>%
  mutate(param_m = map_chr(nms, function(x) x[1])) %>%
  mutate(param_mi = map_chr(nms, function(x) x[2])) %>%
  dplyr::select(param_m, param_mi, comp1, comp2) %>%
  pivot_longer(c('comp1', 'comp2'), names_to = 'comp', values_to = 'err') %>%
  modify_at(c('param_m', 'param_mi', 'err'), as.numeric)

outplot_cer_6h <- ggplot(cer_6h, aes(x = param_m, y = param_mi, size = err, color = err)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_cer_6h) #120 by 120 wins again

#Run final model
diablo_6h_final <- block.splsda(blocks_6h, indY = 3, ncomp = 2, 
                                keepX = list(metab = c(120, 50),
                                             micro = c(130, 75)),
                                design = design)


# Generate ROC Curves -----------------------------------------------------
#Perf Diablo 
perf_diablo_6h <- perf(diablo_6h_final, validation = 'loo', auc = T)

#Generate Class Info
classmat <- as.data.frame(Y_6h) %>%
  mutate(alf = ifelse(Y_6h == 'alf', 1, 0)) %>%
  mutate(broc = ifelse(Y_6h == 'broc', 1, 0)) %>%
  dplyr::select(-Y_6h)

#Microbiome curves
#Comp 1
micro_6h_comp1 <- as.data.frame(perf_diablo_6h$predict$nrep1$micro$comp1)
#Generate ROC Curve
micro_6h_roc_comp1 <- roc(classmat$broc, micro_6h_comp1$broc)
#Get the AUC
auc(micro_6h_roc_comp1)
ggroc(micro_6h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_6h_roc_comp1), 3))) +
  ggtitle('Microbiome Comp 1')

#Comp 2
micro_6h_comp2 <- as.data.frame(perf_diablo_6h$predict$nrep1$micro$comp2)
#Generate ROC Curve
micro_6h_roc_comp2 <- roc(classmat$broc, micro_6h_comp2$broc)
#Get the AUC
auc(micro_6h_roc_comp2)
ggroc(micro_6h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_6h_roc_comp2), 3))) +
  ggtitle('Microbiome Comp 2')
  
#Metabolome curves
#Comp 1
metab_6h_comp1 <- as.data.frame(perf_diablo_6h$predict$nrep1$metab$comp1)
#Generate ROC Curve
metab_6h_roc_comp1 <- roc(classmat$broc, metab_6h_comp1$broc)
#Get the AUC
auc(metab_6h_roc_comp1)
ggroc(metab_6h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_6h_roc_comp1), 3))) +
  ggtitle('Metabolome Comp 1')

#Comp 2
metab_6h_comp2 <- as.data.frame(perf_diablo_6h$predict$nrep1$metab$comp2)
#Generate ROC Curve
metab_6h_roc_comp2 <- roc(classmat$broc, metab_6h_comp2$broc)
#Get the AUC
auc(metab_6h_roc_comp2)
ggroc(metab_6h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_6h_roc_comp2), 3))) +
  ggtitle('Metabolome Comp 2')


plotIndiv(diablo_6h_final, legend = T)
plotArrow(diablo_6h_final)
plotIndiv(diablo_6h_final, legend = T)
plotLoadings(diablo_6h_final)
plotVar(diablo_6h_final)
circosPlot(diablo_6h_final, cutoff = 0.6, showIntraLinks = F, line = T, size.variables = 0.5,
           color.Y = vegPal)

metabLoadings1_6h <- plotLoadings(diablo_6h_final, block = 'metab', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

metabLoadings2_6h <- plotLoadings(diablo_6h_final, block = 'metab', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

microLoadings1_6h <- plotLoadings(diablo_6h_final, block = 'micro', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

microLoadings2_6h <- plotLoadings(diablo_6h_final, block = 'micro', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

diablo_out_6h <- rbind(metabLoadings1_6h, metabLoadings2_6h) %>%
  rbind(., microLoadings1_6h) %>%
  rbind(., microLoadings2_6h)

dcir_6h <- mixOmics::plotVar(diablo_6h_final, style = 'ggplot2', var.names = T, plot = F) %>%
  filter(names %in% diablo_out_6h$feature)

#Circle function to make plotting easier
circleFun <- function(center = c(-1,1),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#Make the circle
circle <- circleFun(c(0,0), 2, npoints = 100)

#Put it all together
dcir2_6h <- dcir_6h %>%
  dplyr::select(x, y, names) %>%
  left_join(diablo_out_6h, by = c(c('names' = 'feature')))

vegPal <- c('orchid4', 'darkgreen')
names(vegPal) <- c('alf', 'broc')

#Plot it
corcir_6h <- ggplot(dcir2_6h, aes(x, y, color = GroupContrib, shape = Block)) + 
  geom_point(size = 4, aes(text = names)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1) +
  xlab('Component 1') +
  ylab('Component 2') +
  ggtitle('Correlation Circle')  +
  #Make the circle independent of other aesthetics
  geom_path(aes(x,y), data = circle, inherit.aes = F)  +
  scale_color_manual(values = vegPal) +
  scale_y_continuous(labels = signs::signs_format()) +
  scale_x_continuous(labels = signs::signs_format()) +
  theme(legend.position = 'none')

plotly::ggplotly(corcir_6h, tooltip = 'text')

# 24h Diablo ---------------------------------------------------------------

#Join to make sure samples are in the same order
full_24h <- left_join(metab_24h, micro_0h)
full_24h2 <- left_join(metab_24h, micro_24h)

#Extract out the individual data blocks of interest
sig_pos_24h <- checkmat_pos %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_pos', Feature))) %>%
  filter(h24_adj < 0.05) %>%
  dplyr::pull(fixed_name)

sig_neg_24h <- checkmat_neg %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_neg', Feature))) %>%
  filter(h24_adj < 0.05) %>%
  dplyr::pull(fixed_name)

X_metab_24h <- full_24h %>%
  dplyr::select(all_of(c(sig_pos_24h, sig_neg_24h)))

X_micro_24h <- full_24h %>%
  dplyr::select(starts_with('ASV'))


Y_24h <- full_24h$veg

X_metab_24h2 <- full_24h2 %>%
  dplyr::select(all_of(c(sig_pos_24h, sig_neg_24h)))

X_micro_24h2 <- full_24h2 %>%
  dplyr::select(starts_with('ASV'))


Y_24h2 <- full_24h2$veg

#Combine into a final data block
blocks_24h <- list(metab = X_metab_24h, micro = X_micro_24h, treatment = Y_24h)
blocks_24h2 <- list(metab = X_metab_24h2, micro = X_micro_24h2, treatment = Y_24h2)

## Model Tuning ------------------------------------------------------------
#Tune the model looking at correlations between components
#Run multiple models with different number of features kept to evaluate the sparsity
metab_test <- c(seq(5,100, 4))
micro_test <- c(seq(5,300, 5))
finalstats <- matrix(ncol = 4)
colnames(finalstats) <- c('comp1', 'comp2', 'param_m', 'param_mi')
for(dt in metab_test){
  for(mt in micro_test){
    keepX_temp <- list(metab = rep(dt, 2),
                       micro = rep(mt, 2))
    tempmod <- block.splsda(X = blocks_24h2, indY =  3, ncomp = 2, keepX = keepX_temp, design = design)
  outstats <- cbind(t(diag(cor(tempmod$variates$metab, tempmod$variates$micro))),
      data.frame(param_m = dt,
                 param_mi = mt))
  finalstats <- rbind(finalstats, outstats)
  }
}

finalclean_24h <- finalstats %>%
  drop_na() %>%
  pivot_longer(cols = c('comp1', 'comp2'), names_to = 'comp', values_to = 'cor')

outplot_24h <- ggplot(finalclean_24h, aes(x = param_m, y = param_mi, size = cor, color = cor)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_24h) #MvD wants high but MvS wants low - Split the difference

tuneX <- list(metab =  metab_test, micro = micro_test)

tuner_24h <- tune.block.splsda(X = blocks_24h2, indY = 3, ncomp = 2, test.keepX = tuneX, design = design,
                           progressBar = T, validation = 'Mfold')

cer_24h <- as.data.frame(tuner_24h$error.rate) %>%
  rownames_to_column('param') %>% 
  mutate(nms = str_split(param, '_')) %>%
  mutate(param_m = map_chr(nms, function(x) x[1])) %>%
  mutate(param_mi = map_chr(nms, function(x) x[2])) %>%
  dplyr::select(param_m, param_mi, comp1, comp2) %>%
  pivot_longer(c('comp1', 'comp2'), names_to = 'comp', values_to = 'err') %>%
  modify_at(c('param_m', 'param_mi', 'err'), as.numeric)

outplot_cer_24h <- ggplot(cer_24h, aes(x = param_m, y = param_mi, size = err, color = err)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_cer_24h) #120 by 120 wins again

#Run final model
diablo_24h_final <- block.splsda(blocks_24h, indY = 3, ncomp = 2, 
                                keepX = list(metab = c(89, 9),
                                             micro = c(190, 45)),
                                design = design)

diablo_24h_final2 <- block.splsda(blocks_24h2, indY = 3, ncomp = 2, 
                                keepX = list(metab = c(89, 9),
                                             micro = c(190, 45)),
                                design = design)

# Generate ROC Curves -----------------------------------------------------
#Perf Diablo 
perf_diablo_24h <- perf(diablo_24h_final, validation = 'loo', auc = T)
perf_diablo_24h2 <- perf(diablo_24h_final2, validation = 'loo', auc = T)

#Generate Class Info
classmat <- as.data.frame(Y_24h) %>%
  mutate(alf = ifelse(Y_24h == 'alf', 1, 0)) %>%
  mutate(broc = ifelse(Y_24h == 'broc', 1, 0)) %>%
  dplyr::select(-Y_24h)

#Microbiome curves
#Comp 1
micro_24h_comp1 <- as.data.frame(perf_diablo_24h$predict$nrep1$micro$comp1)
#Generate ROC Curve
micro_24h_roc_comp1 <- roc(classmat$broc, micro_24h_comp1$broc)
#Get the AUC
auc(micro_24h_roc_comp1)
ggroc(micro_24h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_24h_roc_comp1), 3))) +
  ggtitle('Microbiome Comp 1')

#Microbiome curves
#Comp 1
micro_24h_comp12 <- as.data.frame(perf_diablo_24h2$predict$nrep1$micro$comp1)
#Generate ROC Curve
micro_24h_roc_comp12 <- roc(classmat$broc, micro_24h_comp12$broc)
#Get the AUC
auc(micro_24h_roc_comp12)
ggroc(micro_24h_roc_comp12, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_24h_roc_comp12), 3))) +
  ggtitle('Microbiome Comp 1')

#Comp 2
micro_24h_comp2 <- as.data.frame(perf_diablo_24h$predict$nrep1$micro$comp2)
#Generate ROC Curve
micro_24h_roc_comp2 <- roc(classmat$broc, micro_24h_comp2$broc)
#Get the AUC
auc(micro_24h_roc_comp2)
ggroc(micro_24h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_24h_roc_comp2), 3))) +
  ggtitle('Microbiome Comp 2')
  
#Metabolome curves
#Comp 1
metab_24h_comp1 <- as.data.frame(perf_diablo_24h$predict$nrep1$metab$comp1)
#Generate ROC Curve
metab_24h_roc_comp1 <- roc(classmat$broc, metab_24h_comp1$broc)
#Get the AUC
auc(metab_24h_roc_comp1)
ggroc(metab_24h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_24h_roc_comp1), 3))) +
  ggtitle('Metabolome Comp 1')

#Comp 2
metab_24h_comp2 <- as.data.frame(perf_diablo_24h$predict$nrep1$metab$comp2)
#Generate ROC Curve
metab_24h_roc_comp2 <- roc(classmat$broc, metab_24h_comp2$broc)
#Get the AUC
auc(metab_24h_roc_comp2)
ggroc(metab_24h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_24h_roc_comp2), 3))) +
  ggtitle('Metabolome Comp 2')


plotIndiv(diablo_24h_final, legend = T)
plotArrow(diablo_24h_final)
plotIndiv(diablo_24h_final, legend = T)
plotLoadings(diablo_24h_final)
plotVar(diablo_24h_final)
circosPlot(diablo_24h_final, cutoff = 0.6, showIntraLinks = F, line = T, size.variables = 0.5,
           color.Y = vegPal)

metabLoadings1_24h <- plotLoadings(diablo_24h_final, block = 'metab', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

metabLoadings2_24h <- plotLoadings(diablo_24h_final, block = 'metab', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

microLoadings1_24h <- plotLoadings(diablo_24h_final, block = 'micro', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

microLoadings2_24h <- plotLoadings(diablo_24h_final, block = 'micro', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

diablo_out_24h <- rbind(metabLoadings1_24h, metabLoadings2_24h) %>%
  rbind(., microLoadings1_24h) %>%
  rbind(., microLoadings2_24h)

dcir_24h <- mixOmics::plotVar(diablo_24h_final, style = 'ggplot2', var.names = T, plot = F) %>%
  filter(names %in% diablo_out_24h$feature)

#Put it all together
dcir2_24h <- dcir_24h %>%
  dplyr::select(x, y, names) %>%
  left_join(diablo_out_24h, by = c(c('names' = 'feature')))

#Plot it
corcir_24h <- ggplot(dcir2_24h, aes(x, y, color = GroupContrib, shape = Block)) + 
  geom_point(size = 4, aes(text = names)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1) +
  xlab('Component 1') +
  ylab('Component 2') +
  ggtitle('Correlation Circle')  +
  #Make the circle independent of other aesthetics
  geom_path(aes(x,y), data = circle, inherit.aes = F)  +
  scale_color_manual(values = vegPal) +
  scale_y_continuous(labels = signs::signs_format()) +
  scale_x_continuous(labels = signs::signs_format()) +
  theme(legend.position = 'none')

plotly::ggplotly(corcir_24h, tooltip = 'text')

#Heatmap





# 48h Diablo --------------------------------------------------------------

#Join to make sure samples are in the same order
full_48h <- left_join(metab_48h, micro_0h)

#Extract out the individual data blocks of interest
sig_pos_48h <- checkmat_pos %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_pos', Feature))) %>%
  filter(h48_adj < 0.05) %>%
  dplyr::pull(fixed_name)

sig_neg_48h <- checkmat_neg %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_neg', Feature))) %>%
  filter(h48_adj < 0.05) %>%
  dplyr::pull(fixed_name)

X_metab_48h <- full_48h %>%
  dplyr::select(all_of(c(sig_pos_48h, sig_neg_48h)))

X_micro_48h <- full_48h %>%
  dplyr::select(starts_with('ASV'))

Y_48h <- full_48h$veg

#Combine into a final data block
blocks_48h <- list(metab = X_metab_48h, micro = X_micro_48h, treatment = Y_48h)

## Model Tuning ------------------------------------------------------------
#Tune the model looking at correlations between components
#Run multiple models with different number of features kept to evaluate the sparsity
metab_test <- c(seq(5,50, 2))
micro_test <- c(seq(5,300, 5))
finalstats <- matrix(ncol = 4)
colnames(finalstats) <- c('comp1', 'comp2', 'param_m', 'param_mi')
for(dt in metab_test){
  for(mt in micro_test){
    keepX_temp <- list(metab = rep(dt, 2),
                       micro = rep(mt, 2))
    tempmod <- block.splsda(X = blocks_48h, indY =  3, ncomp = 2, keepX = keepX_temp, design = design)
  outstats <- cbind(t(diag(cor(tempmod$variates$metab, tempmod$variates$micro))),
      data.frame(param_m = dt,
                 param_mi = mt))
  finalstats <- rbind(finalstats, outstats)
  }
}

finalclean_48h <- finalstats %>%
  drop_na() %>%
  pivot_longer(cols = c('comp1', 'comp2'), names_to = 'comp', values_to = 'cor')

outplot_48h <- ggplot(finalclean_48h, aes(x = param_m, y = param_mi, size = cor, color = cor)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_48h) #MvD wants high but MvS wants low - Split the difference

tuneX <- list(metab =  metab_test, micro = micro_test)

tuner_48h <- tune.block.splsda(X = blocks_48h, indY = 3, ncomp = 2, test.keepX = tuneX, design = design,
                           progressBar = T, validation = 'Mfold')

cer_48h <- as.data.frame(tuner_48h$error.rate) %>%
  rownames_to_column('param') %>% 
  mutate(nms = str_split(param, '_')) %>%
  mutate(param_m = map_chr(nms, function(x) x[1])) %>%
  mutate(param_mi = map_chr(nms, function(x) x[2])) %>%
  dplyr::select(param_m, param_mi, comp1, comp2) %>%
  pivot_longer(c('comp1', 'comp2'), names_to = 'comp', values_to = 'err') %>%
  modify_at(c('param_m', 'param_mi', 'err'), as.numeric)

outplot_cer_48h <- ggplot(cer_48h, aes(x = param_m, y = param_mi, size = err, color = err)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_cer_48h) #120 by 120 wins again

#Run final model
diablo_48h_final <- block.splsda(blocks_48h, indY = 3, ncomp = 2, 
                                keepX = list(metab = c(7, 29),
                                             micro = c(210, 40)),
                                design = design)


# Generate ROC Curves -----------------------------------------------------
#Perf Diablo 
perf_diablo_48h <- perf(diablo_48h_final, validation = 'loo', auc = T)

#Generate Class Info
classmat <- as.data.frame(Y_48h) %>%
  mutate(alf = ifelse(Y_48h == 'alf', 1, 0)) %>%
  mutate(broc = ifelse(Y_48h == 'broc', 1, 0)) %>%
  dplyr::select(-Y_48h)

#Microbiome curves
#Comp 1
micro_48h_comp1 <- as.data.frame(perf_diablo_48h$predict$nrep1$micro$comp1)
#Generate ROC Curve
micro_48h_roc_comp1 <- roc(classmat$broc, micro_48h_comp1$broc)
#Get the AUC
auc(micro_48h_roc_comp1)
ggroc(micro_48h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_48h_roc_comp1), 3))) +
  ggtitle('Microbiome Comp 1')

#Comp 2
micro_48h_comp2 <- as.data.frame(perf_diablo_48h$predict$nrep1$micro$comp2)
#Generate ROC Curve
micro_48h_roc_comp2 <- roc(classmat$broc, micro_48h_comp2$broc)
#Get the AUC
auc(micro_48h_roc_comp2)
ggroc(micro_48h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_48h_roc_comp2), 3))) +
  ggtitle('Microbiome Comp 2')
  
#Metabolome curves
#Comp 1
metab_48h_comp1 <- as.data.frame(perf_diablo_48h$predict$nrep1$metab$comp1)
#Generate ROC Curve
metab_48h_roc_comp1 <- roc(classmat$broc, metab_48h_comp1$broc)
#Get the AUC
auc(metab_48h_roc_comp1)
ggroc(metab_48h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_48h_roc_comp1), 3))) +
  ggtitle('Metabolome Comp 1')

#Comp 2
metab_48h_comp2 <- as.data.frame(perf_diablo_48h$predict$nrep1$metab$comp2)
#Generate ROC Curve
metab_48h_roc_comp2 <- roc(classmat$broc, metab_48h_comp2$broc)
#Get the AUC
auc(metab_48h_roc_comp2)
ggroc(metab_48h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_48h_roc_comp2), 3))) +
  ggtitle('Metabolome Comp 2')


plotIndiv(diablo_48h_final, legend = T)
plotArrow(diablo_48h_final)
plotLoadings(diablo_48h_final)
plotVar(diablo_48h_final)
circosPlot(diablo_48h_final, cutoff = 0.6, showIntraLinks = F, line = T, size.variables = 0.5,
           color.Y = vegPal)

metabLoadings1_48h <- plotLoadings(diablo_48h_final, block = 'metab', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

metabLoadings2_48h <- plotLoadings(diablo_48h_final, block = 'metab', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

microLoadings1_48h <- plotLoadings(diablo_48h_final, block = 'micro', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

microLoadings2_48h <- plotLoadings(diablo_48h_final, block = 'micro', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

diablo_out_48h <- rbind(metabLoadings1_48h, metabLoadings2_48h) %>%
  rbind(., microLoadings1_48h) %>%
  rbind(., microLoadings2_48h)

dcir_48h <- mixOmics::plotVar(diablo_48h_final, style = 'ggplot2', var.names = T, plot = F) %>%
  filter(names %in% diablo_out_48h$feature)

#Put it all together
dcir2_48h <- dcir_48h %>%
  dplyr::select(x, y, names) %>%
  left_join(diablo_out_48h, by = c(c('names' = 'feature')))

#Plot it
corcir_48h <- ggplot(dcir2_48h, aes(x, y, color = GroupContrib, shape = Block)) + 
  geom_point(size = 4, aes(text = names)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1) +
  xlab('Component 1') +
  ylab('Component 2') +
  ggtitle('Correlation Circle')  +
  #Make the circle independent of other aesthetics
  geom_path(aes(x,y), data = circle, inherit.aes = F)  +
  scale_color_manual(values = vegPal) +
  scale_y_continuous(labels = signs::signs_format()) +
  scale_x_continuous(labels = signs::signs_format()) +
  theme(legend.position = 'none')

plotly::ggplotly(corcir_48h, tooltip = 'text')

# 72h Diablo --------------------------------------------------------------

#Join to make sure samples are in the same order
full_72h <- left_join(metab_72h, micro_48h)

#Extract out the individual data blocks of interest
sig_pos_72h <- checkmat_pos %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_pos', Feature))) %>%
  filter(h72_adj < 0.05) %>%
  dplyr::pull(fixed_name)

sig_neg_72h <- checkmat_neg %>% 
  rownames_to_column('Feature') %>%
  mutate(fixed_name = paste0('FT_', gsub('(m/z|n)', '_neg', Feature))) %>%
  filter(h72_adj < 0.05) %>%
  dplyr::pull(fixed_name)

X_metab_72h <- full_72h %>%
  dplyr::select(all_of(c(sig_pos_72h, sig_neg_72h)))

X_micro_72h <- full_72h %>%
  dplyr::select(starts_with('ASV'))

Y_72h <- full_72h$veg

#Combine into a final data block
blocks_72h <- list(metab = X_metab_72h, micro = X_micro_72h, treatment = Y_72h)

## Model Tuning ------------------------------------------------------------
#Tune the model looking at correlations between components
#Run multiple models with different number of features kept to evaluate the sparsity
metab_test <- c(seq(1,10, 1))
micro_test <- c(seq(5,300, 5))
finalstats <- matrix(ncol = 4)
colnames(finalstats) <- c('comp1', 'comp2', 'param_m', 'param_mi')
for(dt in metab_test){
  for(mt in micro_test){
    keepX_temp <- list(metab = rep(dt, 2),
                       micro = rep(mt, 2))
    tempmod <- block.splsda(X = blocks_72h, indY =  3, ncomp = 2, keepX = keepX_temp, design = design)
  outstats <- cbind(t(diag(cor(tempmod$variates$metab, tempmod$variates$micro))),
      data.frame(param_m = dt,
                 param_mi = mt))
  finalstats <- rbind(finalstats, outstats)
  }
}

finalclean_72h <- finalstats %>%
  drop_na() %>%
  pivot_longer(cols = c('comp1', 'comp2'), names_to = 'comp', values_to = 'cor')

outplot_72h <- ggplot(finalclean_72h, aes(x = param_m, y = param_mi, size = cor, color = cor)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_72h) #MvD wants high but MvS wants low - Split the difference

tuneX <- list(metab =  metab_test, micro = micro_test)

tuner_72h <- tune.block.splsda(X = blocks_72h, indY = 3, ncomp = 2, test.keepX = tuneX, design = design,
                           progressBar = T, validation = 'Mfold')

cer_72h <- as.data.frame(tuner_72h$error.rate) %>%
  rownames_to_column('param') %>% 
  mutate(nms = str_split(param, '_')) %>%
  mutate(param_m = map_chr(nms, function(x) x[1])) %>%
  mutate(param_mi = map_chr(nms, function(x) x[2])) %>%
  dplyr::select(param_m, param_mi, comp1, comp2) %>%
  pivot_longer(c('comp1', 'comp2'), names_to = 'comp', values_to = 'err') %>%
  modify_at(c('param_m', 'param_mi', 'err'), as.numeric)

outplot_cer_72h <- ggplot(cer_72h, aes(x = param_m, y = param_mi, size = err, color = err)) +
  geom_point() +
  facet_wrap(~comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot_cer_72h) #120 by 120 wins again

#Run final model
diablo_72h_final <- block.splsda(blocks_72h, indY = 3, ncomp = 2, 
                                keepX = list(metab = c(2, 3),
                                             micro = c(3, 25)),
                                design = design)


# Generate ROC Curves -----------------------------------------------------
#Perf Diablo 
perf_diablo_72h <- perf(diablo_72h_final, validation = 'loo', auc = T)

#Generate Class Info
classmat <- as.data.frame(Y_72h) %>%
  mutate(alf = ifelse(Y_72h == 'alf', 1, 0)) %>%
  mutate(broc = ifelse(Y_72h == 'broc', 1, 0)) %>%
  dplyr::select(-Y_72h)

#Microbiome curves
#Comp 1
micro_72h_comp1 <- as.data.frame(perf_diablo_72h$predict$nrep1$micro$comp1)
#Generate ROC Curve
micro_72h_roc_comp1 <- roc(classmat$broc, micro_72h_comp1$broc)
#Get the AUC
auc(micro_72h_roc_comp1)
ggroc(micro_72h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_72h_roc_comp1), 3))) +
  ggtitle('Microbiome Comp 1')

#Comp 2
micro_72h_comp2 <- as.data.frame(perf_diablo_72h$predict$nrep1$micro$comp2)
#Generate ROC Curve
micro_72h_roc_comp2 <- roc(classmat$broc, micro_72h_comp2$broc)
#Get the AUC
auc(micro_72h_roc_comp2)
ggroc(micro_72h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(micro_72h_roc_comp2), 3))) +
  ggtitle('Microbiome Comp 2')
  
#Metabolome curves
#Comp 1
metab_72h_comp1 <- as.data.frame(perf_diablo_72h$predict$nrep1$metab$comp1)
#Generate ROC Curve
metab_72h_roc_comp1 <- roc(classmat$broc, metab_72h_comp1$broc)
#Get the AUC
auc(metab_72h_roc_comp1)
ggroc(metab_72h_roc_comp1, size = 2)  +
  geom_abline(intercept = 1) + 
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_72h_roc_comp1), 3))) +
  ggtitle('Metabolome Comp 1')

#Comp 2
metab_72h_comp2 <- as.data.frame(perf_diablo_72h$predict$nrep1$metab$comp2)
#Generate ROC Curve
metab_72h_roc_comp2 <- roc(classmat$broc, metab_72h_comp2$broc)
#Get the AUC
auc(metab_72h_roc_comp2)
ggroc(metab_72h_roc_comp2, size = 2)  +
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot()  +
  annotate('text', x = 0.25, y = 0.25, label = paste0('AUC: ', round(auc(metab_72h_roc_comp2), 3))) +
  ggtitle('Metabolome Comp 2')


plotIndiv(diablo_72h_final, legend = T)
plotArrow(diablo_72h_final)
plotLoadings(diablo_72h_final)
plotVar(diablo_72h_final)
circosPlot(diablo_72h_final, cutoff = 0.6, showIntraLinks = F, line = T, size.variables = 0.5,
           color.Y = vegPal)

metabLoadings1_72h <- plotLoadings(diablo_72h_final, block = 'metab', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

metabLoadings2_72h <- plotLoadings(diablo_72h_final, block = 'metab', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Metabolomics') %>%
  filter(abs(importance) >= 0.1)

microLoadings1_72h <- plotLoadings(diablo_72h_final, block = 'micro', comp = 1, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 1) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

microLoadings2_72h <- plotLoadings(diablo_72h_final, block = 'micro', comp = 2, contrib = 'max', method = 'median', plot = F) %>%
  dplyr::select(importance, GroupContrib) %>%
  rownames_to_column('feature') %>%
  magrittr::inset('Component', value = 2) %>%
  magrittr::inset('Block', value = 'Microbiome') %>%
  filter(abs(importance) >= 0.1)

diablo_out_72h <- rbind(metabLoadings1_72h, metabLoadings2_72h) %>%
  rbind(., microLoadings1_72h) %>%
  rbind(., microLoadings2_72h)

dcir_72h <- mixOmics::plotVar(diablo_72h_final, style = 'ggplot2', var.names = T, plot = F) %>%
  filter(names %in% diablo_out_72h$feature)

#Put it all together
dcir2_72h <- dcir_72h %>%
  dplyr::select(x, y, names) %>%
  left_join(diablo_out_72h, by = c(c('names' = 'feature')))

#Plot it
corcir_72h <- ggplot(dcir2_72h, aes(x, y, color = GroupContrib, shape = Block)) + 
  geom_point(size = 4, aes(text = names)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1) +
  xlab('Component 1') +
  ylab('Component 2') +
  ggtitle('Correlation Circle')  +
  #Make the circle independent of other aesthetics
  geom_path(aes(x,y), data = circle, inherit.aes = F)  +
  scale_color_manual(values = vegPal) +
  scale_y_continuous(labels = signs::signs_format()) +
  scale_x_continuous(labels = signs::signs_format()) +
  theme(legend.position = 'none')

plotly::ggplotly(corcir_72h, tooltip = 'text')


#save(diablo_3h_final, diablo_6h_final, diablo_24h_final, diablo_48h_final, diablo_72h_final, file = 'DiabloFinalModel_NoDiet.RData')

