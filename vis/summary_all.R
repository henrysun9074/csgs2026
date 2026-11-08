library(tidyverse)
library(ggpubr)
library(cowplot)

# Master summary across all training settings
#   group      : focal cohort the model is evaluated on (22DBW / 22N9)
#   prediction : training phenotype (Status/Cox)
#   dataset    : solo/pooled
#   model_set  : ML models (python pipeline) / R models
#
# Outputs
#   master_corr  : one row per setting x model (CV corr, focal-group sample corr, gain)
#   master_gebv  : long table of GEBVs for every animal x setting x model

# =============================================================================
# generate master tables from individual model outputs

# --- Settings table (from "Model Summary - training runs") -------------------
# R model paths are given as dir + file prefix (the "*" glob in the summary sheet).
# NOTE: sheet lists gebvs_coxp_corr/sel_22n9, but the folder on disk is sel_n9.
# settings <- tribble(
#   ~group,  ~prediction, ~dataset, ~model_set,  ~path,                                  ~prefix,
#   "22DBW", "Status",    "Solo",   "ML models", "gebvs_v2/wild_22dbw",                  NA,
#   "22DBW", "Status",    "Solo",   "R models",  "gebvs_paul/wildGEBVs",                 "22DBW",
#   "22DBW", "Status",    "Pooled", "ML models", "gebvs_v2/wild_all",                    NA,
#   "22DBW", "Status",    "Pooled", "R models",  "gebvs_paul/wildGEBVs",                 "Wild",
#   "22DBW", "Cox",       "Solo",   "ML models", "gebvs_coxp_corr/wild_22dbw",           NA,
#   "22DBW", "Cox",       "Solo",   "R models",  NA,                                     NA,
#   "22DBW", "Cox",       "Pooled", "ML models", "gebvs_coxp_corr/wild_all",             NA,
#   "22DBW", "Cox",       "Pooled", "R models",  NA,                                     NA,
#   "22N9",  "Status",    "Solo",   "ML models", "gebvs_v2/sel_n9",                      NA,
#   "22N9",  "Status",    "Solo",   "R models",  "gebvs_paul/selectedGEBVs/Status",      "22N9",
#   "22N9",  "Status",    "Pooled", "ML models", "gebvs_v2/sel_all",                     NA,
#   "22N9",  "Status",    "Pooled", "R models",  "gebvs_paul/selectedGEBVs/Status",      "SelectionLines",
#   "22N9",  "Cox",       "Solo",   "ML models", "gebvs_coxp_corr/sel_n9",               NA,
#   "22N9",  "Cox",       "Solo",   "R models",  NA,                                     NA,
#   "22N9",  "Cox",       "Pooled", "ML models", "gebvs_coxp_corr/sel_all",              NA,
#   "22N9",  "Cox",       "Pooled", "R models",  "gebvs_paul/selectedGEBVs/CoxP",        "SelectionLines"
# ) %>%
#   mutate(
#     pheno_file = case_when(
#       group == "22DBW" & dataset == "Solo"   ~ "data/wildlines/pheno/wild_22dbw_pheno.csv",
#       group == "22DBW" & dataset == "Pooled" ~ "data/wildlines/pheno/wild_all_pheno.csv",
#       group == "22N9"  & dataset == "Solo"   ~ "data/selectedlines/pheno/sel_n9_pheno.csv",
#       group == "22N9"  & dataset == "Pooled" ~ "data/selectedlines/pheno/sel_all_pheno.csv"
#     ),
#     setting_id = paste(group, prediction, dataset, model_set, sep = " | ")
#   )

# available <- settings %>% filter(!is.na(path))

# key_cols <- c("group", "prediction", "dataset", "model_set")

# # load GEBVs for ML models
# load_gebv_ml <- function(path) {
#   read.csv(file.path(path, "Final_GEBVs_Summary.csv")) %>%
#     select(ID, matches("_GEBV_(Mean|SD)$")) %>%
#     pivot_longer(-ID, names_to = c("Model", ".value"),
#                  names_pattern = "(.*)_GEBV_(Mean|SD)") %>%
#     rename(GEBV = Mean, GEBV_SD = SD)
# }

# r_trad_file <- function(path, prefix) {
#   f <- list.files(path, pattern = paste0("^", prefix, ".*_TradModels\\.csv$"), full.names = TRUE)
#   stopifnot(length(f) == 1)
#   f
# }

# load_gebv_r <- function(path, prefix) {
#   read.csv(r_trad_file(path, prefix)) %>%
#     select(ID, matches("_(Mean|SD)$")) %>%
#     pivot_longer(-ID, names_to = c("Model", ".value"),
#                  names_pattern = "(.*)_(Mean|SD)") %>%
#     rename(GEBV = Mean, GEBV_SD = SD)
# }

# # LOAD CV CORRELATIONS
# # ML: mean PearsonR over folds within each repeat, then mean over repeats
# load_cv_ml <- function(path) {
#   read.csv(file.path(path, "CrossValidation_Metrics.csv")) %>%
#     mutate(Model = gsub("_GEBV_Mean", "", Model)) %>%
#     group_by(Model, Repeat) %>%
#     summarise(Mean_Acc = mean(PearsonR, na.rm = TRUE), .groups = "drop") %>%
#     group_by(Model) %>%
#     summarise(CV_Correlation = mean(Mean_Acc, na.rm = TRUE), .groups = "drop")
# }

# # R
# load_cv_r <- function(path, prefix, prediction) {
#   target_col <- if (prediction == "Status") "Status" else "CoxP"
#   trad   <- read.csv(r_trad_file(path, prefix))
#   models <- sub("_Mean$", "", grep("_Mean$", names(trad), value = TRUE))
#   pheno  <- trad %>% select(ID, y = all_of(target_col))
#   cv_files <- list.files(path, pattern = paste0("^", prefix, ".*_CrossVal\\.csv$"), full.names = TRUE)

#   map_df(cv_files, function(f) {
#     stem  <- sub("_CrossVal\\.csv$", "", basename(f))
#     model <- models[endsWith(stem, models)]
#     stopifnot(length(model) == 1)
#     cv <- read.csv(f) %>% rename(ID = X) %>% inner_join(pheno, by = "ID")
#     gcols <- grep("^gpred(\\.[0-9]+)?$", names(cv), value = TRUE)
#     rep_cor <- map_dbl(gcols, ~ cor(cv[[.x]], cv$y, use = "complete.obs"))
#     tibble(Model = model, CV_Correlation = mean(rep_cor, na.rm = TRUE))
#   })
# }

# # Build master GEBV table
# pheno_cols <- c("SampleID", "Group", "status_01", "phenotype_gwas")

# master_gebv <- available %>%
#   pmap_df(function(group, prediction, dataset, model_set, path, prefix, pheno_file, ...) {
#     gebv <- if (model_set == "ML models") load_gebv_ml(path) else load_gebv_r(path, prefix)
#     pheno <- read.csv(pheno_file) %>% select(all_of(pheno_cols))
#     gebv %>%
#       left_join(pheno, by = c("ID" = "SampleID")) %>%
#       mutate(group = group, prediction = prediction, dataset = dataset, model_set = model_set,
#              in_focal_group = Group == group)
#   }) %>%
#   select(all_of(key_cols), Model, ID, Group, in_focal_group, GEBV, GEBV_SD, status_01, phenotype_gwas)

# # CV correlations 
# cv_all <- available %>%
#   pmap_df(function(group, prediction, dataset, model_set, path, prefix, ...) {
#     cv <- if (model_set == "ML models") load_cv_ml(path) else load_cv_r(path, prefix, prediction)
#     cv %>% mutate(group = group, prediction = prediction, dataset = dataset, model_set = model_set)
#   })

# # Focal-group sample correlation + genetic gain 
# # Sample_Cor : cor(GEBV, status_01) within focal group (as in summary_status/cox)
# # gain_top30 : mean GEBV of top 30 survivors / mean GEBV of all survivors (focal group)
# focal_metrics <- master_gebv %>%
#   filter(in_focal_group) %>%
#   group_by(across(all_of(key_cols)), Model) %>%
#   summarise(
#     n_focal    = sum(!is.na(GEBV)),
#     Sample_Cor = cor(GEBV, status_01, use = "complete.obs"),
#     gain_top30 = {
#       living <- GEBV[status_01 == 1 & !is.na(GEBV)]
#       mean(head(sort(living, decreasing = TRUE), 30)) / mean(living)
#     },
#     .groups = "drop"
#   )

# # Correlation 
# master_corr <- focal_metrics %>%
#   left_join(cv_all, by = c(key_cols, "Model")) %>%
#   mutate(
#     Correlation = if_else(dataset == "Solo", CV_Correlation, Sample_Cor),
#     across(c(CV_Correlation, Sample_Cor, Correlation, gain_top30), ~ round(.x, 4))
#   ) %>%
#   select(all_of(key_cols), Model, n_focal, Correlation, CV_Correlation, Sample_Cor, gain_top30) %>%
#   mutate(across(all_of(key_cols), as.character)) %>%
#   arrange(group, prediction, dataset, model_set, Model)

# # Settings with no results (NA paths), for reference
# missing_settings <- settings %>% filter(is.na(path)) %>% select(all_of(key_cols))

# # output
# print(as.data.frame(master_corr), row.names = FALSE)
# cat("\nSettings without results:\n")
# print(as.data.frame(missing_settings), row.names = FALSE)

# write.table(master_corr, "vis/master_correlations.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(master_gebv, "vis/master_gebvs.txt",        sep = "\t", row.names = FALSE, quote = FALSE)


### load in outputs
master_corr <- read.table("vis/master_correlations.txt", sep = "\t", header = TRUE)
master_gebv <- read.table("vis/master_correlations.txt", sep = "\t", header = TRUE)
