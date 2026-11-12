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

setwd("/work/hs325/csgs2026")

# =============================================================================
# generate master tables from individual model outputs

# RKHS only exists for 22N9 Cox Pooled but new Cox runs have BRR + BayesB only
# R 22DBW Cox Solo was trained on CoxP from wild_all_pheno (pooled Cox fit) but this has r = 0.996 to solo Cox fit

# settings <- tribble(
#   ~group,  ~prediction, ~dataset, ~model_set,  ~path,                                  ~prefix,
#   "22DBW", "Status",    "Solo",   "ML models", "gebvs_v2/wild_22dbw",                  NA,
#   "22DBW", "Status",    "Solo",   "R models",  "gebvs_paul/wildGEBVs/Status",          "22DBW",
#   "22DBW", "Status",    "Pooled", "ML models", "gebvs_v2/wild_all",                    NA,
#   "22DBW", "Status",    "Pooled", "R models",  "gebvs_paul/wildGEBVs/Status",          "Wild",
#   "22DBW", "Cox",       "Solo",   "ML models", "gebvs_coxp_corr/wild_22dbw",           NA,
#   "22DBW", "Cox",       "Solo",   "R models",  "gebvs_paul/wildGEBVs/CoxP",            "22DBW",
#   "22DBW", "Cox",       "Pooled", "ML models", "gebvs_coxp_corr/wild_all",             NA,
#   "22DBW", "Cox",       "Pooled", "R models",  "gebvs_paul/wildGEBVs/CoxP",            "Wild",
#   "22N9",  "Status",    "Solo",   "ML models", "gebvs_v2/sel_n9",                      NA,
#   "22N9",  "Status",    "Solo",   "R models",  "gebvs_paul/selectedGEBVs/Status",      "22N9",
#   "22N9",  "Status",    "Pooled", "ML models", "gebvs_v2/sel_all",                     NA,
#   "22N9",  "Status",    "Pooled", "R models",  "gebvs_paul/selectedGEBVs/Status",      "SelectionLines",
#   "22N9",  "Cox",       "Solo",   "ML models", "gebvs_coxp_corr/sel_n9",               NA,
#   "22N9",  "Cox",       "Solo",   "R models",  "gebvs_paul/selectedGEBVs/CoxP",        "22N9",
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
master_gebv <- read.table("vis/master_gebvs.txt", sep = "\t", header = TRUE)

# Summary statistics of prediction results, separately for 22DBW and 22N9
corr_tbl <- master_corr

metric_cols <- c("Correlation", "CV_Correlation", "Sample_Cor", "gain_top30")

summarise_metrics <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      n_models = n(),
      across(all_of(metric_cols),
             list(mean = ~ mean(.x, na.rm = TRUE),
                  sd   = ~ sd(.x, na.rm = TRUE),
                  min  = ~ if (all(is.na(.x))) NA_real_ else min(.x, na.rm = TRUE),
                  max  = ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric) & !n_models, ~ round(.x, 4)))
}

for (grp in c("22DBW", "22N9")) {
  grp_df <- corr_tbl %>% filter(group == grp)

  cat("\n\n#################### ", grp, " ####################\n", sep = "")
  print(as.data.frame(summarise_metrics(grp_df)), row.names = FALSE)
  print(as.data.frame(summarise_metrics(grp_df, prediction)), row.names = FALSE)
  print(as.data.frame(summarise_metrics(grp_df, dataset)), row.names = FALSE)
  print(as.data.frame(summarise_metrics(grp_df, model_set)), row.names = FALSE)

  grp_df %>%
    group_by(prediction, dataset, model_set) %>%
    summarise(n_models  = n(),
              mean_corr = round(mean(Correlation, na.rm = TRUE), 4),
              sd_corr   = round(sd(Correlation, na.rm = TRUE), 4),
              best_model = if (all(is.na(Correlation))) NA_character_ else Model[which.max(Correlation)],
              best_corr  = if (all(is.na(Correlation))) NA_real_ else max(Correlation, na.rm = TRUE),
              .groups = "drop") %>%
    as.data.frame() %>%
    print(row.names = FALSE)

}

# compare prediction accuracy within each model from pooling vs not pooling (when available)
pool_vs_solo <- master_corr %>%
  select(group, prediction, model_set, Model, dataset, Correlation, Sample_Cor) %>%
  pivot_wider(names_from = dataset, values_from = c(Correlation, Sample_Cor)) %>%
  filter(!is.na(Correlation_Solo), !is.na(Correlation_Pooled)) %>%
  mutate(
    diff_Correlation = Correlation_Pooled - Correlation_Solo,
    diff_Sample_Cor  = Sample_Cor_Pooled  - Sample_Cor_Solo,
    better_Sample_Cor = if_else(diff_Sample_Cor > 0, "Pooled", "Solo")
  ) %>%
  select(group, prediction, model_set, Model,
         Correlation_Solo, Correlation_Pooled, diff_Correlation,
         Sample_Cor_Solo, Sample_Cor_Pooled, diff_Sample_Cor, better_Sample_Cor) %>%
  arrange(group, prediction, model_set, Model)

cat("\n\n=== Pooled vs Solo, per model (diff = Pooled - Solo) ===\n")
print(as.data.frame(pool_vs_solo), row.names = FALSE)

summarise_pool_gain <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      n_pairs               = n(),
      n_pooled_better       = sum(diff_Sample_Cor > 0),
      mean_diff_Sample_Cor  = round(mean(diff_Sample_Cor), 4),
      median_diff_Sample_Cor = round(median(diff_Sample_Cor), 4),
      mean_diff_Correlation = round(mean(diff_Correlation), 4),
      # paired Wilcoxon signed-rank test
      wilcox_p = if (n() >= 3) round(wilcox.test(Sample_Cor_Pooled, Sample_Cor_Solo,
                                                 paired = TRUE, exact = FALSE)$p.value, 4)
                 else NA_real_,
      .groups = "drop"
    )
}

print(as.data.frame(summarise_pool_gain(pool_vs_solo, group)), row.names = FALSE)
print(as.data.frame(summarise_pool_gain(pool_vs_solo, group, prediction, model_set)), row.names = FALSE)

################ save outputs
write_csv(corr_tbl, "data/corr_tbl_8nov.csv")
