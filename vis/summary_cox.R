library(tidyverse)

data_dir <- "gebvs_coxp"

# --- Load Data ---
sel_n9_gebv     <- read.csv(file.path(data_dir, "sel_n9", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "Selected (N9)")
sel_n9_metrics  <- read.csv(file.path(data_dir, "sel_n9", "CrossValidation_Metrics.csv")) %>% mutate(Dataset = "Selected (N9)")
sel_all_gebv    <- read.csv(file.path(data_dir, "sel_all", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "Selected (All)")
sel_all_metrics <- read.csv(file.path(data_dir, "sel_all", "CrossValidation_Metrics.csv")) %>% mutate(Dataset = "Selected (All)")

wild_22dbw_gebv <- read.csv(file.path(data_dir, "wild_22dbw", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "Wild (22dbw)")
wild_22dbw_metrics <- read.csv(file.path(data_dir, "wild_22dbw", "CrossValidation_Metrics.csv")) %>% mutate(Dataset = "Wild (22dbw)")
wild_all_gebv   <- read.csv(file.path(data_dir, "wild_all", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "Wild (All)")
wild_all_metrics <- read.csv(file.path(data_dir, "wild_all", "CrossValidation_Metrics.csv")) %>% mutate(Dataset = "Wild (All)")

n9_pheno      <- read.csv("data/selectedlines/pheno/sel_n9_pheno.csv")
dbw_pheno     <- read.csv("data/wildlines/pheno/wild_22dbw_pheno.csv")
selall_pheno  <- read.csv("data/selectedlines/pheno/sel_all_pheno.csv")
wildall_pheno <- read.csv("data/wildlines/pheno/wild_all_pheno.csv")

metrics_all <- bind_rows(sel_n9_metrics, sel_all_metrics, wild_22dbw_metrics, wild_all_metrics)
cv_correlations <- metrics_all %>%
  group_by(Dataset, Model, Repeat) %>%
  summarise(Mean_Acc = mean(PearsonR, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Dataset, Model) %>%
  summarise(Correlation = round(mean(Mean_Acc, na.rm = TRUE), 4), .groups = 'drop') %>%
  mutate(Model = gsub("_GEBV_Mean", "", Model))

# --- Genetic Gain Function (Cox Specific Models) ---
compute_genetic_gain_cox <- function(data, group_name, model_pipeline_label) {
  living_data <- data %>% filter(status_01 == 1)
  models <- c("RR_GEBV_Mean", "EN_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")
  map_df(models, function(model_col) {
    gebv_vector <- living_data[[model_col]]
    ratio <- mean(head(sort(gebv_vector, decreasing = TRUE), 30), na.rm = TRUE) / mean(gebv_vector, na.rm = TRUE)
    tibble(Cohort = group_name, Pipeline = model_pipeline_label, 
           Model = gsub("_GEBV_Mean", "", model_col), gain_top30 = round(ratio, 4))
  })
}

# --- Prepare Datasets ---
sel_single_data  <- sel_n9_gebv %>% inner_join(n9_pheno, by = c("ID" = "SampleID"))
sel_pooled_data  <- sel_all_gebv %>% inner_join(selall_pheno, by = c("ID" = "SampleID")) %>% filter(Group == "22N9")
wild_single_data <- wild_22dbw_gebv %>% inner_join(dbw_pheno, by = c("ID" = "SampleID"))
wild_pooled_data <- wild_all_gebv %>% inner_join(wildall_pheno, by = c("ID" = "SampleID")) %>% filter(Group == "22DBW")

# --- Manual Pooled Correlations (Cox Specific Models) ---
get_pooled_cor_cox <- function(df) {
  df %>% summarise(RR = cor(RR_GEBV_Mean, status_01, use = "complete.obs"),
                   EN = cor(EN_GEBV_Mean, status_01, use = "complete.obs"),
                   RF = cor(RF_GEBV_Mean, status_01, use = "complete.obs"),
                   GB = cor(GB_GEBV_Mean, status_01, use = "complete.obs")) %>%
    pivot_longer(everything(), names_to = "Model", values_to = "Pooled_Cor")
}

sel_p_cor_cox <- get_pooled_cor_cox(sel_pooled_data)
wild_p_cor_cox <- get_pooled_cor_cox(wild_pooled_data)

# --- Final Tables ---
cox_selected_table <- bind_rows(compute_genetic_gain_cox(sel_single_data, "22N9", "Single"), 
                                compute_genetic_gain_cox(sel_pooled_data, "22N9", "Pooled")) %>%
  mutate(CV_Match = if_else(Pipeline == "Single", "Selected (N9)", "Selected (All)")) %>%
  left_join(cv_correlations, by = c("Model", "CV_Match" = "Dataset")) %>%
  left_join(sel_p_cor_cox, by = "Model") %>%
  mutate(Correlation = if_else(Pipeline == "Pooled", round(Pooled_Cor, 4), Correlation),
         Dataset = if_else(Pipeline == "Single", "single n9", "pooled")) %>%
  select(Dataset, Model, Correlation, gain_top30)

cox_wild_table <- bind_rows(compute_genetic_gain_cox(wild_single_data, "22DBW", "Single"), 
                            compute_genetic_gain_cox(wild_pooled_data, "22DBW", "Pooled")) %>%
  mutate(CV_Match = if_else(Pipeline == "Single", "Wild (22dbw)", "Wild (All)")) %>%
  left_join(cv_correlations, by = c("Model", "CV_Match" = "Dataset")) %>%
  left_join(wild_p_cor_cox, by = "Model") %>%
  mutate(Correlation = if_else(Pipeline == "Pooled", round(Pooled_Cor, 4), Correlation),
         Dataset = if_else(Pipeline == "Single", "single dbw", "pooled")) %>%
  select(Dataset, Model, Correlation, gain_top30)

print(as.data.frame(cox_selected_table))
print(as.data.frame(cox_wild_table))