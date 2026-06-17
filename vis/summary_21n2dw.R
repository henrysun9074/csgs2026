library(tidyverse)

sel_data_dir <- "gebvs_21n2"
wild_data_dir <- "gebvs_21dw"

sel_21n2_status_gebv    <- read.csv(file.path(sel_data_dir, "csgs_pred2", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "21n2_status")
sel_21n2_status_metrics <- read.csv(file.path(sel_data_dir, "csgs_pred2", "CrossValidation_Metrics.csv")) %>% 
  mutate(Dataset = "21n2_status", Prediction = "status")

wild_21dw_status_gebv    <- read.csv(file.path(wild_data_dir, "csgs_pred2", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "21dw_status")
wild_21dw_status_metrics <- read.csv(file.path(wild_data_dir, "csgs_pred2", "CrossValidation_Metrics.csv")) %>% 
  mutate(Dataset = "21dw_status", Prediction = "status")

sel_21n2_cox_gebv    <- read.csv(file.path(sel_data_dir, "csgs_predCoxP_corr", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "21n2_cox")
sel_21n2_cox_metrics <- read.csv(file.path(sel_data_dir, "csgs_predCoxP_corr", "CrossValidation_Metrics.csv")) %>% 
  mutate(Dataset = "21n2_cox", Prediction = "cox")

wild_21dw_cox_gebv    <- read.csv(file.path(wild_data_dir, "csgs_predCoxP_corr", "Final_GEBVs_Summary.csv")) %>% mutate(Dataset = "21dw_cox")
wild_21dw_cox_metrics <- read.csv(file.path(wild_data_dir, "csgs_predCoxP_corr", "CrossValidation_Metrics.csv")) %>% 
  mutate(Dataset = "21dw_cox", Prediction = "cox")

# phenotype data
selall_pheno  <- read.csv("data/selectedlines/pheno/sel_all_pheno.csv")
wildall_pheno <- read.csv("data/wildlines/pheno/wild_all_pheno.csv")

metrics_all <- bind_rows(
  sel_21n2_status_metrics, sel_21n2_cox_metrics, 
  wild_21dw_status_metrics, wild_21dw_cox_metrics
)

cv_correlations <- metrics_all %>%
  group_by(Dataset, Prediction, Model, Repeat) %>%
  summarise(Mean_Acc = mean(PearsonR, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Dataset, Prediction, Model) %>%
  summarise(CV_Correlation = round(mean(Mean_Acc, na.rm = TRUE), 4), .groups = 'drop') %>%
  mutate(Model = gsub("_GEBV_Mean", "", Model))

compute_genetic_gain <- function(data, dataset_label) {
  living_data <- data %>% filter(status_01 == 1)
  possible_models <- c("LR_GEBV_Mean", "RR_GEBV_Mean", "EN_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")
  models <- intersect(possible_models, colnames(living_data))
  
  map_df(models, function(model_col) {
    gebv_vector <- living_data[[model_col]]
    
    if(length(na.omit(gebv_vector)) >= 30) {
      ratio <- mean(head(sort(gebv_vector, decreasing = TRUE), 30), na.rm = TRUE) / mean(gebv_vector, na.rm = TRUE)
    } else {
      ratio <- NA
    }
    
    tibble(
      Dataset = dataset_label, 
      Model = gsub("_GEBV_Mean", "", model_col), 
      gain_top30 = round(ratio, 4)
    )
  })
}

# Computes overall ground-truth sample correlation vs phenotypic status dynamically
get_sample_correlation <- function(df, dataset_label) {
  df %>% 
    summarise(
      across(
        any_of(c("LR_GEBV_Mean", "RR_GEBV_Mean", "EN_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")),
        ~ cor(.x, status_01, use = "complete.obs"),
        .names = "{gsub('_GEBV_Mean', '', .col)}"
      )
    ) %>%
    pivot_longer(everything(), names_to = "Model", values_to = "Sample_Cor") %>%
    mutate(Dataset = dataset_label, Sample_Cor = round(Sample_Cor, 4))
}

data_groups <- list(
  "21n2_status" = sel_21n2_status_gebv %>% inner_join(selall_pheno, by = c("ID" = "SampleID")),
  "21n2_cox"    = sel_21n2_cox_gebv    %>% inner_join(selall_pheno, by = c("ID" = "SampleID")),
  "21dw_status" = wild_21dw_status_gebv %>% inner_join(wildall_pheno, by = c("ID" = "SampleID")),
  "21dw_cox"    = wild_21dw_cox_gebv    %>% inner_join(wildall_pheno, by = c("ID" = "SampleID"))
)

build_summary_table <- function(label) {
  target_df <- data_groups[[label]]
  
  gain_df <- compute_genetic_gain(target_df, label)
  corr_df <- get_sample_correlation(target_df, label)
  
  gain_df %>%
    left_join(corr_df, by = c("Dataset", "Model")) %>%
    left_join(cv_correlations, by = c("Dataset", "Model")) %>%
    select(Dataset, Model, CV_Correlation, Sample_Cor, gain_top30)
}

# Generate tables
table_21n2_status <- build_summary_table("21n2_status")
table_21n2_cox    <- build_summary_table("21n2_cox")
table_21dw_status <- build_summary_table("21dw_status")
table_21dw_cox    <- build_summary_table("21dw_cox")

# --- Print Outputs ---
print("=== 21n2 Status Models ===")
print(as.data.frame(table_21n2_status))

print("=== 21n2 Cox Models ===")
print(as.data.frame(table_21n2_cox))

print("=== 21dw Status Models ===")
print(as.data.frame(table_21dw_status))

print("=== 21dw Cox Models ===")
print(as.data.frame(table_21dw_cox))