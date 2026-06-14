library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggsignif)
library(multcompView)
library(ggsci)
library(cowplot)
library(corrplot)
library(ggridges)
library(patchwork)

data_dir <- "gebvs_v2"

# --- Selected Lines ---

sel_n9_gebv    <- read.csv(file.path(data_dir, "sel_n9", "Final_GEBVs_Summary.csv"))
sel_n9_metrics <- read.csv(file.path(data_dir, "sel_n9", "CrossValidation_Metrics.csv"))

sel_all_gebv    <- read.csv(file.path(data_dir, "sel_all", "Final_GEBVs_Summary.csv"))
sel_all_metrics <- read.csv(file.path(data_dir, "sel_all", "CrossValidation_Metrics.csv"))

# --- Wild Lines ---
wild_22dbw_gebv    <- read.csv(file.path(data_dir, "wild_22dbw", "Final_GEBVs_Summary.csv"))
wild_22dbw_metrics <- read.csv(file.path(data_dir, "wild_22dbw", "CrossValidation_Metrics.csv"))

wild_all_gebv    <- read.csv(file.path(data_dir, "wild_all", "Final_GEBVs_Summary.csv"))
wild_all_metrics <- read.csv(file.path(data_dir, "wild_all", "CrossValidation_Metrics.csv"))

# Add dataset for merging
sel_n9_metrics     <- sel_n9_metrics     %>% mutate(Dataset = "Selected (N9)")
sel_all_metrics    <- sel_all_metrics    %>% mutate(Dataset = "Selected (All)")
wild_22dbw_metrics <- wild_22dbw_metrics %>% mutate(Dataset = "Wild (22dbw)")
wild_all_metrics   <- wild_all_metrics   %>% mutate(Dataset = "Wild (All)")

sel_n9_gebv     <- sel_n9_gebv     %>% mutate(Dataset = "Selected (N9)")
sel_all_gebv    <- sel_all_gebv    %>% mutate(Dataset = "Selected (All)")
wild_22dbw_gebv <- wild_22dbw_gebv %>% mutate(Dataset = "Wild (22dbw)")
wild_all_gebv   <- wild_all_gebv   %>% mutate(Dataset = "Wild (All)")

metrics_all <- bind_rows(sel_n9_metrics, sel_all_metrics, wild_22dbw_metrics, wild_all_metrics)
gebv_all <- bind_rows(sel_n9_gebv, sel_all_gebv, wild_22dbw_gebv, wild_all_gebv)

### plot average from each iteration of 5fold CV (more stable)
repeat_summary <- metrics_all %>%
  group_by(Dataset, Model, Repeat) %>%
  summarise(Mean_Accuracy = mean(PearsonR, na.rm = TRUE), .groups = 'drop')

cv_correlations <- repeat_summary %>%
  group_by(Dataset, Model) %>%
  summarise(Correlation = round(mean(Mean_Accuracy, na.rm = TRUE), 4), .groups = 'drop')

# Standardize the Model names to match your genetic gain output (LR, RF, GB)
cv_correlations <- cv_correlations %>%
  mutate(Model = case_when(
    Model == "LR_GEBV_Mean" ~ "LR",
    Model == "RF_GEBV_Mean" ~ "RF",
    Model == "GB_GEBV_Mean" ~ "GB",
    TRUE ~ Model
  ))

# calculate pearson correlation from combined datasets for just oysters in n9 and dbw
n9_pheno <- read.csv("data/selectedlines/pheno/sel_n9_pheno.csv")
n9_survival <- n9_pheno$status_01
dbw_pheno <- read.csv("data/wildlines/pheno/wild_22dbw_pheno.csv")
dbw_survival <- dbw_pheno$status_01

# compute correlation between sel_all_gebv$LR, RF, GB_GEBV_Mean and n9_survival
## ONLY for animals with sel_all_gebv$ID in n9_pheno$SampleID

sel_merged <- inner_join(
  sel_all_gebv,
  n9_pheno,
  by = c("ID" = "SampleID")
)
sel_cor_lr <- cor(sel_merged$LR_GEBV_Mean, sel_merged$status_01, use = "complete.obs")
sel_cor_rf <- cor(sel_merged$RF_GEBV_Mean, sel_merged$status_01, use = "complete.obs")
sel_cor_gb <- cor(sel_merged$GB_GEBV_Mean, sel_merged$status_01, use = "complete.obs")

wild_merged <- inner_join(
  wild_all_gebv,
  dbw_pheno,
  by = c("ID" = "SampleID")
)

selall_pheno  <- read.csv("data/selectedlines/pheno/sel_all_pheno.csv")
wildall_pheno <- read.csv("data/wildlines/pheno/wild_all_pheno.csv")

sel_pooled_data <- sel_all_gebv %>% 
  filter(Dataset == "Wild (All)") %>% 
  inner_join(selall_pheno, by = c("ID" = "SampleID")) %>% 
  filter(Group == "22N9")

wild_pooled_data <- wild_all_gebv %>% 
  filter(Dataset == "W (All)") %>% 
  inner_join(selall_pheno, by = c("ID" = "SampleID")) %>% 
  filter(Group == "22DBW")

compute_genetic_gain <- function(data, group_name, model_pipeline_label) {
  # Filter down to survivors only
  living_data <- data %>% filter(status_01 == 1)
  models <- c("LR_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")

  map_df(models, function(model_col) {
    gebv_vector <- living_data[[model_col]]
    avg_all_living <- mean(gebv_vector, na.rm = TRUE)
    top_30_living  <- head(sort(gebv_vector, decreasing = TRUE), 30)
    avg_top_30     <- mean(top_30_living, na.rm = TRUE)
    gebv_ratio     <- avg_top_30 / avg_all_living
    tibble(
      Cohort       = group_name,
      Pipeline     = model_pipeline_label,
      Model        = gsub("_GEBV_Mean", "", model_col),
      Avg_All      = round(avg_all_living, 4),
      Avg_Top_30   = round(avg_top_30, 4),
      GEBV_Ratio   = round(gebv_ratio, 4)
    )
  })
}

sel_single_data <- sel_n9_gebv %>%
  filter(Dataset == "Selected (N9)") %>%
  inner_join(n9_pheno, by = c("ID" = "SampleID"))

sel_pooled_data <- sel_all_gebv %>%
  filter(Dataset == "Selected (All)") %>%
  inner_join(selall_pheno, by = c("ID" = "SampleID")) %>%
  filter(Group == "22N9")

sel_single_gain <- compute_genetic_gain(sel_single_data, "Selected (22N9)", "Single")
sel_pooled_gain <- compute_genetic_gain(sel_pooled_data, "Selected (22N9)", "Pooled")

wild_single_data <- wild_22dbw_gebv %>%
  filter(Dataset == "Wild (22dbw)") %>%
  inner_join(dbw_pheno, by = c("ID" = "SampleID"))

wild_pooled_data <- wild_all_gebv %>%
  filter(Dataset == "Wild (All)") %>%
  inner_join(wildall_pheno, by = c("ID" = "SampleID")) %>%
  filter(Group == "22DBW")

wild_single_gain <- compute_genetic_gain(wild_single_data, "Wild (22DBW)", "Single")
wild_pooled_gain <- compute_genetic_gain(wild_pooled_data, "Wild (22DBW)", "Pooled")

final_gain_summary <- bind_rows(
  sel_single_gain, sel_pooled_gain,
  wild_single_gain, wild_pooled_gain
)
print(as.data.frame(final_gain_summary), row.names = FALSE)

#########################################################
sel_pooled_cor <- sel_pooled_data %>%
  summarise(
    LR = cor(LR_GEBV_Mean, status_01, use = "complete.obs"),
    RF = cor(RF_GEBV_Mean, status_01, use = "complete.obs"),
    GB = cor(GB_GEBV_Mean, status_01, use = "complete.obs")
  ) %>%
  pivot_longer(cols = everything(), names_to = "Model", values_to = "Pooled_Cor")

wild_pooled_cor <- wild_pooled_data %>%
  summarise(
    LR = cor(LR_GEBV_Mean, status_01, use = "complete.obs"),
    RF = cor(RF_GEBV_Mean, status_01, use = "complete.obs"),
    GB = cor(GB_GEBV_Mean, status_01, use = "complete.obs")
  ) %>%
  pivot_longer(cols = everything(), names_to = "Model", values_to = "Pooled_Cor")

selected_table <- bind_rows(sel_single_gain, sel_pooled_gain) %>%
  mutate(CV_Dataset_Match = if_else(Pipeline == "Single", "Selected (N9)", "Selected (All)")) %>%
  left_join(cv_correlations, by = c("CV_Dataset_Match" = "Dataset", "Model" = "Model")) %>%
  left_join(sel_pooled_cor, by = "Model") %>%
  mutate(Correlation = if_else(Pipeline == "Pooled", round(Pooled_Cor, 4), Correlation)) %>%
  select(
    Dataset = Pipeline,
    Model,
    Correlation,
    gain_top30 = GEBV_Ratio
  ) %>%
  mutate(Dataset = if_else(Dataset == "Single", "single n9", "pooled"))

cat("\n=== SELECTED LINES SUMMARY TABLE ===\n")
cat("Note: Pooled correlation is calculated only for Group 22N9\n")
print(as.data.frame(selected_table), row.names = FALSE)

wild_table <- bind_rows(wild_single_gain, wild_pooled_gain) %>%
  mutate(CV_Dataset_Match = if_else(Pipeline == "Single", "Wild (22dbw)", "Wild (All)")) %>%
  # Join CV correlations for Single
  left_join(cv_correlations, by = c("CV_Dataset_Match" = "Dataset", "Model" = "Model")) %>%
  # Join the specific Pooled correlations calculated above
  left_join(wild_pooled_cor, by = "Model") %>%
  mutate(Correlation = if_else(Pipeline == "Pooled", round(Pooled_Cor, 4), Correlation)) %>%
  select(
    Dataset = Pipeline,
    Model,
    Correlation,
    gain_top30 = GEBV_Ratio
  ) %>%
  mutate(Dataset = if_else(Dataset == "Single", "single dbw", "pooled"))

print(as.data.frame(wild_table), row.names = FALSE)
