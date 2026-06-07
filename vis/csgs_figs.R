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


### plot across all folds
p1 <- ggplot(metrics_all, aes(x = Model, y = PearsonR, fill = Model)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.2) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3.5, color = "darkturquoise") +
  facet_wrap(~Dataset, scales = "free_y") +
  labs(x = "Model", y = "Correlation Accuracy") +
  theme_pubr() + 
  theme(legend.position = "none")
p1

### plot average from each iteration of 5fold CV (more stable)
repeat_summary <- metrics_all %>%
  group_by(Dataset, Model, Repeat) %>%
  summarise(Mean_Accuracy = mean(PearsonR, na.rm = TRUE), .groups = 'drop')
p2 <- ggplot(repeat_summary, aes(x = Model, y = Mean_Accuracy, fill = Model)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~Dataset, scales = "free_y") +
  labs(x = "Model", y = "Correlation Accuracy") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
  theme_pubr() +
  theme(legend.position = "none")
p2

# make density plot
gebv_long <- gebv_all %>%
  pivot_longer(
    cols = c(LR_GEBV_Mean, RF_GEBV_Mean, GB_GEBV_Mean),
    names_to = "Model",
    values_to = "GEBV_Value"
  ) %>%
  mutate(Model = case_when(
    Model == "LR_GEBV_Mean" ~ "Logistic Regression (LR)",
    Model == "RF_GEBV_Mean" ~ "Random Forest (RF)",
    Model == "GB_GEBV_Mean" ~ "Gradient Boosting (GB)",
    TRUE ~ Model
  ))

p3_updated <- ggplot(gebv_long, aes(x = GEBV_Value, fill = Dataset)) +
  geom_density(alpha = 0.6) +
  facet_grid(Dataset ~ Model, scales = "free") +
  labs(
    x = "GEBV", 
    y = "Density"
  ) +
  theme_pubr() +
  scale_fill_jco() +
  theme(
    legend.position = "right",
    panel.spacing.x = unit(1.5, "lines"), 
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_blank() 
  )
p3_updated

# combine plots
final_plot <- (p1 | p2) / p3_updated
final_plot + plot_layout(guides = 'collect')
ggsave("vis/plot.jpg", final_plot, width = 10, height = 8, dpi = 300)

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
wild_cor_lr <- cor(wild_merged$LR_GEBV_Mean, wild_merged$status_01, use = "complete.obs")
wild_cor_rf <- cor(wild_merged$RF_GEBV_Mean, wild_merged$status_01, use = "complete.obs")
wild_cor_gb <- cor(wild_merged$GB_GEBV_Mean, wild_merged$status_01, use = "complete.obs")

cat("=== Selected Lines (All GEBV vs N9 Survival) ===\n")
cat(paste("Linear Regression (LR):", round(sel_cor_lr, 4), "\n"))
cat(paste("Random Forest (RF):    ", round(sel_cor_rf, 4), "\n"))
cat(paste("Gradient Boosting (GB):", round(sel_cor_gb, 4), "\n\n"))

cat("=== Wild Lines (All GEBV vs DBW Survival) ===\n")
cat(paste("Linear Regression (LR):", round(wild_cor_lr, 4), "\n"))
cat(paste("Random Forest (RF):    ", round(wild_cor_rf, 4), "\n"))
cat(paste("Gradient Boosting (GB):", round(wild_cor_gb, 4), "\n"))


### compute genetic gain simple function

# for n9 gebv only
sel_merged <- inner_join(
  sel_n9_gebv, 
  n9_pheno, 
  by = c("ID" = "SampleID")
)

living_n9 <- sel_merged %>% 
  filter(status_01 == 1)

models <- c("LR_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")
results <- map_df(models, function(model_col) {
  gebv_vector <- living_n9[[model_col]]
    avg_all_living <- mean(gebv_vector, na.rm = TRUE)
  top_30_living <- head(sort(gebv_vector, decreasing = TRUE), 30)
  avg_top_30 <- mean(top_30_living, na.rm = TRUE)
  gebv_ratio <- avg_top_30 / avg_all_living
  tibble(
    Model = gsub("_GEBV_Mean", "", model_col),
    Avg_Top_30_Living = avg_top_30,
    Avg_All_Living = avg_all_living,
    GEBV_Ratio = gebv_ratio
  )
})
print(results)


### compute genetic gain simple function

# for n9 gebv only
wild_merged <- inner_join(
  wild_22dbw_gebv, 
  dbw_pheno, 
  by = c("ID" = "SampleID")
)

living_dbw <- wild_merged %>% 
  filter(status_01 == 1)

models <- c("LR_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")
results <- map_df(models, function(model_col) {
  gebv_vector <- living_dbw[[model_col]]
  avg_all_living <- mean(gebv_vector, na.rm = TRUE)
  top_30_living <- head(sort(gebv_vector, decreasing = TRUE), 30)
  avg_top_30 <- mean(top_30_living, na.rm = TRUE)
  gebv_ratio <- avg_top_30 / avg_all_living
  tibble(
    Model = gsub("_GEBV_Mean", "", model_col),
    Avg_Top_30_Living = avg_top_30,
    Avg_All_Living = avg_all_living,
    GEBV_Ratio = gebv_ratio
  )
})
print(results)


### separate live vs dead
pheno_lookup <- bind_rows(
  n9_pheno  %>% select(SampleID, status_01),
  dbw_pheno %>% select(SampleID, status_01)
) %>% 
  distinct(SampleID, .keep_all = TRUE) # Ensure unique IDs

gebv_long_status <- gebv_long %>%
  left_join(pheno_lookup, by = c("ID" = "SampleID")) %>%
  filter(!is.na(status_01)) %>% 
  mutate(Survival_Status = case_when(
    status_01 == 0 ~ "Dead (0)",
    status_01 == 1 ~ "Alive (1)",
    TRUE ~ as.character(status_01)
  ))


#### plotting to separate live and dead
p4 <- ggplot(gebv_long_status, aes(x = GEBV_Value, fill = Survival_Status)) +
  geom_density(alpha = 0.3) + 
  facet_grid(Dataset ~ Model, scales = "free") +
  labs(
    x = "GEBV", 
    y = "Density",
    fill = "Survival Status"
  ) +
  theme_pubr() +
  scale_fill_manual(values = c("Dead (0)" = "#E41A1C", "Alive (1)" = "#377EB8")) + 
  theme(
    legend.position = "right",
    panel.spacing.x = unit(1.5, "lines"), 
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 10) 
  )
p4
ggsave("vis/plot2.jpg", p4, width = 10, height = 8, dpi = 300)