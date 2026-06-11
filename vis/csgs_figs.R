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


### do this for all other groups 
selall_pheno <- read.csv("data/selectedlines/pheno/sel_all_pheno.csv")
selall_survival <- selall_pheno$status_01
wildall_pheno <- read.csv("data/wildlines/pheno/wild_all_pheno.csv")
wildall_survival <- wildall_pheno$status_01

sel_group_merged <- inner_join(
  sel_all_gebv, 
  selall_pheno, 
  by = c("ID" = "SampleID")
)
sel_group_correlations <- sel_group_merged %>%
  group_by(Group) %>%
  summarize(
    N_Animals = n(),
    Cor_LR = round(cor(LR_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_RF = round(cor(RF_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_GB = round(cor(GB_GEBV_Mean, status_01, use = "complete.obs"), 4),
    .groups = "drop"
  )

wild_group_merged <- inner_join(
  wild_all_gebv, 
  wildall_pheno, 
  by = c("ID" = "SampleID")
)
wild_group_correlations <- wild_group_merged %>%
  group_by(Group) %>%
  summarize(
    N_Animals = n(),
    Cor_LR = round(cor(LR_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_RF = round(cor(RF_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_GB = round(cor(GB_GEBV_Mean, status_01, use = "complete.obs"), 4),
    .groups = "drop"
  )

ggplot(sel_group_merged, aes(x = LR_GEBV_Mean, y = status_01, color = Group)) +
  geom_point(alpha = 0.3) +
  scale_colour_brewer(palette = "Set2") +
  geom_smooth(method = "lm", se = FALSE) +  # Within-group trendlines
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", size = 1.2) + # Overall trendline
  labs(subtitle = "Within-Group Trends vs. Overall Trend",
       x = "Random Forest GEBV", y = "Survival Status") +
  theme_pubr()
# in general within group predictions are rough :(
# but the correlation accuracy is better across the overall because the groups differ in average GEBV

ggplot(wild_group_merged, aes(x = LR_GEBV_Mean, y = status_01, color = Group)) +
  geom_point(alpha = 0.3) +
  scale_colour_brewer(palette = "Set3") +
  geom_smooth(method = "lm", se = FALSE) +  # Within-group trendlines
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", size = 1.2) + # Overall trendline
  labs(subtitle = "Within-Group Trends vs. Overall Trend",
       x = "LR GEBV", y = "Survival Status") +
  theme_pubr()


### compute genetic gain simple function
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


### separate live vs dead in density plots
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

## but how different are GEBV rankings between single vs pooled run
sel_single_22n9 <- sel_n9_gebv %>% filter(Dataset == "Selected (N9)")
sel_pooled_22n9 <- sel_all_gebv %>%
  filter(Dataset == "Selected (All)") %>%
  inner_join(selall_pheno, by = c("ID" = "SampleID")) %>%
  filter(Group == "22N9") 

rank_single <- sel_single_22n9 %>%
  mutate(Rank_Single = rank(desc(RF_GEBV_Mean))) %>%
  select(ID, Rank_Single)

rank_pooled <- sel_pooled_22n9 %>%
  mutate(Rank_Pooled = rank(desc(RF_GEBV_Mean))) %>%
  select(ID, Rank_Pooled)

sel_22n9_compare <- inner_join(rank_single, rank_pooled, by = "ID")
top_50_ids <- sel_22n9_compare %>%
  filter(Rank_Single <= 50 | Rank_Pooled <= 50) %>%
  pull(ID)

heatmap_data_sel <- sel_22n9_compare %>%
  filter(ID %in% top_50_ids) %>%
  pivot_longer(cols = c(Rank_Single, Rank_Pooled), 
               names_to = "Run_Type", 
               values_to = "Rank") %>%
  mutate(Run_Type = recode(Run_Type, "Rank_Single" = "Single Run (N9)", "Rank_Pooled" = "Pooled Run (Within-Cohort)"))

heatmap_data_sel$ID <- factor(heatmap_data_sel$ID, 
                              levels = sel_22n9_compare %>% 
                                filter(ID %in% top_50_ids) %>% 
                                arrange(Rank_Single) %>% 
                                pull(ID))
ggplot(heatmap_data_sel, aes(x = Run_Type, y = ID, fill = Rank)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_viridis_c(direction = -1, option = "mako", name = "Cohort Rank") +
  theme_minimal() +
  labs(title = "GEBV Ranks (22N9)",
       x = "Training Data", y = "Oyster ID") +
  theme(axis.text.y = element_text(size = 6), panel.grid = element_blank())

## average difference in gebvs and average rank shift for top 50
sel_single_22n9 <- sel_n9_gebv %>% filter(Dataset == "Selected (N9)")
sel_pooled_22n9 <- sel_all_gebv %>%
  filter(Dataset == "Selected (All)") %>%
  inner_join(selall_pheno, by = c("ID" = "SampleID")) %>%
  filter(Group == "22N9")

sel_compare <- inner_join(
  sel_single_22n9 %>% mutate(Rank_Single = rank(desc(RF_GEBV_Mean))) %>% select(ID, Rank_Single, GEBV_Single = RF_GEBV_Mean),
  sel_pooled_22n9 %>% mutate(Rank_Pooled = rank(desc(RF_GEBV_Mean))) %>% select(ID, Rank_Pooled, GEBV_Pooled = RF_GEBV_Mean),
  by = "ID"
)
sel_top50_metrics <- sel_compare %>%
  filter(Rank_Single <= 50) %>%
  mutate(
    GEBV_Raw_Diff = GEBV_Pooled - GEBV_Single,
    GEBV_Abs_Diff = abs(GEBV_Pooled - GEBV_Single),
    Rank_Shift    = abs(Rank_Pooled - Rank_Single)
  ) %>%
  summarize(
    Cohort = "Selected (22N9) - Top 50 Only",
    Oysters_Evaluated = n(),
    Avg_Raw_GEBV_Diff = mean(GEBV_Raw_Diff),
    Avg_Abs_GEBV_Diff = mean(GEBV_Abs_Diff),
    Avg_Rank_Shift    = mean(Rank_Shift),
    Max_Rank_Shift    = max(Rank_Shift)
  )
sel_top50_metrics

## but how different are GEBV rankings between single vs pooled run for wild
wild_single_dbw <- wild_22dbw_gebv %>% filter(Dataset == "Wild (22dbw)")

wild_pooled_dbw <- wild_all_gebv %>%
  filter(Dataset == "Wild (All)") %>%
  inner_join(wildall_pheno, by = c("ID" = "SampleID")) %>%
  filter(Group == "22DBW") 

rank_single_wild <- wild_single_dbw %>%
  mutate(Rank_Single = rank(desc(RF_GEBV_Mean))) %>%
  select(ID, Rank_Single)

rank_pooled_wild <- wild_pooled_dbw %>%
  mutate(Rank_Pooled = rank(desc(RF_GEBV_Mean))) %>%
  select(ID, Rank_Pooled)

wild_dbw_compare <- inner_join(rank_single_wild, rank_pooled_wild, by = "ID")
top_50_wild_ids <- wild_dbw_compare %>%
  filter(Rank_Single <= 50 | Rank_Pooled <= 50) %>%
  pull(ID)

heatmap_data_wild <- wild_dbw_compare %>%
  filter(ID %in% top_50_wild_ids) %>%
  pivot_longer(cols = c(Rank_Single, Rank_Pooled), 
               names_to = "Run_Type", 
               values_to = "Rank") %>%
  mutate(Run_Type = recode(Run_Type, "Rank_Single" = "Single Run (22DBW)", "Rank_Pooled" = "Pooled Run (Within-Cohort)"))

heatmap_data_wild$ID <- factor(heatmap_data_wild$ID, 
                               levels = wild_dbw_compare %>% 
                                 filter(ID %in% top_50_wild_ids) %>% 
                                 arrange(Rank_Single) %>% 
                                 pull(ID))

ggplot(heatmap_data_wild, aes(x = Run_Type, y = ID, fill = Rank)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_viridis_c(direction = -1, option = "plasma", name = "Cohort Rank") +
  theme_minimal() +
  labs(title = "GEBV Ranks (22DBW)",
       x = "Model Framework", y = "Oyster ID") +
  theme(axis.text.y = element_text(size = 6), panel.grid = element_blank())

wild_single_dbw <- wild_22dbw_gebv %>% filter(Dataset == "Wild (22dbw)")
wild_pooled_dbw <- wild_all_gebv %>%
  filter(Dataset == "Wild (All)") %>%
  inner_join(wildall_pheno, by = c("ID" = "SampleID")) %>%
  filter(Group == "22DBW")

wild_compare <- inner_join(
  wild_single_dbw %>% mutate(Rank_Single = rank(desc(RF_GEBV_Mean))) %>% select(ID, Rank_Single, GEBV_Single = RF_GEBV_Mean),
  wild_pooled_dbw %>% mutate(Rank_Pooled = rank(desc(RF_GEBV_Mean))) %>% select(ID, Rank_Pooled, GEBV_Pooled = RF_GEBV_Mean),
  by = "ID"
)

wild_top50_metrics <- wild_compare %>%
  filter(Rank_Single <= 50) %>%
  mutate(
    GEBV_Raw_Diff = GEBV_Pooled - GEBV_Single,
    GEBV_Abs_Diff = abs(GEBV_Pooled - GEBV_Single),
    Rank_Shift    = abs(Rank_Pooled - Rank_Single)
  ) %>%
  summarize(
    Cohort = "Wild (22DBW) - Top 50 Only",
    Oysters_Evaluated = n(),
    Avg_Raw_GEBV_Diff = mean(GEBV_Raw_Diff),
    Avg_Abs_GEBV_Diff = mean(GEBV_Abs_Diff),
    Avg_Rank_Shift    = mean(Rank_Shift),
    Max_Rank_Shift    = max(Rank_Shift)
  )

all_gebv_metrics <- rbind(wild_top50_metrics, sel_top50_metrics)
all_gebv_metrics$Avg_Raw_GEBV_Diff <- NULL
all_gebv_metrics$Oysters_Evaluated <- NULL
all_gebv_metrics
