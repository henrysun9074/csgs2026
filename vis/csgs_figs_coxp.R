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

data_dir <- "gebvs_coxp"

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
    cols = c(EN_GEBV_Mean, RR_GEBV_Mean, RF_GEBV_Mean, GB_GEBV_Mean),
    names_to = "Model",
    values_to = "GEBV_Value"
  ) %>%
  mutate(Model = case_when(
    Model == "EN_GEBV_Mean" ~ "Elastic Net (EN)",
    Model == "RR_GEBV_Mean" ~ "Ridge Regression (RR)",
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
ggsave("vis/plot_coxp.jpg", final_plot, width = 10, height = 8, dpi = 300)


# calculate pearson correlation from combined datasets for just oysters in n9 and dbw
n9_pheno <- read.csv("data/selectedlines/pheno/sel_n9_pheno.csv")
n9_survival <- n9_pheno$status_01
dbw_pheno <- read.csv("data/wildlines/pheno/wild_22dbw_pheno.csv")
dbw_survival <- dbw_pheno$status_01

# compute correlation between sel_all_gebv$LR, RF, GB_GEBV_Mean and n9_survival
## ONLY for animals with sel_all_gebv$ID in n9_pheno$SampleID
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
    Cor_RR = round(cor(RR_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_EN = round(cor(EN_GEBV_Mean, status_01, use = "complete.obs"), 4),
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
    Cor_RR = round(cor(RR_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_EN = round(cor(EN_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_RF = round(cor(RF_GEBV_Mean, status_01, use = "complete.obs"), 4),
    Cor_GB = round(cor(GB_GEBV_Mean, status_01, use = "complete.obs"), 4),
    .groups = "drop"
  )

sel_group_correlations
wild_group_correlations

ggplot(sel_group_merged, aes(x = RF_GEBV_Mean, y = status_01, color = Group)) +
  geom_point(alpha = 0.3) +
  scale_colour_brewer(palette = "Set2") +
  geom_smooth(method = "lm", se = FALSE) +  # Within-group trendlines
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", size = 1.2) + # Overall trendline
  labs(subtitle = "Within-Group Trends vs. Overall Trend",
       x = "Random Forest GEBV", y = "Survival Status") +
  theme_pubr()
# in general within group predictions are rough :(
# but the correlation accuracy is better across the overall because the groups differ in average GEBV

ggplot(wild_group_merged, aes(x = RF_GEBV_Mean, y = status_01, color = Group)) +
  geom_point(alpha = 0.3) +
  scale_colour_brewer(palette = "Set3") +
  geom_smooth(method = "lm", se = FALSE) +  # Within-group trendlines
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", size = 1.2) + # Overall trendline
  labs(subtitle = "Within-Group Trends vs. Overall Trend",
       x = "Random Forest GEBV", y = "Survival Status") +
  theme_pubr()


# genetic gain metric
compute_genetic_gain <- function(data, group_name, model_pipeline_label) {
  # Filter down to survivors only
  living_data <- data %>% filter(status_01 == 1)
  
  models <- c("RR_GEBV_Mean", "EN_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")
  
  map_df(models, function(model_col) {
    gebv_vector <- living_data[[model_col]]
    
    avg_all_living <- mean(gebv_vector, na.rm = TRUE)
    top_30_living  <- head(sort(gebv_vector, decreasing = FALSE), 30)
    avg_top_30     <- mean(top_30_living, na.rm = TRUE)
    gebv_ratio     <- abs(avg_top_30 / avg_all_living)
    
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
final_gain_summary <- final_gain_summary[order(final_gain_summary$GEBV_Ratio), decreasing = TRUE]
print(as.data.frame(final_gain_summary), row.names = FALSE)

### calculate CoxP GEBV ranks
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

wild_compare <- wild_compare %>%
  inner_join(wildall_pheno, by = c("ID" = "SampleID")) 
sel_compare <- sel_compare %>%
  inner_join(selall_pheno, by = c("ID" = "SampleID")) 

### reformat to compare between this and Status
sel_v1_single  <- sel_n9_gebv %>% filter(Dataset == "Selected (N9)")
sel_v1_pooled  <- sel_all_gebv %>% 
  filter(Dataset == "Selected (All)") %>% 
  inner_join(selall_pheno, by = c("ID" = "SampleID")) %>% 
  filter(Group == "22N9")

wild_v1_single <- wild_22dbw_gebv %>% filter(Dataset == "Wild (22dbw)")
wild_v1_pooled <- wild_all_gebv %>% 
  filter(Dataset == "Wild (All)") %>% 
  inner_join(wildall_pheno, by = c("ID" = "SampleID")) %>% 
  filter(Group == "22DBW")

sel_v1_compare <- inner_join(
  sel_v1_single %>% mutate(Rank_Sing_v1 = rank(RF_GEBV_Mean)) %>% select(ID, Rank_Sing_v1, GEBV_Sing_v1 = RF_GEBV_Mean),
  sel_v1_pooled %>% mutate(Rank_Pool_v1 = rank(RF_GEBV_Mean)) %>% select(ID, Rank_Pool_v1, GEBV_Pool_v1 = RF_GEBV_Mean, Group),
  by = "ID"
)

wild_v1_compare <- inner_join(
  wild_v1_single %>% mutate(Rank_Sing_v1 = rank(RF_GEBV_Mean)) %>% select(ID, Rank_Sing_v1, GEBV_Sing_v1 = RF_GEBV_Mean),
  wild_v1_pooled %>% mutate(Rank_Pool_v1 = rank(RF_GEBV_Mean)) %>% select(ID, Rank_Pool_v1, GEBV_Pool_v1 = RF_GEBV_Mean, Group),
  by = "ID"
)

all_v1_ranks <- rbind(wild_v1_compare, sel_v1_compare)

### compare to Status GEBVs
all_v2_ranks_raw <- readRDS("gebvs_v2/statusGEBVranks.rds")
all_v2_ranks <- all_v2_ranks_raw %>%
  select(ID, 
         Rank_Sing_v2 = Rank_Single, GEBV_Sing_v2 = GEBV_Single,
         Rank_Pool_v2 = Rank_Pooled, GEBV_Pool_v2 = GEBV_Pooled)

cross_trait_master <- inner_join(all_v1_ranks, all_v2_ranks, by = "ID")

cross_trait_metrics <- cross_trait_master %>%
  filter(Rank_Sing_v1 <= 50) %>%
  group_by(Group) %>%
  summarize(
    Single_Run_Rank_Cor = cor(Rank_Sing_v1, Rank_Sing_v2, method = "spearman"),
    Pooled_Run_Rank_Cor = cor(Rank_Pool_v1, Rank_Pool_v2, method = "spearman"),
    Avg_Single_Rank_Shift = mean(abs(Rank_Sing_v2 - Rank_Sing_v1)),
    Avg_Pooled_Rank_Shift = mean(abs(Rank_Pool_v2 - Rank_Pool_v1)),
    .groups = "drop"
  )
print(as.data.frame(cross_trait_metrics))
