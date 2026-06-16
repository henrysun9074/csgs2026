library(tidyverse)

# Status Phenotype
sel_all_gebv_status  <- read.csv("gebvs_v2/sel_all/Final_GEBVs_Summary.csv")
wild_all_gebv_status <- read.csv("gebvs_v2/wild_all/Final_GEBVs_Summary.csv")

# Cox Phenotype
sel_all_gebv_cox     <- read.csv("gebvs_coxp/sel_all/Final_GEBVs_Summary.csv")
wild_all_gebv_cox    <- read.csv("gebvs_coxp/wild_all/Final_GEBVs_Summary.csv")

# Survival data
selall_pheno         <- read.csv("data/selectedlines/pheno/sel_all_pheno.csv")
wildall_pheno        <- read.csv("data/wildlines/pheno/wild_all_pheno.csv")

compute_pooled_group_metrics <- function(gebv_df, pheno_df, pred_type, model_cols) {
  
  # Clean up model naming schema (e.g., "RF_GEBV_Mean" -> "RF")
  clean_names <- gsub("_GEBV_Mean", "", model_cols)
  
  # Merge GEBVs with Phenotypes to get the "Group" column for every animal
  merged_data <- gebv_df %>% 
    inner_join(pheno_df, by = c("ID" = "SampleID"))
  
  # Loop through every unique group found in the phenotype file
  map_df(unique(merged_data$Group), function(grp) {
    grp_data <- merged_data %>% filter(Group == grp)
    
    # Isolate survivors within this group for the genetic gain metric
    living_grp_data <- grp_data %>% filter(status_01 == 1)
    
    # Calculate metrics for each model architecture
    map2_df(model_cols, clean_names, function(m_col, m_name) {
      
      # PearsonR
      pearson_r <- cor(grp_data[[m_col]], grp_data$status_01, use = "complete.obs")
      
      # Genetic gain
      gebv_vector <- living_grp_data[[m_col]]
      
      # Can adjust count 
      n_select <- min(30, length(gebv_vector))
      
      gain <- if(n_select > 0) {
        mean(head(sort(gebv_vector, decreasing = TRUE), n_select), na.rm = TRUE) / mean(gebv_vector, na.rm = TRUE)
      } else {
        NA
      }
      
      tibble(
        prediction_variable = pred_type,
        model               = m_name,
        group               = grp,
        pearsonR            = round(pearson_r, 4),
        gain                = round(gain, 4)
      )
    })
  })
}

status_models <- c("LR_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")
cox_models    <- c("RR_GEBV_Mean", "EN_GEBV_Mean", "RF_GEBV_Mean", "GB_GEBV_Mean")

sel_status_metrics <- compute_pooled_group_metrics(sel_all_gebv_status, selall_pheno, "Status", status_models)
sel_cox_metrics    <- compute_pooled_group_metrics(sel_all_gebv_cox,    selall_pheno, "Cox",    cox_models)

wild_status_metrics <- compute_pooled_group_metrics(wild_all_gebv_status, wildall_pheno, "Status", status_models)
wild_cox_metrics    <- compute_pooled_group_metrics(wild_all_gebv_cox,    wildall_pheno, "Cox",    cox_models)

master_selected_table <- bind_rows(sel_status_metrics, sel_cox_metrics) %>% arrange(group, prediction_variable, model)
master_wild_table     <- bind_rows(wild_status_metrics, wild_cox_metrics) %>% arrange(group, prediction_variable, model)

print(as.data.frame(master_selected_table), row.names = FALSE)
print(as.data.frame(master_wild_table), row.names = FALSE)

# final_selected_filtered <- master_selected_table %>% 
#   filter(group != "22N9")
# final_wild_filtered <- master_wild_table %>% 
#   filter(group != "22DBW")

write.table(master_selected_table, "vis/master_selected_filtered.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(master_wild_table, "vis/master_wild_filtered.txt", sep="\t", row.names=FALSE, quote=FALSE)

######## recalculate 21dw
wild_all_gebv_status <- read.csv("gebvs_coxp/wild_all/Final_GEBVs_Summary.csv")
wildall_pheno        <- read.csv("data/wildlines/pheno/wild_all_pheno.csv")

pheno_21DW <- wildall_pheno %>%
  filter(Group == "21DW") %>%
  select(SampleID, status_01)

# Join GEBVs to phenotype status
gebv_21DW <- wild_all_gebv_status %>%
  inner_join(pheno_21DW, by = c("ID" = "SampleID"))

# Pearson correlations with status_01
cor_results <- gebv_21DW %>%
  summarise(
    RR_cor = cor(RR_GEBV_Mean, status_01, use = "complete.obs", method = "pearson"),
    EN_cor = cor(EN_GEBV_Mean, status_01, use = "complete.obs", method = "pearson"),
    RF_cor = cor(RF_GEBV_Mean, status_01, use = "complete.obs", method = "pearson"),
    GB_cor = cor(GB_GEBV_Mean, status_01, use = "complete.obs", method = "pearson")
  )
cor_results
