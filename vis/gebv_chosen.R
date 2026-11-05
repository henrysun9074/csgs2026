library(tidyverse)
library(readxl)

available <- read_excel("data/spawn/LiveOystersWithTags.xlsx")
chosen <- read.csv("data/spawn/OystersUsedJun17.csv")
chosen <- chosen %>% filter(obtained == "y" & sex != "ng")
chosen_wild <- chosen %>% filter(group == "22DBWU" | group == "22DBWD")
chosen_sel <- chosen %>% filter(group == "22N9U" | group == "22N9D")

wild_avail <- available$"22DBW"
wild_avail <- wild_avail[!is.na(wild_avail)]
sel_avail <- available$"22N9"
sel_avail <- sel_avail[!is.na(sel_avail)]

## phenotype
wild_pheno <- read.csv("data/wildlines/pheno/wild_22dbw_pheno.csv")
wild_pheno_alive <- wild_pheno %>% filter(status_01 == 1)
sel_pheno <- read.csv("data/selectedlines/pheno/sel_n9_pheno.csv")
sel_pheno_alive <- sel_pheno %>% filter(status_01 == 1)

sel_gebv <- read.csv("gebvs_coxp/sel_all/Final_GEBVs_Summary.csv") #pooled cox
wild_gebv <- read.csv("gebvs_v2/wild_all/Final_GEBVs_Summary.csv") #pooled status

# join wild_pheno and sel_gebv on "SampleID" and "ID"
# make a new df for sel and wild with only rows with status_01 == 1
# compute a "rank_up" column for sel and wild_gebv$GB_GEBV_Mean from highest to lowest
# compute a "rank_down" for GB_GEBV_Mean from lowest to highest
# join the rank_up, rank_down, GB_GEBV_Mean, GB_GEBV_SD, phenotype_gwas, status_01 columns
wild_combined <- wild_gebv %>% 
  inner_join(wild_pheno, by = c("ID" = "SampleID")) %>% 
  filter(status_01 == 1) %>% 
  mutate(
    rank_up   = min_rank(desc(GB_GEBV_Mean)), # Highest to lowest
    rank_down = min_rank(GB_GEBV_Mean)        # Lowest to highest
  )

sel_combined <- sel_gebv %>% 
  inner_join(sel_pheno, by = c("ID" = "SampleID")) %>% 
  filter(status_01 == 1) %>% 
  mutate(
    rank_up   = min_rank(desc(GB_GEBV_Mean)), # Highest to lowest
    rank_down = min_rank(GB_GEBV_Mean)        # Lowest to highest
  )

# Clean IDs in the reference datasets
wild_combined <- wild_combined %>% mutate(ID_clean = sub("_.*", "", ID))
sel_combined  <- sel_combined  %>% mutate(ID_clean = sub("_.*", "", ID))
chosen_wild <- chosen_wild %>% mutate(ID_clean = sub("_.*", "", ID))
chosen_sel  <- chosen_sel  %>% mutate(ID_clean = sub("_.*", "", ID))

target_cols <- c("ID_clean", "rank_up", "rank_down", "GB_GEBV_Mean", "GB_GEBV_SD", "phenotype_gwas", "status_01")


# figure out what is getting dropped here... should be 39
final_chosen_wild <- chosen_wild %>% 
  left_join(
    wild_combined %>% select(any_of(target_cols)), 
    by = "ID_clean"
  )
final_chosen_wild <- final_chosen_wild %>% select(-phenotype_gwas)

final_chosen_sel <- chosen_sel %>% 
  left_join(
    sel_combined %>% select(any_of(target_cols)), 
    by = "ID_clean"
  )

# inner join phenotype_gwas column from sel_pheno and wild_pheno to final chosen wild/sel  
wild_gwas_lookup <- wild_pheno %>%
  mutate(ID_clean = sub("_.*", "", SampleID)) %>%
  select(ID_clean, phenotype_gwas)

sel_gwas_lookup <- sel_pheno %>%
  mutate(ID_clean = sub("_.*", "", SampleID)) %>%
  select(ID_clean, phenotype_gwas)

final_chosen_wild <- final_chosen_wild %>%
  inner_join(wild_gwas_lookup, by = "ID_clean")

final_chosen_sel <- final_chosen_sel %>%
  inner_join(sel_gwas_lookup, by = "ID_clean")

final_chosen <- rbind(final_chosen_wild, final_chosen_sel)

growth <- read_excel("data/spawn/growth_data.xlsx")
growth$Tag <- toupper(growth$Tag)
colnames(growth)[4] <- "ID_clean"
growth <- growth %>% 
  mutate(ID_clean = str_pad(as.character(ID_clean), width = 3, pad = "0"))

final_chosen <- final_chosen %>%
  left_join(
    growth %>% select(ID_clean, Height, Weight), 
    by = "ID_clean"
  )

write.csv(final_chosen, "data/spawn/chosen_oysters.csv")

############# analysis TBD
# compare size/weight between groups
# look at difference in size/weight from last workover between up and downselect
# compare growth rates from previous workovers
# look at breeding values if we had done phenotypic selection (top 30 highest/lowest weight + growth)

# look where model predictions would have differed along three axes: 
# 1) solo vs pooled training
# 2) status vs cox
# 3) different model architectures - need to read in separate files

# get survival data - compare model accuracies from last workover to spawn: how many mortalities? 
# did those oysters that died have low GEBV?