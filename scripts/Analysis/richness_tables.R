
# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(reactablefmtr)

# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 


plants_full <- fread("data/wcvp_accepted_merged.txt")
richness_patterns_raw <- fread("data/richness_patterns.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

continent_names_vec <- unique(continent_names$LEVEL1_NAM)

richness_patterns <- richness_patterns_raw %>% 
  left_join(continent_names, by = "LEVEL3_COD")



#Summary statistics analysis----------------------------------------------------


richness_patterns_con <- dist_native %>%
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  mutate(LEVEL_COD = as.character(continent_code_l1)) %>%  
  mutate(ID = "con")


index_endemic <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n == 1) # 56%

richness_patterns_con_en <- dist_native %>%
  filter(plant_name_id %in% index_endemic$plant_name_id) %>% 
  group_by(continent_code_l1) %>% 
  summarise(rich_en = n_distinct(plant_name_id))%>% 
  mutate(LEVEL_COD = as.character(continent_code_l1)) %>%  
  mutate(ID = "con")

index_widespread <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) # 43.9%

richness_patterns_con_nonen <- dist_native %>%
  filter(plant_name_id %in% index_widespread$plant_name_id) %>% 
  group_by(continent_code_l1) %>% 
  summarise(rich_wide = n_distinct(plant_name_id))%>% 
  mutate(LEVEL_COD = as.character(continent_code_l1)) %>% 
  mutate(ID = "con")


richvsamplesize <- threshold %>%
  filter(!continent == "GLOBAL") %>% 
  left_join(richness_patterns, by = join_by("id" == "LEVEL_COD")) %>%
  left_join(richness_patterns_con_en, by = join_by("id" == "LEVEL_COD")) %>%
  left_join(richness_patterns_con_nonen, by = join_by("id" == "LEVEL_COD")) %>% 
  select(continent, sp, richness, rich_en, rich_wide) %>% 
  mutate(prop_en = rich_en / (rich_en + rich_wide), 
         prop_wide = rich_wide / (rich_en + rich_wide)) %>% 
  ungroup() %>%
  mutate(sp = sp,
         ratio =   log10(sp / richness),
         rank_rich = rank(richness),
         rank_en = rank(rich_en),
         rank_wide = rank(rich_wide),
         rank_threshold =  rank(sp))



GLOBAL_patterns <- dist_native %>% 
  summarise(continent = "GLOBAL", 
            sp = as.numeric(threshold[8,2]),
            richness = n_distinct(plant_name_id),
            rich_en = n_distinct(index_endemic$plant_name_id),
            rich_wide = n_distinct(index_widespread$plant_name_id),
            prop_en = rich_en / (rich_en + rich_wide), 
            prop_wide = rich_wide / (rich_en + rich_wide),
            ratio =  log10(sp / richness)) %>%
  bind_rows(richvsamplesize) %>% 
  replace(is.na(.), 0) 

rich_summary <- GLOBAL_patterns %>% 
  rename(threshold = sp, ratio_rich_tresh = ratio) 

rich_table <- reactable(
  rich_summary,
  defaultSorted = "threshold",
  defaultSortOrder = "desc",
  defaultColDef = colDef(
    cell = data_bars(rich_summary, 
                     #number_fmt = scales::number_format(accuracy = 0.01),
                     fill_color = rev(viridis::viridis(100)), 
                     #background = "grey",
                     text_position = "above", 
                     bar_height = 6, 
                     #fill_gradient = TRUE
                     box_shadow = T), align = "center", width = 150
    
  ),
  columns = list(
    continent = colDef(show = T),
    threshold = colDef( show = T, format = colFormat(digits = 0)),
    richness = colDef( show = T, format = colFormat(digits = 0)), 
    rich_en = colDef( show = T, format = colFormat(digits = 0)),
    rich_wide = colDef( show = T, format = colFormat(digits = 0)),
    prop_en = colDef( show = T, format = colFormat(digits = 2)),
    prop_wide = colDef( show = T, format = colFormat(digits = 2)),
    ratio_rich_thresh = colDef( show = T, format = colFormat(digits = 1)),
    rank_rich = colDef( show = T, format = colFormat(digits = 0)),
    rank_en = colDef( show = T, format = colFormat(digits = 0)),
    rank_wide = colDef( show = T, format = colFormat(digits = 0)),
    rank_threshold = colDef( show = T, format = colFormat(digits = 0))
  )
)

rich_table
#save_reactable(rich_table, "rich__summary_table.png", vwidth = 1900)
