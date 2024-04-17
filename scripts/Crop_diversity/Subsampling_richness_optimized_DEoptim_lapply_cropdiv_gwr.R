library(DEoptim)
library(tidyverse)
library(data.table)
library(sf)
library(foreach)
library(doSNOW)
library(GWmodel)


output <- function(species=NULL,best_corr=NULL){
  me <- data.frame(
    Species = species,
    best_corr = best_corr
  )
  
  return(me)
}
objective <- function(subset_indices) {
  species <- slice_sample(n = round(subset_indices),plantlist_names)
  dist <- plantlist_dist %>%
    filter(plant_name_id %in% species$plant_name_id)
  
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
 
   rich_rel_shp <- rich_overall_bru %>% 
    right_join(sample_rich_bru, by = "LEVEL_COD") %>%
    left_join(midpoints_red, by = c("LEVEL_COD" = "LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(species)) %>% 
    st_sf() %>% 
    as("Spatial")
   
   
   bw <- try(bw.gwr(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel_shp,
                    approach = "AIC",
                    adaptive = T)) 
   
   if (!is.numeric(bw)) {
     bw <- bw_global
   }
   
   m.gwr <- try(gwr.basic(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel_shp,
                          adaptive = T,
                          bw = bw, 
                          cv = T, 
                          longlat = T)) 
  
   
   if (class(m.gwr)=="try-error") {
     error_data <- 10000
     return(error_data)
   }
   else {
     # m.gwr <- gwr(richness ~ richness_sample, data=rich_rel_shp,
     #             adapt = bw_spgwr, 
     #                  gweight = gwr.bisquare) 
     
    res_extract <- as.data.frame(m.gwr$SDF@data)
    cor.gwr <- mean(res_extract$Local_R2)
   }
    
    if (cor.gwr > max(maxima)) {
      assign("maxima", cor.gwr, envir = .GlobalEnv)
      #maxima <- list(maxima, as.list(cor.gwr))
      ids <- unique(species$plant_name_id) 
      write.table(c(cor.gwr,ids), paste0(length(ids),"_target_population.txt"))
      } 

  return(-cor.gwr)  # Return negative correlation (since DEoptim minimizes)
  #return(ifelse(correlation >= 0.95, correlation, 0))
}
deoptim_lapply <- function(sp_num){
  assign("maxima", 0, envir = .GlobalEnv)
  de_result <- DEoptim(objective, lower = sp_num, upper = sp_num,
                       DEoptim.control(VTR = -Inf, strategy = 2,
                                       bs = FALSE, itermax = 10, CR = 0.5, F = 0.8, 
                                       storepopfreq = 1, trace = T, c = 0.5, steptol = 100))
  
  return(output(de_result$optim$bestmem, de_result$optim$bestval*-1))
}

# Working directory ------------------------------------------------------------
# Import data ------------------------------------------------------------------
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

richness_patterns <- fread("data/richness_patterns.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")

midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints") %>% 
  dplyr::select(LEVEL3_COD,geometry)

# creating objects for the function

plantlist_names <- plants_full %>% 
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)


richness_patterns_con <- plantlist_dist %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- plantlist_dist %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- plantlist_dist %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

richness_patterns1 <- richness_patterns_con %>% 
  rename( LEVEL_COD = LEVEL1_COD)
richness_patterns2 <- richness_patterns_reg %>% 
  rename( LEVEL_COD = LEVEL2_COD)
richness_patterns3 <- richness_patterns_bru %>% 
  rename( LEVEL_COD = LEVEL3_COD)

richness_patterns <- rbind(richness_patterns1,richness_patterns2, 
                           richness_patterns3)

rich_overall_bru_bw <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  mutate(richness2 = richness) %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  st_sf() %>% 
  as("Spatial") 


bw_global <- bw.gwr(sqrt(richness) ~ sqrt(richness2), data=rich_overall_bru_bw,
                    approach = "AIC",
                    adaptive = T)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>%  
  dplyr::select(-geometry)



# Create a cluster
max_sp <- 5123

solutions_listed <- mclapply(seq(1, max_sp, 1), mc.cores = 30, deoptim_lapply)

deoptim_lapply(100)

solutions <- do.call(rbind, solutions_listed)
write.table(solutions, "data/optimised/DEoptim_results.txt")
