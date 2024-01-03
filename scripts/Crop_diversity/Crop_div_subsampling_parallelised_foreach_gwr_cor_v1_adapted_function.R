
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(parallel)
library(doParallel)
library(sf)
library(GWmodel)    # to undertake the GWR
library(foreach)
library(parallelly)
library(doFuture)





# Defining functions -----------------------------------------------------------
std.error <- function(x) sd(x)/sqrt(length(x))
subsampling.plants <- function(spec_n) {
  cumulative_namelist <- c()
  
  plantlist_names_left <- plantlist_names %>% 
    filter(!plant_name_id %in%  cumulative_namelist)
  if (nrow(plantlist_names_left) > spec_n) {
    species_sample <- sample_n(plantlist_names_left, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names_left
  }
  cumulative_namelist <- c(cumulative_namelist, species_sample$plant_name_id)
  species <-plantlist_names %>% 
    filter(plant_name_id %in%  cumulative_namelist)
  
  dist <- dist_native %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # richness patterns across brus
  rich_rel_shp <- dist %>% 
    group_by(area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3) %>% 
    right_join(rich_overall_bru, by = "LEVEL_COD") %>%
    left_join(midpoints_red, by = c("LEVEL_COD" = "LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist),
           richness_sample = sqrt(richness_sample),
           richness = sqrt(richness)) %>% 
    st_sf() %>% 
    as("Spatial")
  
  
  bw <- try(bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
                   approach = "AIC",
                   adaptive = T))
  
  if (!is.numeric(bw)) {
    bw <- bw_global
  }
  
  m.gwr <- try(gwr.basic(richness ~ richness_sample, data=rich_rel_shp,
                         adaptive = T,
                         bw = bw, 
                         cv = T, 
                         longlat = T))
  
  stats <- try(gwss(rich_rel_shp, vars = c("richness", "richness_sample"), adaptive = T, bw = bw))
  
  if (class(m.gwr)=="try-error" | class(stats)=="try-error") {
    error_data <- as.data.frame(matrix(nrow = 10, ncol = 9))
    names(error_data) <- c("id","cor.gwr", "cor.sp","sd.cor.gwr", "se.cor.gwr","sample_coef","sample_se","sp", "shapiro")
    error_data$sp <- length(cumulative_namelist)
    error_data$id <- c("1","2","3","4","5","6","7","8","9","overall")
    df <- error_data %>% replace(is.na(.), 0)
    return(df)
  }
  
  else {
    res_extract <- as.data.frame(m.gwr$SDF@data) %>% 
      dplyr::select(richness_sample_PRED = richness_sample,Local_R2, richness_sample_SE = richness_sample_SE, residual) %>% 
      mutate(LEVEL3_COD = rich_rel_shp$LEVEL_COD) %>% 
      left_join(tdwg_codes, by = "LEVEL3_COD") 
    res_extract$cor.sp <- stats$SDF$Spearman_rho_richness.richness_sample
    
    continents <- res_extract %>% 
      group_by(LEVEL1_COD) %>% 
      summarise(cor.gwr = round(mean(Local_R2), digits = 3),
                cor.sp = round(mean(cor.sp), digits = 3),
                sd.cor.gwr = round(sd(Local_R2), digits = 3),
                se.cor.gwr = round(std.error(Local_R2), digits = 3),
                sample_coef = round(mean(richness_sample_PRED), digits = 3), 
                sample_se = round(mean(richness_sample_SE), digits = 3)) %>% 
      rename(id = LEVEL1_COD)
    
    overall <- res_extract %>% 
      summarise(cor.gwr =   round(mean(Local_R2), digits = 3),
                cor.sp = round(mean(cor.sp), digits = 3),
                sd.cor.gwr  = round(sd(Local_R2), digits = 3),
                se.cor.gwr  = round(std.error(Local_R2), digits = 3),
                sample_coef = round(mean(richness_sample_PRED), digits = 3), 
                sample_se = round(mean(richness_sample_SE), digits = 3), 
                id = "overall")
    
    output_file <- rbind(continents, overall) %>% 
      mutate(sp = length(cumulative_namelist),
             shapiro = shapiro.test(res_extract$residual)$p.value)
    
    
    return(output_file)
  }
}


# Working directory ------------------------------------------------------------
message("Number of CPU cores in R: ", parallelly::availableCores())

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 


tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints") %>% 
  dplyr::select(LEVEL3_COD,geometry)

# manipulate data --------------------------------------------------------------

wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

wcvp_sp <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(!taxon_rank %in% c("Genus")) 

wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) %>% 
  filter(species_hybrid == "") %>% 
  filter(genus_hybrid =="")


dist_native <- dist_raw %>% 
  filter(introduced == 0) %>%
  filter(!location_doubtful == 1) %>% 
  filter(!area == "") %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)

dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))

plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")

plants_full_extinct_norange <- wcvp_accepted %>% 
  filter(!plant_name_id %in% plants_full$plant_name_id)

# baseline richness ------------------------------------------------------------
# calculating overall richness patterns 

richness_patterns_con <- dist_native %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- dist_native %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- dist_native %>% 
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


# Analysis ---------------------------------------------------------------------
plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)



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



######
######

# Define the number of cores you want to use
#num_cores <- 10  # Adjust this based on your system

# Create a cluster
#cl <- makeCluster(num_cores)


# Register the cluster with foreach
#registerDoParallel(cl)

plan(multisession)


Sys.time()

samples <- foreach(i = 1:4, .options.future = list(seed = TRUE, scheduling = 1.0)) %dofuture% { 
   lapply(seq(1,1000,1), subsampling.plants)
}

Sys.time()

xx <- rbindlist(samples)
write.table(xx, paste0("data/fullsamples_test/Samples_iteration_gwr_future_for",sample(1:10000000, 1, replace=TRUE),".txt")) 


stopCluster(cl)


