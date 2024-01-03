
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
std.error <- function(x) sd(x)/sqrt(length(x)) # function to caluclate st
subsampling.plants <- function(spec_n) {

  #extract random sample of species from the overall species list
  species_sample <- sample_n(plantlist_names, spec_n)
  
  # this file contains the distribution of each species with one row corresponding
  # to one botanical country
  dist <- dist_native %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  # here I calculate the number of distinct species per botanical country 
  rich_rel_shp <- dist %>% 
    group_by(area_code_l3) %>%  # botanical country code grouping
    summarise(richness_sample = n_distinct(plant_name_id)) %>%  # richness calculation
    rename(LEVEL_COD = area_code_l3) %>% 
    right_join(rich_overall_bru, by = "LEVEL_COD") %>% # join to dataset with original patterns
    left_join(midpoints_red, by = c("LEVEL_COD" = "LEVEL3_COD")) %>% # join to dataset containing midpoints of each country in order to take spatial relations into account 
    replace(is.na(.), 0) %>% # replace NAs with 0s to also take into account countries with no species sampled yet
    mutate(sp = nrow(species_sample),
           richness_sample = sqrt(richness_sample),
           richness = sqrt(richness)) %>%  # transform to achieve linearity 
    st_sf() %>% 
    as("Spatial")
  
  # at low species numbers the function breaks if it can not find the best 
  # possible bandwith or if there are multiple options, I use the try to prevent 
  # that from happening 
  bw <- try(bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
                   approach = "AIC",
                   adaptive = T))
  
  # if the bandwidth can not be calculated I use the global bandwidth as a proxy
  if (!is.numeric(bw)) {
    bw <- bw_global
  }
  
  # this calculates the GWR model 
  m.gwr <- try(gwr.basic(richness ~ richness_sample, data=rich_rel_shp,
                         adaptive = T,
                         bw = bw, 
                         cv = T, 
                         longlat = T))
  # this calculations correlation coefficicents 
  stats <- try(gwss(rich_rel_shp, vars = c("richness", "richness_sample"), adaptive = T, bw = bw))
  
  # again, at low sample numbers correlations might not be possible to calculate
  # the code below gives me an error dataframe with the same structure as the normal output 
  # at the moment this is set to give me 0 instead of NA's 
  if (class(m.gwr)=="try-error" | class(stats)=="try-error") {
    error_data <- as.data.frame(matrix(nrow = 10, ncol = 9))
    names(error_data) <- c("id","cor.gwr", "cor.sp","sd.cor.gwr", "se.cor.gwr","sample_coef","sample_se","sp", "shapiro")
    error_data$sp <- nrow(species_sample)
    error_data$id <- c("1","2","3","4","5","6","7","8","9","overall")
    df <- error_data %>% replace(is.na(.), 0)
    return(df)
  }
  # if it is possible to calculate the correlations the code below extracts the results 
  else {
    res_extract <- as.data.frame(m.gwr$SDF@data) %>% 
      dplyr::select(richness_sample_PRED = richness_sample,Local_R2, richness_sample_SE = richness_sample_SE, residual) %>% 
      mutate(LEVEL3_COD = rich_rel_shp$LEVEL_COD) %>% 
      left_join(tdwg_codes, by = "LEVEL3_COD") 
    res_extract$cor.sp <- stats$SDF$Spearman_rho_richness.richness_sample
  
  # this is for all continents separately   
  continents <- res_extract %>% 
      group_by(LEVEL1_COD) %>% 
      summarise(cor.gwr = round(mean(Local_R2), digits = 3),
                cor.sp = round(mean(cor.sp), digits = 3),
                sd.cor.gwr = round(sd(Local_R2), digits = 3),
                se.cor.gwr = round(std.error(Local_R2), digits = 3),
                sample_coef = round(mean(richness_sample_PRED), digits = 3), 
                sample_se = round(mean(richness_sample_SE), digits = 3)) %>% 
      rename(id = LEVEL1_COD)
    
  # here the overall patterns 
    overall <- res_extract %>% 
      summarise(cor.gwr =   round(mean(Local_R2), digits = 3),
                cor.sp = round(mean(cor.sp), digits = 3),
                sd.cor.gwr  = round(sd(Local_R2), digits = 3),
                se.cor.gwr  = round(std.error(Local_R2), digits = 3),
                sample_coef = round(mean(richness_sample_PRED), digits = 3), 
                sample_se = round(mean(richness_sample_SE), digits = 3), 
                id = "overall")
   
    # creating the output file by merging the two datasets  
    output_file <- rbind(continents, overall) %>% 
      mutate(sp = nrow(species_sample),
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

midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints") %>% 
  dplyr::select(LEVEL3_COD,geometry)

# manipulate data --------------------------------------------------------------
# remove names without accepted name
wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

# remove all hybrids and genus names 
wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) %>% 
  filter(species_hybrid == "") %>% 
  filter(genus_hybrid =="")

# filter distrbutions 
dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% # not introduced 
  filter(!location_doubtful == 1) %>%  # not doubtful occurrence 
  filter(!area == "") %>% # not unknown
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)

# summarise information
dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))

# create combined dataset
plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")

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


plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native %>% 
  filter(plant_name_id %in% plantlist_names$plant_name_id) %>% 
  dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)


# create spatial object to calculate global bandwidth
rich_overall_bru_bw <- richness_patterns %>% 
  filter(ID =="bru") %>%  
  mutate(richness2 = richness) %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  st_sf() %>% 
  as("Spatial") 

# global bandwidth 
bw_global <- bw.gwr(sqrt(richness) ~ sqrt(richness2), data=rich_overall_bru_bw,
       approach = "AIC",
       adaptive = T)

# global richness patterns as data frame 
rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>%  
  dplyr::select(-geometry)


# calculate correlations with a random sample at each step using lapply to do so sequentially
Sys.time()
samples <- lapply(seq(1,nrow(plantlist_names),1), subsampling.plants) 
Sys.time()

# export results as a .txt file 
xx <- rbindlist(samples)
write.table(xx, paste0("data/fullsamples_test/Samples_iteration_gwr_future_for",sample(1:10000000, 1, replace=TRUE),".txt")) 



