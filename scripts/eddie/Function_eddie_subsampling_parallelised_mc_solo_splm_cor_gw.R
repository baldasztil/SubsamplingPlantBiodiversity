# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(spmodel)
library(vegan)
library(GWmodel)



# Defining functions -----------------------------------------------------------
splm.continents <-  function(x, rich_rel)  {
  
  continent_data <- rich_rel %>% 
    filter(LEVEL1_NAM == x) 
  
  output_glance_model <- try(splm(richness ~ richness_sample, data=continent_data,
                                  spcov_type =  c("exponential")))
  
  if (class(output_glance_model)[1] =="try-error" & continent_data$sp[1] < 50000)  {
    output_glance <- global %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        continent = unique(continent_data$LEVEL1_NAM),
        sp = unique(continent_data$sp),
        model_id = 0, 
        coeff = 0,
        coeffse = 0
      )
    return(output_glance)
  }
  
  if (class(output_glance_model)[1] =="try-error" & continent_data$sp[1] > 50000)  {
    
    output_glance_model <- try(splm(richness ~ richness_sample, data=continent_data,
                                    spcov_type =  c("gravity")))
    
    if (class(output_glance_model)[1] =="try-error") {
      
      output_glance <- global %>% 
        replace(!is.numeric(.) ==T, 0) %>% 
        mutate( 
          continent = unique(continent_data$LEVEL1_NAM),
          sp = unique(continent_data$sp),
          model_id = 0, 
          coeff = 1,
          coeffse = 0,
          pseudo.r.squared = 1
        )
    } else
      output_glance <- glance(output_glance_model) %>% 
        mutate(continent = x, 
               sp = unique(continent_data$sp), 
               model_id = 2,
               coeff = output_glance_model$coefficients$fixed[2], 
               coeffse = sqrt(diag(vcov(output_glance_model)))[2]) 
    
    
    return(output_glance)
  }
  
  else {
    
    output_glance <- glance(output_glance_model) %>% 
      mutate(continent = x, 
             sp = unique(continent_data$sp), 
             model_id = 1,
             coeff = output_glance_model$coefficients$fixed[2], 
             coeffse = sqrt(diag(vcov(output_glance_model)))[2])
    
  }
  return(output_glance)
}
cor.gw <- function(rich_rel_shp) {
  bw <- try(bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
                   adaptive = T, 
                   longlat = T))
  
  if (!is.numeric(bw)) {
    bw <- bw_global
  }
  
  stats <- try(as.data.frame(gwss(rich_rel_shp, rich_rel_shp, 
                                  vars = c("richness", "richness_sample"),
                                  adaptive = T, bw = bw, longlat = T)$SDF))
  
  if (class(stats)[1]=="try-error") {
    error_data <- output_model %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate(id = output_model$id, 
             sp = nrow(species_sample))
    return(error_data)
  }
  
  else {
    
    stats_total <-  cbind(stats, rich_rel_shp) %>% 
      dplyr::select(cor.sp_gwr =Spearman_rho_richness.richness_sample, LEVEL3_COD, 
                    continent = LEVEL1_NAM) %>% 
      replace(is.na(.), 0) %>% 
      setDT()
    
    funs <- lst(mean, sd)
    
    global <- stats_total %>%  
      summarise(across(
        .cols = where(is.numeric), 
        .fns = list(mean = mean, sd = sd, med = median),  
        .names = "{col}_{fn}"
      )) %>% 
      mutate(continent = "GLOBAL") 
    
    continents <- stats_total %>% 
      group_by(continent) %>% 
      summarise(across(
        .cols = where(is.numeric), 
        .fns = list(mean = mean, sd = sd, med =median), 
        .names = "{col}_{fn}"
      )) 
    
    output_file <- rbind(continents, global) 
    return(output_file)
  }
}

subsampling.plants <- function(spec_n) {
  
  if (nrow(plantlist_names) > spec_n) {
    species_sample <- sample_n(plantlist_names, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names
  }
  
  dist <- plantlist_dist %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  table_data_sample <- table(dist$LEVEL3_COD, dist$growth_form)
  shannon_div_sample <- as.data.frame(diversity(table_data_sample, index = "shannon")) 
  shannon_div_sample$C <- rownames(shannon_div_sample)
  names(shannon_div_sample) <- c("richness_sample", "LEVEL3_COD")
  
  #plot(rich_rel$richness_sample ~ rich_rel$richness)
  
  rich_rel <- shannon_div_sample  %>% 
    right_join(rich_overall_bru, by = c("LEVEL3_COD")) %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = nrow(species_sample),
           richness_sample = richness_sample,
           richness = richness) %>% 
    st_sf()  
  #filter(!LEVEL1_NAM == "ANTARCTICA")
  
  rich_rel_shp <- rich_rel %>% 
    as("Spatial")
  
  model <- try(splm(richness ~ richness_sample, data=rich_rel,
                    spcov_type =  c("exponential")))
  
  if (class(model)[1] =="try-error")  {
    global_er <- test %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        sp = nrow(species_sample),
        model_id = -1) 
    return(global_er)
  } else 
    
    global <- as.data.table(glance(model)) %>% 
    mutate(continent = "GLOBAL", 
           sp = length(species_sample$plant_name_id), 
           model_id = 1, 
           coeff = model$coefficients$fixed[2],
           coeffse = sqrt(diag(vcov(model)))[2])
  
  rich_rel_preds <- as.data.table(model$obdata)[, c("LEVEL1_NAM", "LEVEL3_COD", "richness", "richness_sample")] %>% 
    mutate(preds =  fitted.values(model),  
           residuals =   residuals(model)) 
  
  model_error <- suppressWarnings(cube(rich_rel_preds, j = list(sum(abs(richness - preds))/length(richness),
                                                                mean(preds), sum(residuals^2), mean(residuals), 
                                                                sd(residuals), min(residuals), max(residuals)), 
                                       by = c("LEVEL1_NAM")))
  
  colnames(model_error) <- c("continent", "mae", "mpred","sumres2", "mres", "sdres", "minres", "maxres")
  model_error[is.na(model_error)] <- "GLOBAL"
  
  model_comp <- rbindlist(lapply(continent_names_vec, splm.continents, 
                                 rich_rel = rich_rel )) %>%
    rbind(global) 
  
  correlation <- cor.gw(rich_rel_shp = rich_rel_shp)
  
  model_comp_out <- model_comp %>% 
    left_join(correlation, by = "continent") %>% 
    left_join(model_error, by = "continent") %>% 
    mutate_if(is.numeric,
              round,
              digits = 4) %>% 
    dplyr::select(-n, -p, -npar, -value, -AIC)
  
  #cumulative pattern
  if (nrow(species_sample) %in% seq(1,nrow(plantlist_names),50000)){
    tf1 <- paste0("/exports/eddie/scratch/s2117440/data/functional/checkpoints/checkpoint_fun_gbif",
                  nrow(species_sample),".txt")
    fwrite(model_comp_out, tf1) 
  }
  
  return(model_comp_out)
}

# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full <- fread("data/wcvp_accepted_merged.txt")
richness_patterns_raw <- fread("data/richness_patterns.txt")

growth <- read.csv("data/growths_rf.csv") %>% 
  dplyr::select(growth_form = apd_plant_growth_form, plant_name_id)

plants_full <- plants_full %>% 
  left_join(growth, by = "plant_name_id")

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()



dist_native_calc <- dist_native %>% 
  left_join(plants_full, by = "plant_name_id")

table_data <- table(dist_native_calc$area_code_l3, dist_native_calc$growth_form)


shannon_div <- as.data.frame(diversity(table_data, index = "shannon")) 
shannon_div$C <- rownames(shannon_div)
names(shannon_div) <- c("richness", "LEVEL3_COD")

richness_patterns <- shannon_div %>% 
  left_join(continent_names, by = "LEVEL3_COD") #%>% 
  #filter(!LEVEL1_NAM == "ANTARCTICA")

continent_names_vec <- unique(richness_patterns$LEVEL1_NAM)


# manipulate data -------------------------------------------------------- ------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)

plantlist_dist <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3, growth_form)


rich_overall_bru <- richness_patterns

rich_overall_bru_mid <- richness_patterns %>% 
  left_join(midpoints_red, by = "LEVEL3_COD")  %>% 
  mutate(richness2 = richness) %>% 
  st_as_sf()

rich_overall_bru_shp <- richness_patterns %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()

rich_overall_bru_mid_shp <- rich_overall_bru_mid %>% 
  as("Spatial") 

bw_global <-bw.gwr(richness ~ richness2, rich_overall_bru_mid_shp, 
                   approach = "AIC",
                   adaptive = T)

global <- glance(splm(richness ~ 1, data = rich_overall_bru_shp, spcov_type = "exponential"))
test <- subsampling.plants(30000)

# compare models ---------------------------------------------------------------
setDT(plantlist_names)
setDT(dist_native)
setDT(rich_overall_bru)
setDT(midpoints_red)


start <- Sys.time()
subsampling.plants(349113)
Sys.time() - start

start <- Sys.time()
subsampling.plants(80000)
Sys.time() - start

start <- Sys.time()
subsampling.plants(1)
Sys.time() - start

result_list <- lapply(seq(1,nrow(plantlist_names),1), 
                      subsampling.plants)

results <- rbindlist(result_list)

res_file <- paste0("/exports/eddie/scratch/s2117440/data/functional/Samples_iteration_splm_fun",
                   sample(1:10000000, 1, replace=TRUE),".txt.gz")
fwrite(results, file=res_file,  compress = "gzip") 

