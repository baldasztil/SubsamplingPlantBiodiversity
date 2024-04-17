
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------

test <- "Lophopyxidaceae" 
# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(sf)
library(spmodel)
library(GWmodel)
library(feather)

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
        sp = unique(continent_data$sp)[1],
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
          sp = unique(continent_data$sp)[1],
          model_id = 0, 
          coeff = 1,
          coeffse = 0,
          pseudo.r.squared = 1
        )
    } else
      output_glance <- glance(output_glance_model) %>% 
        mutate(continent = x, 
               sp = unique(continent_data$sp)[1], 
               model_id = 2,
               coeff = output_glance_model$coefficients$fixed[2], 
               coeffse = sqrt(diag(vcov(output_glance_model)))[2]) 
    
    
    return(output_glance)
  }
  
  else {
    
    output_glance <- glance(output_glance_model) %>% 
      mutate(continent = x, 
             sp = unique(continent_data$sp)[1], 
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
             sp = unique(continent_data$sp)[1])
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

subsampling.plants <- function(fam, source) {
  number <- which(families == fam)

  
  print(paste0("This is the ", fam, " family"," [", number, "/", length(families), "]"))
        
        
  species_sample <- plants_full %>% 
    filter(family == fam)

  
  sp_origin <- dist_native %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>% 
    setDT() %>% 
    cube(., j= length(unique(plant_name_id)), by = "continent") %>% 
    mutate(continent = ifelse(is.na(continent), "GLOBAL", continent)) %>% 
    dplyr::select(continent, sp_origin = V1)
  
  dist <- source %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>% 
    left_join(continent_names, by = "LEVEL3_COD") %>% 
    rename(continent = LEVEL1_NAM)
  
  # richness patterns across brus
  rich_rel <- dist %>% 
    group_by(continent) %>% 
    mutate(sp = n_distinct(plant_name_id)) %>% 
    ungroup() %>%  
    group_by(LEVEL3_COD) %>% 
    summarise(richness_sample = n_distinct(plant_name_id), 
              sp = unique(sp)) %>% 
    right_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    right_join(rich_overall_bru, by = "LEVEL3_COD") %>%
    replace(is.na(.), 0) %>% 
    mutate(
           richness_sample = richness_sample,
           richness = richness) %>% 
    st_sf() 
  
  rich_rel_shp <- rich_rel %>% 
    as("Spatial")
  
  model <- try(splm(richness ~ richness_sample, data=rich_rel,
                    spcov_type =  c("exponential")))
  
  if (class(model)[1] =="try-error")  {
    global_er <- test %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        sp = length(unique(dist$plant_name_id)),
        model_id = -1, 
        family = fam) 
    return(global_er)
  } else 
    
    global <- as.data.table(glance(model)) %>% 
    mutate(continent = "GLOBAL", 
           sp = length(unique(dist$plant_name_id)), 
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
    left_join(sp_origin, by = "continent") %>% 
    
    mutate_if(is.numeric,
              round,
              digits = 4) %>% 
    dplyr::select(-n, -p, -npar, -value, -AIC) %>% 
    mutate(family = fam,
           sp_miss = sp_origin - sp) 
    
  
  #cumulative pattern 
  
  return(model_comp_out)
}

# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 

gbif_dist <- fread("data/dist_gbif.txt") 
gbif_dist_inat <- fread("data/ranges_gbif_inat_only.txt") 
gbif_dist_herb <- fread("data/ranges_gbif_herb_only.txt") 


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


# manipulate data --------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  left_join(plants_full, by = c("plant_name_id")) %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3)


plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)

nrow(plantlist_names)

rich_overall_bru <- richness_patterns

rich_overall_bru_mid <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  mutate(richness2 = richness) %>% 
  left_join(midpoints_red, by = "LEVEL3_COD") %>% 
  st_as_sf()

rich_overall_bru_shp <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>%
  mutate(richness2 = richness) %>% 
  st_as_sf() 

rich_overall_bru_mid_shp <- rich_overall_bru_mid %>% 
  as("Spatial") 

bw_global <-bw.gwr(richness ~ richness2, rich_overall_bru_mid_shp, 
                   approach = "AIC",
                   adaptive = T)


global <- glance(splm(richness ~ 1, data = rich_overall_bru_shp, spcov_type = "exponential"))

setDT(plantlist_names)
setDT(plantlist_dist)
setDT(gbif_dist)
setDT(rich_overall_bru)
setDT(midpoints_red)

families <- unique(plants_full$family)
names(families) <- "fam"

#families_inat <- plants_full %>% 
 # filter(plant_name_id %in% gbif_dist_inat$plant_name_id) %>% 
  #group_by(family) %>% 
  #summarise() 

#families %>% 
 # filter(!fam %in% families_inat$family)
  

test <- subsampling.plants("Poaceae", source = gbif_dist)
test
# compare models ---------------------------------------------------------------




result_list <- lapply(families, 
                      subsampling.plants, source = plantlist_dist)

results <- rbindlist(result_list) %>% 
  mutate(source = "wcvp")

res_file <- paste0("data/Samples_iteration_splm_rich_family_global.txt")
fwrite(results, file=res_file) 




result_list_gbif <- lapply(families, 
                           subsampling.plants, source = gbif_dist)

results_gbif <- rbindlist(result_list_gbif) %>% 
  mutate(source = "gbif")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif_global.txt")
fwrite(results_gbif, file=res_file)




result_list_gbif_inat <- lapply(families, 
                           subsampling.plants, source = gbif_dist_inat)

results_gbif_inat <- rbindlist(result_list_gbif_inat) %>% 
  mutate(source = "gbif_inat")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif_inat_global.txt")
fwrite(results_gbif_inat, file=res_file)




result_list_gbif_herb <- lapply(families, 
                                subsampling.plants, source = gbif_dist_herb)

results_gbif_herb <- rbindlist(result_list_gbif_herb) %>% 
  mutate(source = "gbif_herb")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif_herb_global.txt")
fwrite(results_gbif_herb, file=res_file)




results_full_raw <- rbind(results, results_gbif, results_gbif_inat, results_gbif_herb)
res_file <- paste0("data/Samples_iteration_splm_rich_family_full_global.txt")
fwrite(results_full_raw, file=res_file)

results_full <- results_full_raw %>% 
  filter(!model_id == 0 & !model_id == -1 )
  
threshold <- results_full %>% 
  group_by(source) %>% 
  filter(continent == "GLOBAL") %>% 
  slice_max(cor.sp_gwr_mean, n = 5)


ggplot(results_full, aes(x =cor.sp_gwr_mean, y = source, fill = source)) +
  geom_boxplot() +
  facet_wrap(~continent) +
  coord_flip()

ggplot(results_full, aes(x = log10(sp) , y = source, fill = source)) +
  geom_boxplot() +
  facet_wrap(~continent)

ggplot(results_full, aes(x = log10(sp_miss) , y = source, fill = source)) +
  geom_boxplot()+
  facet_wrap(~continent) +
  coord_flip()

library(purrr)

results_gbif_inat_join <- results_gbif_inat %>% 
  rename_with(~paste0(., "_inat"), -c(5, 20)) %>% 
  mutate(family = family_inat)

results_gbif_herb_join <- results_gbif_herb %>% 
  rename_with(~paste0(., "_herb"), -c(5, 20)) %>% 
  mutate(family = family_herb)

results_gbif_join <- results_gbif %>% 
  rename_with(~paste0(., "_gbif"), -c(5, 20)) %>% 
  mutate(family = family_gbif)



joined <- results %>% 
  left_join(results_gbif_inat_join, by = c("continent", "family")) %>% 
  left_join(results_gbif_herb_join, by = c("continent", "family")) %>% 
  left_join(results_gbif_join, by = c("continent", "family")) %>% 
  
  summarise(
    diff_all = abs(cor.sp_gwr_mean - cor.sp_gwr_mean_gbif), 
    
    diff_inat = abs(cor.sp_gwr_mean - cor.sp_gwr_mean_inat), 
    diff_herb = abs(cor.sp_gwr_mean - cor.sp_gwr_mean_herb),
    
    diff_sp_inat = sp - sp_inat, 
    diff_sp_herb = sp - sp_herb,
    diff_sp_gbif = sp - sp_gbif,
    
    family = family, 
    continent = continent, 
    sp_origin = results$sp_origin, 
    prop_diff_gbif = (1 - sp_gbif / sp_origin) *-1, 
  ) 

joined_glob <- joined %>% 
  filter(continent == "GLOBAL")

#http://127.0.0.1:44133/graphics/plot_zoom_png?width=1200&height=900

ggplot(joined, aes(y = diff_inat , x = continent, fill = continent)) +
  geom_boxplot() 

ggplot(joined, aes(y = diff_herb , x = continent, fill = continent)) +
  geom_boxplot()

ggplot(joined, aes(y = diff_all , x = continent, fill = continent)) +
  geom_boxplot()

ggplot(joined, aes(y = diff_sp_herb , x = continent, fill = continent)) +
  geom_boxplot()

ggplot(joined, aes(y = log10(diff_sp_gbif) , x = continent, fill = continent)) +
  geom_boxplot()


ggplot(joined, aes(y = log10(diff_sp_inat) )) +
  geom_boxplot()

ggplot(joined, aes(x = diff_sp_herb)) +
  geom_density()

results_full %>%                               # Summary by group using purrr
  split(.$source) %>%
  map(summary)

results_wide <- results %>% 
  left_join(results_gbif, by = c("family", "continent"))

sum_res <- results_wide %>%  
  group_by(family, continent) %>% 
  mutate(diff = abs(cor.sp_gwr_mean.x - cor.sp_gwr_mean.y)) 

stats <-  sum_res %>%
  ungroup() %>% 
  #group_by(continent) %>% 
  filter(sp_origin.x > 30) %>% 
  slice_max(diff, n = 25)

ggplot(sum_res, aes(x = diff, y = continent, fill = continent)) +
  geom_violin()




mean(sum_res$diff, na.rm = T)
sd(sum_res$diff, na.rm = T)
max(sum_res$diff, na.rm = T)




sum_res_glob <- sum_res %>% 
filter(continent == "GLOBAL") 
  
mean(sum_res_glob$diff, na.rm = T)
sd(sum_res_glob$diff, na.rm = T)
max(sum_res_glob$diff, na.rm = T)



results_full_glob <- results_full %>% 
  filter(continent == "GLOBAL")
results_full$pseudo.r.squared

boxplot(cor.sp_gwr_mean ~ source, data = results_full_glob)

boxplot(pseudo.r.squared ~ source, data = results_full_glob)
boxplot(mae ~ source, data = results_full_glob)


xx <- aov(cor.sp_gwr_mean ~ source, data = results_full_glob)
TukeyHSD(xx)

stats <-  results_full %>% 
  ungroup() %>% 
  #ilter(continent == "GLOBAL") %>% 
  group_by(continent, source) %>% 
  slice_max(cor.sp_gwr_mean,n = 3)


stats <-  results_full %>% 
filter(source == "wcvp") %>% 
  group_by(continent, source) %>% 
  filter(sp > 100) %>% 
  slice_max(cor.sp_gwr_mean,n = 3)

stats <-  results_full %>% 
  filter(source == "wcvp") %>% 
  group_by(continent, source) %>% 
  slice_max(pseudo.r.squared,n = 3)

stats <-  results_full %>% 
  filter(source == "wcvp") %>% 
  filter(cor.sp_gwr_mean > 0.9) %>% 
  filter(!continent == "ANTARCTICA")






