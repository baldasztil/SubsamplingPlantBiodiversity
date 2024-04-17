
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
  
  if (class(output_glance_model)[1] =="try-error" & sum(continent_data$richness) < 5000)  {
    output_glance <- global %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        continent = unique(continent_data$LEVEL1_NAM),
       
        model_id = 0, 
        coeff = 0,
        coeffse = 0
      )
    return(output_glance)
  }
  
  if (class(output_glance_model)[1] =="try-error" & sum(continent_data$richness) > 5000)  {
    
    output_glance_model <- try(splm(richness ~ richness_sample, data=continent_data,
                                    spcov_type =  c("gravity")))
    
    if (class(output_glance_model)[1] =="try-error") {
      
      output_glance <- global %>% 
        replace(!is.numeric(.) ==T, 0) %>% 
        mutate( 
          continent = unique(continent_data$LEVEL1_NAM),
         
          model_id = 0, 
          coeff = 1,
          coeffse = 0,
          pseudo.r.squared = 1
        )
    } else
      output_glance <- glance(output_glance_model) %>% 
        mutate(continent = x, 
               
               model_id = 2,
               coeff = output_glance_model$coefficients$fixed[2], 
               coeffse = sqrt(diag(vcov(output_glance_model)))[2]) 
    
    
    return(output_glance)
  }
  
  else {
    
    output_glance <- glance(output_glance_model) %>% 
      mutate(continent = x, 
             
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


subsampling.plants.fam <- function(fam, source) {
  number <- which(families == fam)
  
  
  print(paste0("This is the ", fam, " family"," [", number, "/", length(families), "]"))
  
  
  species_sample <- plants_full %>% 
    filter(family == fam)
  
  dist <- source %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>%  
    mutate(family = fam)
  
  dist_sp <- source %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>% 
    left_join(continent_names, by = "LEVEL3_COD") %>% 
    rename(continent = LEVEL1_NAM) %>% 
    group_by(continent) %>% 
    summarise(sp_dataset = n_distinct(plant_name_id))
  
  dist_nat <- dist_native %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>% 
    group_by(continent) %>% 
    summarise(sp_origin = n_distinct(plant_name_id)) %>% 
    left_join(dist_sp, by = "continent") %>% 
    replace(is.na(.), 0) 
  

  
  rich_overall_fam  <-rich_overall_bru %>% 
    filter(family == fam) %>% 
    dplyr::select(c(-ID))
     
  # richness patterns across brus
  rich_rel <- dist %>% 
    group_by(LEVEL3_COD, family) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    right_join(rich_overall_fam, by = c("LEVEL3_COD", "family")) %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    left_join(continent_names, by = c("LEVEL3_COD")) %>% 
    mutate(family = family) %>% 
    replace(is.na(.), 0) %>% 
    mutate(
      richness_sample = richness_sample,
      richness = richness, 
      dif = richness - richness_sample, 
      diff_prop = (1 - richness_sample / richness) * -1, 
      countries = nrow(.)) %>% 
    st_sf() 
  
  rich_rel_shp <- rich_rel %>% 
    as("Spatial")
  
  model <- try(splm(richness ~ richness_sample, data=rich_rel,
                    spcov_type =  c("exponential")))
  
  if (class(model)[1] =="try-error")  {
    global_er <- test %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        sp_origin = length(unique(dist$plant_name_id)),
        model_id = -1, 
        family = fam) 
    return(global_er)
  } else 
    
    global <- as.data.table(glance(model)) %>% 
    mutate(continent = "GLOBAL", 
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
  
  model_comp <- rbindlist(lapply(dist_sp$continent, splm.continents, 
                                 rich_rel = rich_rel )) %>%
    rbind(global) 
  
  correlation <- cor.gw(rich_rel_shp = rich_rel_shp)
  
  model_comp_out <- model_comp %>% 
    left_join(correlation, by = "continent") %>% 
    left_join(model_error, by = "continent") %>% 
    mutate_if(is.numeric,
              round,
              digits = 4) %>% 
    dplyr::select(-n, -p, -npar, -value, -AIC) %>% 
    mutate(family = fam) %>% 
    left_join(dist_nat, by = c("continent"))
  model_comp_out$sp_dataset[model_comp_out$continent == "GLOBAL"]  <- length(unique(dist$plant_name_id))
  model_comp_out$sp_origin[model_comp_out$continent == "GLOBAL"]  <- length(unique(species_sample$plant_name_id))
  #cumulative pattern 
  
  return(model_comp_out)
}
subsampling.plants.fam.mapping <- function(fam, source) {
  number <- which(families == fam)
  
  
  print(paste0("This is the ", fam, " family"," [", number, "/", length(families), "]"))
  
  
  species_sample <- plants_full %>% 
    filter(family == fam)
  
  dist <- source %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>%  
    mutate(family = fam)
  
  rich_overall_fam  <-rich_overall_bru %>% 
    filter(family == fam) %>% 
    dplyr::select(c(-ID))
  
  # richness patterns across brus
  rich_rel <- dist %>% 
    group_by(LEVEL3_COD, family) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    right_join(rich_overall_fam, by = c("LEVEL3_COD", "family")) %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    left_join(continent_names, by = c("LEVEL3_COD")) %>% 
    mutate(family = fam) %>% 
    replace(is.na(.), 0) %>% 
    mutate(
      richness_sample = richness_sample,
      richness = richness, 
      dif = richness - richness_sample, 
      diff_prop = (1 - richness_sample / richness) * -1, 
      countries = nrow(.))  
  
  
  return(rich_rel)
}



# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 

gbif_dist <- fread("data/dist_gbif.txt") 
gbif_dist_inat <- fread("data/ranges_gbif_inat_only.txt") 
gbif_dist_herb <- fread("data/ranges_gbif_herb_only.txt") 

richness_patterns <- fread("data/richness_patterns_fam.txt")

plants_full <- fread("data/wcvp_accepted_merged.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

continent_names_vec <- unique(continent_names$LEVEL1_NAM)


plantlist_dist <- dist_native %>% 
  left_join(plants_full, by = c("plant_name_id")) %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3)

midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)

nrow(plantlist_names)

rich_overall_bru  <- richness_patterns
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
rich_overall_bru_mid_shp_test <- rich_overall_bru_mid %>% 
  filter(family == "Poaceae") %>% 
  as("Spatial") 



bw_global <-bw.gwr(richness ~ richness2, rich_overall_bru_mid_shp_test, 
                   approach = "AIC",
                   adaptive = T)


global <- glance(splm(richness ~ 1, data = rich_overall_bru_shp %>% 
                        filter(family == "Poaceae"), spcov_type = "exponential"))

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


test <- subsampling.plants.fam("Fabaceae", source = gbif_dist)
test2 <- subsampling.plants.fam.mapping("Fabaceae", source = gbif_dist)
# compare models ---------------------------------------------------------------


result_list_gbif <- lapply(families, 
                           subsampling.plants.fam, source = gbif_dist)

results_gbif <- rbindlist(result_list_gbif) %>% 
  mutate(source = "gbif")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif.txt")
fwrite(results_gbif, file=res_file)


result_list_gbif_inat <- lapply(families, 
                                subsampling.plants.fam, source = gbif_dist_inat)

results_gbif_inat <- rbindlist(result_list_gbif_inat) %>% 
  mutate(source = "gbif_inat")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif_inat.txt")
fwrite(results_gbif_inat, file=res_file)




result_list_gbif_herb <- lapply(families, 
                                subsampling.plants.fam, source = gbif_dist_herb)

results_gbif_herb <- rbindlist(result_list_gbif_herb) %>% 
  mutate(source = "gbif_herb")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif_herb.txt")
fwrite(results_gbif_herb, file=res_file)





results_full_raw <- rbind(results_gbif, results_gbif_inat, results_gbif_herb)
res_file <- paste0("data/results_full_raw_familieswithin.txt")

fwrite(results_full_raw, file=res_file)


boxplot(results_full_raw$cor.sp_gwr_mean ~ results_full_raw$source)
boxplot(results_full_raw$cor.sp_gwr_sd ~ results_full_raw$source)

results_full_rem<- results_full_raw %>% 
  filter(!model_id == -1 )

aa <- results_full_rem %>% 
  filter(!model_id == -1) %>% 
  filter(continent == "GLOBAL") %>% 
  mutate(diff = sp_dataset - sp_origin) %>% 
  filter(sp_origin > 100) %>% 
  group_by(source) %>% 
  slice_min(cor.sp_gwr_mean, n = 10)

plot(aa$cor.sp_gwr_mean ~ sqrt(aa$diff*-1)) 

library(ggdist)
ggplot(results_full_rem, aes(y = cor.sp_gwr_mean, x = continent)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha = 0.25, aes(col = continent)) +
  #stat_pointinterval() +
  facet_wrap(~source) +
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  coord_flip()

ggplot(results_full_rem, aes(y = log10(sp_origin - sp_dataset), x = continent)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha = 0.25, aes(col = continent)) +
  #stat_pointinterval() +
  facet_wrap(~source) +
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  coord_flip()
  

boxplot(results_full_raw$cor.sp_gwr_mean ~ sqrt(aa$diff*-1))




result_list_gbif <- lapply(families, 
                          subsampling.plants.fam.mapping, source = gbif_dist)

results_gbif <- rbindlist(result_list_gbif) %>% 
  mutate(source = "gbif")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif.txt")
fwrite(results_gbif, file=res_file)




result_list_gbif_inat <- lapply(families, 
                               subsampling.plants.fam.mapping, source = gbif_dist_inat)

results_gbif_inat <- rbindlist(result_list_gbif_inat) %>% 
  mutate(source = "gbif_inat")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif_inat.txt")
fwrite(results_gbif_inat, file=res_file)




result_list_gbif_herb <- lapply(families, 
                               subsampling.plants.fam.mapping, source = gbif_dist_herb)

results_gbif_herb <- rbindlist(result_list_gbif_herb) %>% 
  mutate(source = "gbif_herb")

res_file <- paste0("data/Samples_iteration_splm_rich_family_gbif_herb.txt")
fwrite(results_gbif_herb, file=res_file)




results_full_raw <- rbind(results_gbif, results_gbif_inat, results_gbif_herb)
res_file <- paste0("data/results_full_raw_familieswithin_mapping.txt")

fwrite(results_full_raw, file=res_file)


boxplot(log10(results_full_raw$dif) ~ results_full_raw$LEVEL1_NAM)
boxplot(results_full_raw$diff_prop ~ results_full_raw$LEVEL1_NAM)
hist(results_full_raw$dif)


aa <- results_full_raw %>% 
  filter(source == "gbif") %>% 
  mutate(diff = dif*-1) %>% 
  group_by(LEVEL1_NAM, family) %>% 
  summarise(mean = mean(diff), 
            median = median(diff), 
            sd = sd(diff), 
            median_prop = median(diff_prop))



aa <- results_full_raw %>% 
  filter(source == "gbif") %>% 
  mutate(diff = dif*-1) %>% 
  group_by(LEVEL1_NAM) %>% 
  summarise(mean = mean(diff), 
            median = median(diff), 
            sd = sd(diff), 
            median_prop = mean(diff_prop))



library(plotrix)
library(tmap)

map <- results_full_raw %>% 
  filter(!source == "wcvp") %>% 
  group_by(LEVEL3_COD, source) %>% 
  summarise(mean = mean(dif), 
            se = std.error(dif), 
            median = median(dif), 
            sd = sd(dif), 
            median_prop = median(diff_prop)) %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()

missing_prop <- tm_shape(map) +
  tm_fill(col = "median_prop", palette = "YlOrRd") +
  tm_borders() +
  tm_layout(
            asp = 2) + 
  tm_facets(by = "source",  free.scales = T, ncol = 3, nrow = 1) + 
  tm_format("World")
  

tmap_save(missing_prop, "median_prop_sp_missing_gbif_families.png", dpi = 600)


missing_se <- tm_shape(map) +
  tm_fill(col = "se", palette = "YlOrRd") +
  tm_borders() +
  tm_layout(
    asp = 2) + 
  tm_facets(by = "source",  free.scales = T, ncol = 3, nrow = 1) + 
  tm_format("World")



tmap_save(missing_se, "se_sp_missing_gbif_families.png", dpi = 600)


missing_mean <- tm_shape(map) +
  tm_fill(col = "mean", palette = "YlOrRd") +
  tm_borders() +
  tm_layout(
    asp = 2) +  
  tm_facets(by = "source",  free.scales = T, ncol = 3, nrow = 1) + 
  tm_format("World")



tmap_save(missing_mean, "mean_sp_missing_gbif_families.png", dpi = 600)



missing_median <- tm_shape(map) +
  tm_fill(col = "median", palette = "YlOrRd") +
  tm_borders() +
  tm_layout(
    asp = 2) + 
  tm_facets(by = "source",  free.scales = T, ncol = 3, nrow = 1) + 
  tm_format("World")



tmap_save(missing_median, "median_sp_missing_gbif_families.png", dpi = 600)



