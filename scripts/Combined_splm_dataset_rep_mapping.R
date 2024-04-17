
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
library(ape)
library(phyloregion)
library(GWmodel)
library(feather)
library(vegan)
library(gmodels)


library(tmap)

# Defining functions -----------------------------------------------------------

spatial_lm <- function(x)  {
  model <- splm(richness ~ richness_sample, data = x,  
                spcov_type =  c("exponential"))
  model$dim <- x$dim
  return(model)
}
spatial_lm_grav <- function(x)  {
  model <- splm(richness ~ richness_sample, data = x,  
                spcov_type =  c("gravity"))
  model$dim <- x$dim
  return(model)
}
model_preds <-  function(model) {
  temp <- as.data.table(model$obdata)[, c("LEVEL1_NAM", "LEVEL3_COD", "richness", "richness_sample")] %>%
    mutate(preds =  fitted.values(model),  
           residuals =   residuals(model),
           dim = unique(model$dim)) 
  return(temp)
}


glance_splm <-  function(model) {
  temp <-  glance(model) %>% 
    mutate(continent = "GLOBAL", 
           model_id = 1, 
           coeff = model$coefficients$fixed[2],
           coeffse = sqrt(diag(vcov(model)))[2], 
           dim = unique(model$dim))
  
  return(temp)
}

splm.continents <-  function(x, rich_rel)  {
  
  continent_data <- rich_rel %>% 
    filter(LEVEL1_NAM == x) 
  
  output_glance_model <-try(lapply(split(continent_data,continent_data$dim),spatial_lm))
  
  if (class(output_glance_model)[1] =="try-error" & continent_data$sp[1] < 50000)  {
    output_glance <- global %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        continent = unique(continent_data$LEVEL1_NAM),
        sp = unique(continent_data$sp),
        model_id = 0, 
        coeff = 0,
        coeffse = 0, 
        dim = global$dim
      )
    return(output_glance)
  }
  
  if (class(output_glance_model)[1] =="try-error" & continent_data$sp[1] > 50000)  {
    
    output_glance_model <- try(lapply(split(continent_data,continent_data$dim),spatial_lm_grav))
    
    if (class(output_glance_model)[1] =="try-error") {
      
      output_glance <- global %>% 
        replace(!is.numeric(.) ==T, 0) %>% 
        mutate( 
          continent = unique(continent_data$LEVEL1_NAM),
          sp = unique(continent_data$sp),
          model_id = 0, 
          coeff = 1,
          coeffse = 0,
          pseudo.r.squared = 1, 
          dim = global$dim
        )
    } else
      output_glance <- rbindlist(lapply(output_glance_model, glance_splm)) %>% 
        mutate(continent = x, 
               sp = unique(continent_data$sp), 
               model_id = 2) 
    
    
    return(output_glance)
  }
  
  else {
    
    output_glance <- rbindlist(lapply(output_glance_model, glance_splm)) %>% 
      mutate(continent = x, 
             sp = unique(continent_data$sp), 
             model_id = 1)
    
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
             sp = nrow(species_sample), 
             dim = output_model$dim)
    return(error_data)
  }
  
  else {
    stats_total <-  cbind(stats, rich_rel_shp) %>% 
      dplyr::select(cor.sp_gwr = Spearman_rho_richness.richness_sample, LEVEL3_COD, 
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
    
    output_file <- rbind(continents, global) %>% 
      mutate(dim = unique(rich_rel_shp$dim))
    return(output_file)
  }
}

subsampling.plants <- function(x) {
  

  species_sample <- x
  
  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  # phylo richness
  subset_matrix <- dist.mat[ ,colnames(dist.mat) %in% dist$label]
  prune_tree <- match_phylo_comm(output_tree,subset_matrix)
  
  # doesnt work right now with old tree
  rich_rel_pd <- as.data.frame(PD(subset_matrix, prune_tree$phy)) %>% 
    mutate(LEVEL3_COD = rownames(.)) %>% 
    rename(richness_phy = `PD(subset_matrix, prune_tree$phy)`)
  
  table_data_sample <- table(dist$LEVEL3_COD, dist$growth_form)
  rich_rel_fun <- as.data.frame(diversity(table_data_sample, index = "shannon")) 
  rich_rel_fun$C <- rownames(rich_rel_fun)
  names(rich_rel_fun) <- c("richness_fun", "LEVEL3_COD")
  
  rich_rel_rich_tax <- dist %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(richness = length(unique(plant_name_id)), 
              richness_tax = length(unique(family)))
  
  richness_rel <- rich_rel_rich_tax %>% 
    left_join(rich_rel_pd, by = "LEVEL3_COD") %>% 
    left_join(rich_rel_fun, by = "LEVEL3_COD")
    
  
  rich_rel_long <- richness_rel %>% 
    pivot_longer(cols = starts_with("richness"),
                 names_to = "dim",
                 values_to = "richness_sample") %>% 
    mutate(richness_sample_scaled = scale(richness_sample)[,1]) 
  
  rich_rel <-  rich_rel_long %>% 
    right_join(rich_overall_bru, by = c("LEVEL3_COD", "dim")) %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = nrow(species_sample),
           richness_sample = richness_sample,
           richness = richness) %>% 
    st_sf() 
  
  
  rich_rel_shp <- lapply(split(rich_rel,rich_rel$dim), FUN = function(x)  as(x, "Spatial"))
  
  model <- try(lapply(split(rich_rel,rich_rel$dim),spatial_lm))
  
  if (class(model)[1] =="try-error")  {
    global_er <- test %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        sp = nrow(species_sample),
        model_id = -1) 
    return(global_er)
  } else 
    
  global_mod <- lapply(model, glance_splm) %>% 
    rbindlist() %>% 
    mutate(sp = length(species_sample$plant_name_id))
  
  
  rich_rel_preds  <- lapply(model, model_preds) %>% 
    rbindlist() 
  
  model_error <- suppressWarnings(cube(rich_rel_preds, j = list(sum(abs(richness - 
                                                                          preds))/
                                                                  length(richness),
                                                                mean(preds), sum(residuals^2), 
                                                                mean(residuals), 
                                                                sd(residuals), min(residuals), 
                                                                max(residuals)), 
                                       by = c("LEVEL1_NAM", "dim"))) %>% 
    filter(!is.na(dim))
  
  colnames(model_error) <- c("continent", "dim", "mae", "mpred","sumres2",
                             "mres", "sdres", "minres", "maxres")
  model_error[is.na(model_error)] <- "GLOBAL"

  model_comp <- rbindlist(lapply(continent_names_vec, splm.continents, 
                                 rich_rel = rich_rel )) %>%
    rbind(global_mod) 
  
  correlation <- lapply(rich_rel_shp, cor.gw) %>% 
    rbindlist()
  
  model_comp_out <- model_comp %>% 
    left_join(correlation, by = c("continent", "dim")) %>% 
    left_join(model_error, by = c("continent", "dim")) %>% 
    mutate_if(is.numeric,
              round,
              digits = 4) %>% 
    dplyr::select(-n, -p, -npar, -value, -AIC)
  
  #cumulative pattern
  if (nrow(species_sample) %in% seq(10,nrow(plantlist_names),50000)){
    tf1 <- paste0("/exports/eddie/scratch/s2117440/data/phylogenetic/checkpoints/checkpoint_phylo_",nrow(species_sample),".txt")
    fwrite(test, file= tf1) 
  }
  return(model_comp_out)
}

mapping.plants <- function(x) {
  
  
  species_sample <- x
  
  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  # phylo richness
  subset_matrix <- dist.mat[ ,colnames(dist.mat) %in% dist$label]
  prune_tree <- match_phylo_comm(output_tree,subset_matrix)
  
  # doesnt work right now with old tree
  rich_rel_pd <- as.data.frame(PD(subset_matrix, prune_tree$phy)) %>% 
    mutate(LEVEL3_COD = rownames(.)) %>% 
    rename(richness_phy = `PD(subset_matrix, prune_tree$phy)`)
  
  table_data_sample <- table(dist$LEVEL3_COD, dist$growth_form)
  rich_rel_fun <- as.data.frame(diversity(table_data_sample, index = "shannon")) 
  rich_rel_fun$C <- rownames(rich_rel_fun)
  names(rich_rel_fun) <- c("richness_fun", "LEVEL3_COD")
  
  rich_rel_rich_tax <- dist %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(richness = length(unique(plant_name_id)), 
              richness_tax = length(unique(family)))
  
  richness_rel <- rich_rel_rich_tax %>% 
    left_join(rich_rel_pd, by = "LEVEL3_COD") %>% 
    left_join(rich_rel_fun, by = "LEVEL3_COD")
  
  
  rich_rel_long <- richness_rel %>% 
    pivot_longer(cols = starts_with("richness"),
                 names_to = "dim",
                 values_to = "richness_sample") 
  return(rich_rel_long)
}
mapping.plants.random <- function(spec_n) {
  
  
  if (nrow(plantlist_names) > spec_n) {
    species_sample <- sample_n(plantlist_names, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names
  }
  
  
  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  # phylo richness
  subset_matrix <- dist.mat[ ,colnames(dist.mat) %in% dist$label]
  prune_tree <- match_phylo_comm(output_tree,subset_matrix)
  
  # doesnt work right now with old tree
  rich_rel_pd <- as.data.frame(PD(subset_matrix, prune_tree$phy)) %>% 
    mutate(LEVEL3_COD = rownames(.)) %>% 
    rename(richness_phy = `PD(subset_matrix, prune_tree$phy)`)
  
  table_data_sample <- table(dist$LEVEL3_COD, dist$growth_form)
  rich_rel_fun <- as.data.frame(diversity(table_data_sample, index = "shannon")) 
  rich_rel_fun$C <- rownames(rich_rel_fun)
  names(rich_rel_fun) <- c("richness_fun", "LEVEL3_COD")
  
  rich_rel_rich_tax <- dist %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(richness = length(unique(plant_name_id)), 
              richness_tax = length(unique(family)))
  
  richness_rel <- rich_rel_rich_tax %>% 
    left_join(rich_rel_pd, by = "LEVEL3_COD") %>% 
    left_join(rich_rel_fun, by = "LEVEL3_COD")
  
  
  rich_rel_long <- richness_rel %>% 
    pivot_longer(cols = starts_with("richness"),
                 names_to = "dim",
                 values_to = "richness_sample") 
  return(rich_rel_long)
}


# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full_raw <- fread("data/wcvp_accepted_merged.txt")

output_tree <- read.tree("data/output_tree_20012024.tre")

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")

midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

growth <- read.csv("data/growths_rf.csv") %>% 
  dplyr::select(growth_form = apd_plant_growth_form, plant_name_id)

plants_full <- plants_full_raw %>% 
  left_join(growth, by = "plant_name_id")

dist_native_calc <- dist_native %>% 
  left_join(plants_full, by = "plant_name_id")

table_data <- table(dist_native_calc$area_code_l3, dist_native_calc$growth_form)


shannon_div <- as.data.frame(diversity(table_data, index = "shannon")) 
shannon_div$C <- rownames(shannon_div)
names(shannon_div) <- c("richness_fun", "LEVEL3_COD")

richness_patterns_fun <- shannon_div 

continent_names_vec <- unique(continent_names$LEVEL1_NAM)

# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)

plantlist_dist_phylo_growth <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, taxon_name, LEVEL3_COD = area_code_l3, growth_form, family) %>% 
  mutate(label = gsub(" ", "_", taxon_name))

dist.mat <- long2sparse(plantlist_dist_phylo_growth, grids = "LEVEL3_COD", species = "label")
prune_tree <- match_phylo_comm(output_tree,dist.mat)


richness_patterns_phy <- as.data.frame(PD(dist.mat, prune_tree$phy), 
                                   col.names =c("richness")) %>% 
  mutate(LEVEL3_COD = rownames(.)) %>% 
  rename(richness_phy = `PD(dist.mat, prune_tree$phy)`) 
  

richness_patterns_rich <- dist_native_calc %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = length(unique(plant_name_id))) %>% 
  rename(LEVEL3_COD = area_code_l3)

richness_patterns_tax <- dist_native_calc %>% 
  group_by(area_code_l3) %>% 
  summarise(richness_tax = length(unique(family))) %>% 
  rename(LEVEL3_COD = area_code_l3)


richness_patterns <-richness_patterns_rich %>% 
  left_join(richness_patterns_tax, by = "LEVEL3_COD") %>% 
  left_join(richness_patterns_phy, by = "LEVEL3_COD") %>% 
  left_join(richness_patterns_fun, by = "LEVEL3_COD") %>% 
  left_join(continent_names, by = "LEVEL3_COD") 

rich_patterns_long <- richness_patterns %>% 
  pivot_longer(cols = starts_with("richness"),
               names_to = "dim",
               values_to = "richness") %>% 
  mutate(richness_scaled = scale(richness)[,1])


rich_overall_bru <- rich_patterns_long


rich_overall_bru_mid <- rich_patterns_long %>% 
  left_join(midpoints_red, by = "LEVEL3_COD")  %>% 
  mutate(richness2 = richness) %>% 
  st_as_sf()

rich_overall_bru_shp <- rich_patterns_long %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()

rich_overall_bru_mid_shp <- rich_overall_bru_mid %>% 
  mutate(richness_sample = richness) %>% 
  as("Spatial")

rich_overall_bru_mid_test <- split(rich_overall_bru_mid_shp,rich_overall_bru_mid_shp$dim)


bw_global <-bw.gwr(richness ~ richness2, rich_overall_bru_mid_shp, 
                   approach = "AIC",
                   adaptive = T)


setDT(plantlist_names)
setDT(dist_native)
setDT(rich_overall_bru)
setDT(midpoints_red)


model_global <- lapply(split(rich_overall_bru_shp,rich_overall_bru_shp$dim), FUN = function(x) {
  model <- splm(richness ~ 1, data = x,  
                spcov_type =  c("exponential"))
  model$dim <- x$dim
  return(model)
})
 
global <- lapply(model_global, glance_splm) %>% 
  rbindlist()



output_model <- lapply(rich_overall_bru_mid_test, cor.gw) %>% 
  rbindlist()

test <-  subsampling.plants(plantlist_names[1:2639])



# compare models ---------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)

gbif <- fread("data/dist_gbif.txt") 
srli <- read.csv("data/red/srli/srli_full.csv")
redlist <- read.csv("data/red/redlist_full_syn.csv")

gbif_names <- plants_full %>% 
  filter(plant_name_id %in% gbif$plant_name_id)

srli_names <- plants_full %>% 
  filter(taxon_name %in% srli$taxon_name)

redlist_names <- plants_full %>% 
  filter(taxon_name %in% redlist$taxon_name)

gbif_rep <- subsampling.plants(gbif_names)
srli_rep <- subsampling.plants(srli_names)
redlist_rep <- subsampling.plants(redlist_names)


test2 <-  subsampling.plants(gbif_names[1:2639])



test$diff <- test$cor.sp_gwr_mean - test2$cor.sp_gwr_mean
test$diff_r <- test$pseudo.r.squared - test2$pseudo.r.squared



gbif_mapping <- mapping.plants(gbif_names)

rep <- rep(sample_n(plantlist_names, nrow(srli_names)), 100)



srli_mapping <- mapping.plants(srli_names) 

srli_mapping_random <- lapply(rep(nrow(srli_names),100), #nrow(plantlist_names)
                              mapping.plants.random) %>% 
  rbindlist() %>% 
  group_by(LEVEL3_COD, dim) %>%
  summarise(mean = ci(richness_sample)[1], 
            sd = sd(richness_sample),
            lower_ci =  ci(richness_sample)[2], 
            upper_ci =  ci(richness_sample)[3], 
            se =  ci(richness_sample)[4])

srli_diff <-  srli_mapping_random %>% 
  left_join(srli_mapping, by = c("dim", "LEVEL3_COD")) %>% 
  mutate(diff = richness_sample - mean, 
         bias = diff / sd
         ) %>% 
  mutate(significance = case_when(
    richness_sample > lower_ci  & richness_sample < upper_ci ~ "no_bias", 
    richness_sample < lower_ci ~ "under_estimation",
    richness_sample > upper_ci~ "over_estimation", 
                                TRUE ~ "not applicable"), 
    
    bias_level = case_when(abs(bias) > 5 ~ "high", 
                                abs(bias) > 2  ~ "moderate",
                                abs(bias) < 2 ~ "low", 
                                TRUE ~ "not applicable")) %>% 
  right_join(tdwg_3, by = c("LEVEL3_COD"), keep = T) %>% 
  st_as_sf()


tm_shape(srli_diff) +
  tm_fill(col = "significance") +
  tm_borders() +
  tm_facets(by = "dim",  free.scales = TRUE)


redlist_mapping <- mapping.plants(redlist_names) 


redlist_mapping_random <- lapply(rep(nrow(redlist_names),100), #nrow(plantlist_names)
                              mapping.plants.random) %>% 
  rbindlist() %>% 
  group_by(LEVEL3_COD, dim) %>% 
  group_by(LEVEL3_COD, dim) %>%
  summarise(mean = ci(richness_sample)[1], 
            sd = sd(richness_sample),
            lower_ci =  ci(richness_sample)[2], 
            upper_ci =  ci(richness_sample)[3], 
            se =  ci(richness_sample)[4])

redlist_diff <-  redlist_mapping_random %>% 
  left_join(redlist_mapping, by = c("dim", "LEVEL3_COD")) %>% 
  mutate(diff = richness_sample - mean, 
         bias = diff / sd
  ) %>% 
  mutate(significance = case_when(
    richness_sample > upper_ci~ "over_estimation", 
    richness_sample > lower_ci  & richness_sample < upper_ci ~ "no_bias", 
    richness_sample < lower_ci ~ "under_estimation",
    TRUE ~ "not applicable"), 
    
    bias_level = case_when(abs(bias) > 5 ~ "high", 
                           abs(bias) > 2  ~ "moderate",
                           abs(bias) < 2 ~ "low", 
                           TRUE ~ "not applicable")) %>% 
  right_join(tdwg_3, by = c("LEVEL3_COD"), keep = T) %>% 
  st_as_sf()


sig_red <- tm_shape(redlist_diff) +
  tm_fill(col = "significance", palette = "Accent") +
  tm_borders() +
  tm_facets(by = "dim",  free.scales = F, ncol = 2, nrow = 2) +
  tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
            asp = 2) +
  tm_format("World")

tmap_save(sig_red, "significance_red_maps.png", dpi = 600)


bias_red <- tm_shape(redlist_diff) +
  tm_fill(col = "bias", palette = "PuOr") +
  tm_borders() +
  tm_facets(by = "dim",  free.scales = T, ncol = 2, nrow = 2) +
  #tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
   #         asp = 2) +
  tm_format("World")

tmap_save(bias_red, "bias_level_red_maps.png", dpi = 600)

diff_red <- tm_shape(redlist_diff) +
  tm_fill(col = "diff", palette = "PRGn") +
  tm_borders() +
  tm_facets(by = "dim",  free.scales = T, ncol = 2, nrow = 2) +
  #tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
  #         asp = 2) +
  tm_format("World")

tmap_save(diff_red, "diff_level_red_maps.png", dpi = 600)





gbif_mapping <- mapping.plants(gbif_names) 

gbif_mapping %>%  
  group_by(dim) %>% 
  summarise(sum = sum(richness_sample))

gbif_mapping_random_all <- lapply(rep(nrow(gbif_names),100), #nrow(plantlist_names)
                                 mapping.plants.random) %>% 
  rbindlist()

gbif_mapping_random <- gbif_mapping_random_all%>% 
  group_by(LEVEL3_COD, dim) %>%
  summarise(mean = ci(richness_sample)[1], 
            sd = sd(richness_sample),
            lower_ci =  ci(richness_sample)[2], 
            upper_ci =  ci(richness_sample)[3], 
            se =  ci(richness_sample)[4])

gbif_mapping_random %>%  
  group_by(dim) %>% 
  summarise(sum = sum(mean))



gbif_diff <-  gbif_mapping_random %>% 
  left_join(gbif_mapping, by = c("dim", "LEVEL3_COD")) %>% 
  mutate(diff = richness_sample - mean, 
         bias = diff / sd
  ) %>% 
  mutate(significance = case_when(
    richness_sample > lower_ci  & richness_sample < upper_ci ~ "no_bias", 
    richness_sample < lower_ci ~ "under_estimation",
    richness_sample > upper_ci~ "over_estimation", 
    TRUE ~ "not applicable"), 
    
    bias_level = case_when(abs(bias) > 5 ~ "high", 
                           abs(bias) > 2  ~ "moderate",
                           abs(bias) < 2 ~ "low", 
                           TRUE ~ "not applicable")) %>% 
  right_join(tdwg_3, by = c("LEVEL3_COD"), keep = T) %>% 
  st_as_sf()

sum(gbif_diff$richness_sample)

rich_diff <- gbif_diff %>% 
  st_drop_geometry() %>% 
  filter(dim == "richness") %>% 
  left_join(richness_patterns_rich, by = c("LEVEL3_COD.x" = "LEVEL3_COD")) %>% 
  dplyr::select(LEVEL3_COD.x, richness, richness_sample, mean, sd, se, diff, bias, LEVEL1_NAM)

cpp <- gbif_mapping_random_all %>%  
  filter(dim == "richness") %>% 
  filter(LEVEL3_COD == "CPP")

sig_gbif <- tm_shape(gbif_diff) +
  tm_fill(col = "significance", palette = "Accent") +
  tm_borders() +
  tm_facets(by = "dim",  free.scales = F, ncol = 2, nrow = 2) +
  tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
            asp = 2) +
  tm_format("World")

tmap_save(sig_gbif, "significance_gbif_maps.png", dpi = 600)



bias_gbif <- tm_shape(gbif_diff) +
  tm_fill(col = "bias", palette = "PuOr") +
  tm_borders() +
  tm_facets(by = "dim",  free.scales = T, ncol = 2, nrow = 2) +
  #tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
  #         asp = 2) +
  tm_format("World")

tmap_save(bias_gbif, "bias_level_gbif_maps.png", dpi = 600)

diff_gbif <- tm_shape(gbif_diff) +
  tm_fill(col = "diff", palette = "PRGn") +
  tm_borders() +
  tm_facets(by = "dim",  free.scales = T, ncol = 2, nrow = 2) +
  #tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
  #         asp = 2) +
  tm_format("World")

tmap_save(diff_gbif, "diff_level_gbif_maps.png", dpi = 600)


under <- gbif_diff %>% 
  filter(significance == "under_estimation")

boxplot(gbif_diff$bias)

##### create null model for each country with quantiles of 100 randomisations 
# see if value falls outside 