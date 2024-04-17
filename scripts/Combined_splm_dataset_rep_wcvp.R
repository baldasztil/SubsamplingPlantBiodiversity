
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

library(datawizard)

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
  mutate(richness_scaled = scale(richness)[,1], 
         dataset = "wcvp_all")


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

spec_n <- 82841

sample1 <-  lapply(rep(spec_n,100), #nrow(plantlist_names)
                 mapping.plants.random) %>% 
  rbindlist() %>% 
  mutate(dataset = paste0("sample_size_", spec_n))


spec_n <- 2639

sample2 <-  lapply(rep(spec_n,100), #nrow(plantlist_names)
                 mapping.plants.random) %>% 
  rbindlist() %>% 
  mutate(dataset = paste0("sample_size_", spec_n))

sample <- rbind(sample1, sample2)

overall_joined <- sample %>% 
  group_by(dim) %>% 
  mutate(richness = richness_sample, 
         prop_richness = richness / spec_n) %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_drop_geometry() %>% 
  dplyr::select(LEVEL3_COD, dim,  prop_richness, richness, dataset, continent = LEVEL1_NAM)

combined_rich <- rich_overall_bru %>% 
  dplyr::select(LEVEL3_COD, dim, richness, dataset, continent = LEVEL1_NAM) %>% 
  mutate(prop_richness  = richness / 349113) %>%  
  rbind(overall_joined) 




comparison <- 
  ggplot(combined_rich, aes(x = reorder(continent,richness), y =richness, fill = dataset)) +
  #geom_jitter(color = "black", alpha = 0.05) +
  labs(x = "Continent", y = "richness") +
  scale_fill_viridis_d() +
  stat_boxplot(geom = "errorbar", width = 0.5, alpha = 0.95) +
  facet_wrap(vars(dim, dataset), scales = "free", ncol = 2) + 
  geom_boxplot(notch=F,alpha=0.95, width =0.5, color = "grey10", outlier.alpha = 0.05, outlier.size = 0.5) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw()
comparison
ggsave(paste0("comparison_wcvp_sample_", spec_n, ".png"), 
       comparison, width = 10, height = 14, dpi = 600)

summary_stats  <- combined_rich %>% 
  group_by(LEVEL3_COD, dim, dataset) %>%
  summarise(mean = ci(richness)[1], 
            sd = sd(richness),
            lower_ci =  ci(richness)[2], 
            upper_ci =  ci(richness)[3], 
            se =  ci(richness)[4], 
            coef = sd / mean * 100) %>% 
  ungroup() %>% 
  group_by(dim, dataset) %>% 
  mutate(ranked = 369 - ranktransform(mean)) %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()


overall  <- rich_overall_bru %>% 
  group_by(dim, ) %>% 
  mutate(ranked = 369 - ranktransform(mean)) %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()

rich_sample <- tm_shape(summary_stats) +
  tm_fill(col = "mean", palette = "YlOrRd") +
  tm_borders() +
  tm_facets(by = c("dim", "dataset"),  free.scales = T, ncol = 3, nrow = 4) +
  #tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
   #        asp = 2) +
  tm_format("World_wide")

tmap_save(rich_sample, paste0("richness_sample_maps_", ".png"),  dpi = 600)


ranked_sample <- tm_shape(summary_stats) +
  tm_fill(col = "ranked", palette = "-YlOrRd", n = 37) +
  tm_borders() +
  tm_facets(by = c("dim", "dataset"),  free.scales = F, ncol = 3, nrow = 4) +
  tm_layout(legend.outside = T, legend.position	= c("center", "TOP"), 
           asp = 2) +
  tm_format("World")

tmap_save(ranked_sample, paste0("ranked_richness_sample_maps_", ".png"),  dpi = 600)


##### create null model for each country with quantiles of 100 randomisations 
# see if value falls outside 

rank_overall <- rich_overall_bru %>% 
  dplyr::select(LEVEL3_COD, dim, richness, wcvp_rank) %>% 
  group_by(dim) %>% 
  mutate(wcvp_rank  = 369 - ranktransform(richness))


samples_stats <- summary_stats %>% 
  filter(!dataset == "wcvp_all") %>% 
  left_join(rank_overall, by = c("dim", "LEVEL3_COD")) %>% 
  mutate(rank_diff = ranked - wcvp_rank, 
         log_coef = log10(coef)) 

samples_stats$log_coef[!is.finite(samples_stats$log_coef)] <- 0
samples_stats$coef[samples_stats$coef > 100 ] <- 100


coef_sample <- tm_shape(samples_stats) +
  tm_fill(col = "coef", palette = "YlOrRd", n = 7, midpoint = NA) +
  tm_borders() +
  tm_facets(by = c("dim", "dataset"),  free.scales = F, ncol = 4, nrow = 2) +
  tm_layout(
           asp = 2) +
  tm_format("World")

tmap_save(coef_sample, paste0("coef_sample_maps_",".png"),  dpi = 600)


rank_sample <- tm_shape(samples_stats) +
  tm_fill(col = "rank_diff", palette = "PRGn", n = 10, midpoint = 0) +
  tm_borders() +
  tm_facets(by = c("dim", "dataset"),  free.scales = F, ncol = 4, nrow = 2) +
  tm_layout( 
           asp = 2) +
  tm_format("World")

tmap_save(rank_sample, paste0("rank_diff_sample_maps_", ".png"),  dpi = 600)

full_comp <- tmap_arrange(rank_sample, coef_sample)

tmap_save(full_comp, paste0("rank_diff_log_coef_sample_maps", ".png"),  dpi = 600, 
         asp = 2)

