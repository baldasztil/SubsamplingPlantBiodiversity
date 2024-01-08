
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
library(maptools)
#library(spgwr)
library(maps)
library(modelr)
library(broom)
library(tmap)
library(GWmodel)    # to undertake the GWR
library(fitdistrplus)
library(gamlss)
library(spdep)
library(spmodel)
library(patchwork)
library(rWCVP)
library(hrbrthemes)
library(ggsci)

# function ---------------------------------------------------------------------
subsampling.plants <- function(spec_n) {
  
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
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
}

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 



tdwg_3 <- rWCVP::wgsrpd3
aa <- st_make_valid(tdwg_3)

midpoints_raw <- st_centroid(aa)
midpoints_red <- midpoints_raw %>% 
  dplyr::select(LEVEL3_COD,geometry)

midpoints <- midpoints_raw %>% 
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(LEVEL3_NAM, LEVEL3_COD, lon,lat) %>% 
  left_join(tdwg_codes, by = "LEVEL3_COD" ) 



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

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints, by =c("LEVEL_COD" = "LEVEL3_COD"))

######


# Model comparison -----------------------------------------------------------

species <- subsampling.plants(701)
  
length(species$plant_name_id)
dist <- dist_native %>% 
  filter(plant_name_id %in% species$plant_name_id)

# richness patterns across brus
sample_rich_bru <- dist %>% 
  group_by(area_code_l3) %>% 
  summarise(richness_sample = n_distinct(plant_name_id)) %>% 
  rename(LEVEL_COD = area_code_l3)

rich_rel <- rich_overall_bru %>% 
  left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>% 
  left_join(tdwg_3, by = c("LEVEL_COD" ="LEVEL3_COD")) %>% 
  st_sf()

rich_rel2 <- rich_overall_bru %>% 
  left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>% 
  left_join(midpoints_red, by = c("LEVEL_COD" ="LEVEL3_COD")) %>%  # tdwg_3_red
  #distinct(LEVEL_COD, richness, .keep_all = T) %>%
  st_sf() %>% 
  as("Spatial") 


rich_rel_shp <- dist %>% 
  group_by(area_code_l3) %>% 
  summarise(richness_sample = n_distinct(plant_name_id)) %>% 
  rename(LEVEL_COD = area_code_l3) %>% 
  right_join(rich_overall_bru, by = "LEVEL_COD") %>%
  left_join(midpoints_red, by = c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id),
         richness_sample = richness_sample,
         richness = richness) %>% 
  st_sf() %>% 
  as("Spatial")

rich_rel_shp2 <- dist %>% 
  group_by(area_code_l3) %>% 
  summarise(richness_sample = n_distinct(plant_name_id)) %>% 
  rename(LEVEL_COD = area_code_l3) %>% 
  right_join(rich_overall_bru, by = "LEVEL_COD") %>%
  left_join(midpoints_red, by = c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id),
         richness_sample = richness_sample,
         richness = richness) %>% 
  st_sf() %>% 
  dplyr::select(geometry, richness, richness_sample) %>% 
  as("Spatial")


rich_rel3 <- rich_overall_bru %>% 
  left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>% 
  left_join(tdwg_3, by = c("LEVEL_COD" ="LEVEL3_COD")) %>%  # tdwg_3_red
  #distinct(LEVEL_COD, richness, .keep_all = T) %>%
  #dplyr::select(geometry) %>% 
  st_sf()

pred_point <- splm(richness ~ richness_sample, rich_rel, "exponential")
glance(pred_point)
pred_lm <- lm(richness ~ richness_sample, rich_rel)
glance(pred_lm)

pred_lm_RF <- splmRF(richness ~ richness_sample, rich_rel, sp_cov_type = "exponential")
glance(pred_lm_RF$splm)
summary(pred_lm_RF)



bw <- bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
             approach = "AIC",
             adaptive = T)

m.gwr <- gwr.basic(richness ~ richness_sample, data=rich_rel_shp,
                   adaptive = T,
                   bw = bw, 
                   cv = T, 
                   longlat = T)

spar <- spautor(richness ~ richness_sample, data = rich_rel, spcov_type = "car")
glance(spar)
#spar2 <- spautor(richness ~ 1, data = rich_rel, spcov_type = "car")
rich_rel3$preds_sar <- augment(spar, newdata = spar$newdata)$.fitted

rich_rel3$preds_RF <- predict(pred_lm_RF, newdata = rich_rel3)
rich_rel3$preds_point <- predict(pred_point, newdata = rich_rel3)
rich_rel3$preds_lm <- predict(pred_lm, newdata = rich_rel3)

se_preds_point <- do.call(cbind,predict(pred_point, newdata = rich_rel3, se.fit = T))
se_preds_lm <- do.call(cbind,predict(pred_lm, newdata = rich_rel3, se.fit = T))



pred_gwr <- gwr.predict(richness ~ richness_sample, rich_rel_shp,
                        rich_rel_shp2, bw = bw, kernel="bisquare",adaptive=T,
                        longlat=T)

pred_gwr_ext <- as.data.frame(pred_gwr$SDF$prediction)
names(pred_gwr_ext) <- "preds_gwr"
pred_gwr_ext$LEVEL_COD <- rich_rel_shp$LEVEL_COD


rich_rel_pred_plot <- rich_rel3 %>% 
  left_join(pred_gwr_ext, by = "LEVEL_COD") %>% 
  mutate(
    preds_RF_diff = richness - preds_RF,
    preds_sar_diff = richness - preds_sar,
    preds_point_diff = richness - preds_point, 
    preds_lm_diff = richness - preds_lm, 
    preds_gwr_diff = richness - preds_gwr) %>% 
  dplyr::select(richness, richness_sample, 
                preds_sar, preds_point, preds_gwr, preds_RF, preds_lm,
                preds_sar_diff, preds_point_diff, preds_gwr_diff, preds_RF_diff, preds_lm_diff, 
                LEVEL3_NAM.x, LEVEL1_COD.x) 

summary(rich_rel_pred_plot)

prediction_plots <- rich_rel_pred_plot %>% 
  dplyr::select(richness, 
                preds_sar, preds_point, preds_gwr,preds_RF, preds_lm,
                LEVEL3_NAM.x, LEVEL1_COD.x) %>% 
  pivot_longer(cols = c("richness", "preds_sar", "preds_point", "preds_gwr", "preds_RF", "preds_lm"), 
               values_to = "prediction", 
               names_to = "type")

diff_plots <- rich_rel_pred_plot %>% 
  dplyr::select(preds_sar_diff, preds_point_diff, preds_gwr_diff, preds_lm_diff, preds_RF_diff, 
                LEVEL3_NAM.x, LEVEL1_COD.x) %>% 
  pivot_longer(cols = c("preds_sar_diff", "preds_point_diff", "preds_gwr_diff", 
                        "preds_RF_diff", "preds_lm_diff"), 
               values_to = "prediction_diff", 
               names_to = "type") %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(type = fct_reorder(type, prediction_diff))

library(ggridges)
library(ggpubr)
library(see)
library(ggplot2)
library(ggsci)
library(patchwork)
library(hrbrthemes)
library(modelbased)
library(spmodel)
library(colorspace)

aa<- material_colors()
colour_scale <- as.vector(material_colors()[seq(1,18,3)])
colour_scale <- c(colour_scale[1:5], "grey")
see::palette_material()

pred_plotting <- ggplot(prediction_plots, aes(y = prediction, x = type, fill = type)) +
  geom_violindot(dots_size = 5, binwidth = 10) +
  scale_fill_manual(values = colour_scale) +
  #geom_violindot(data = richness, aes(y = prediction, x = type, fill = "midnightblue")) +
  theme_modern() +
  coord_flip()
pred_plotting

diff_plotting <- ggplot(diff_plots, aes(x = prediction_diff, y= type, fill = type)) +
  geom_density_ridges(alpha=0.8, panel_scaling = F, scale = 0.75) +
  #geom_violin(alpha = 0.7) + 
  scale_fill_manual(values = colour_scale[1:5]) +
  #scale_fill_viridis_d() +
  geom_vline(xintercept = 0, col = "red", alpha = 0.7, lwd = 1, linetype = 2) +
  theme_modern() 
  #coord_flip()
diff_plotting
pred_plotting +diff_plotting

ggsave("model_comparison_linear.svg", width = 15, height = 8)


a <-  ggplot(prediction_plots, aes(fill = prediction)) +
  geom_sf(size = 2.5) +
  scale_fill_viridis_c() +
  theme_modern(base_size = 12) +
  facet_wrap(~type)

b <-  ggplot(diff_plots, aes(fill = log10(prediction_diff))) +
  geom_sf(size = 2.5) +
  scale_fill_continuous_diverging("Blue-Red") +
  theme_modern(base_size = 12) +
  facet_wrap(~type)
