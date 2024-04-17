
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
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
library(modelbased)

library(rWCVP)

library(ggridges)
library(ggpubr)
library(see)
library(ggsci)
library(patchwork)
library(hrbrthemes)
library(spmodel)
library(colorspace)




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
dist_native <- dist_native <- fread("data/dist_native.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
richness_patterns <- fread("data/richness_patterns.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")


# manipulate data --------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id,continent_code_l1,area_code_l3)

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)


rich_overall_bru <- richness_patterns

rich_overall_bru_mid <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by = "LEVEL3_COD") 

rich_overall_bru_shp <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") 

# Analysis ---------------------------------------------------------------------

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
  rename(LEVEL3_COD = area_code_l3)

rich_rel_shp <- rich_overall_bru_shp %>% 
  left_join(sample_rich_bru, by = "LEVEL3_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>%
  st_sf()

rich_rel_shp2 <- rich_overall_bru_shp %>% 
  left_join(sample_rich_bru, by = "LEVEL3_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>%
  st_sf() %>% 
  as("Spatial")

rich_rel_point <- rich_overall_bru_mid %>% 
  left_join(sample_rich_bru, by = "LEVEL3_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>% 
  #distinct(LEVEL3_COD, richness, .keep_all = T) %>%
  st_sf() 

rich_rel_point2 <- rich_overall_bru_mid %>% 
  left_join(sample_rich_bru, by = "LEVEL3_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>% 
  #distinct(LEVEL3_COD, richness, .keep_all = T) %>%
  st_sf() %>% 
  as("Spatial")

rich_rel3 <- rich_overall_bru_shp %>% 
  left_join(sample_rich_bru, by = "LEVEL3_COD") %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(species$plant_name_id)) %>% 
  st_as_sf()

pred_point <- splm(richness ~ richness_sample, rich_rel_point, "gravity")
glance(pred_point)

#pred_point_gen <- spglm(richness ~ richness_sample, family = "poisson", rich_rel_point, spcov_type =  "gravity")
#glance(pred_point_gen)

pred_lm <- lm(richness ~ richness_sample, rich_rel3)
glance(pred_lm)


spar <- spautor(richness ~ richness_sample, data = rich_rel_shp, spcov_type = "car")
glance(spar)

pred_lm_RF <- splmRF(richness ~ richness_sample, rich_rel_point, sp_cov_type = "gravity")
glance(pred_lm_RF$splm)
pred_lm_RF_model <- pred_lm_RF$splm



bw <- bw.gwr(richness ~ richness_sample, data=rich_rel_shp2,
             approach = "AIC",
             adaptive = T)

m.gwr <- gwr.basic(richness ~ richness_sample, data=rich_rel_shp2,
                   adaptive = T,
                   bw = bw, 
                   cv = T, 
                   longlat = T)

#predict(spar, newdata = rich_rel3, se.fit = T)

#spar2 <- spautor(richness ~ 1, data = rich_rel, spcov_type = "car")
rich_rel3$preds_sar <- augment(spar, newdata = spar$newdata)$.fitted
rich_rel3$preds_RF <- predict(pred_lm_RF, newdata = rich_rel3)
rich_rel3$preds_point <- predict(pred_point, newdata = rich_rel3)
rich_rel3$preds_lm <- predict(pred_lm, newdata = rich_rel3)

pred_gwr <- gwr.predict(richness ~ richness_sample, rich_rel_shp2,
                        rich_rel_shp2, bw = bw, kernel="bisquare",adaptive=T,
                        longlat=T)

pred_gwr_ext <- as.data.frame(pred_gwr$SDF$prediction)
names(pred_gwr_ext) <- "preds_gwr"
pred_gwr_ext$LEVEL3_COD <- rich_rel_shp$LEVEL3_COD


se_preds_point <- do.call(cbind, predict(pred_point, newdata = rich_rel3, se.fit = T))
se_preds_lm <- do.call(cbind, predict(pred_lm, newdata = rich_rel3, se.fit = T))
se_preds_RF <- do.call(cbind, predict(pred_lm_RF_model, newdata = rich_rel3, se.fit = T))

se <- as.data.frame(cbind(se_preds_point, se_preds_RF, se_preds_lm))


rich_rel_pred_plot <- rich_rel3 %>% 
  left_join(pred_gwr_ext, by = "LEVEL3_COD") %>% 
  mutate(
    preds_RF_diff = richness - preds_RF,
    preds_sar_diff = richness - preds_sar,
    preds_point_diff = richness - preds_point, 
    preds_lm_diff = richness - preds_lm, 
    preds_gwr_diff = richness - preds_gwr) %>% 
  dplyr::select(richness, richness_sample, 
                preds_sar, preds_point, preds_gwr, preds_RF, preds_lm,
                preds_sar_diff, preds_point_diff, preds_gwr_diff, preds_RF_diff, preds_lm_diff, 
                LEVEL3_NAM, LEVEL1_NAM)

summary(rich_rel_pred_plot)


prediction_plots <- rich_rel_pred_plot %>% 
  dplyr::select(richness, 
                preds_sar, preds_point, preds_gwr,preds_RF, preds_lm,
                LEVEL3_NAM, LEVEL1_NAM) %>% 
  pivot_longer(cols = c("richness", "preds_sar", "preds_point", "preds_gwr", "preds_RF", "preds_lm"), 
               values_to = "prediction", 
               names_to = "type") %>% 
  mutate(type = as.factor(type)) 

prediction_plots$type <- factor(prediction_plots$type, levels=c( 
                                    "preds_lm",
                                    "preds_gwr",
                                    "preds_sar",
                                    "preds_point", 
                                    "preds_RF",
                                    "richness"))


diff_plots <- rich_rel_pred_plot %>% 
  dplyr::select(preds_sar_diff, preds_point_diff, preds_gwr_diff, preds_lm_diff, preds_RF_diff, 
                LEVEL3_NAM, LEVEL1_NAM) %>% 
  pivot_longer(cols = c("preds_sar_diff", "preds_point_diff", "preds_gwr_diff", 
                        "preds_RF_diff", "preds_lm_diff"), 
               values_to = "prediction_diff", 
               names_to = "type") %>% 
  mutate(type = as.factor(type)) 

diff_plots$type <- factor(diff_plots$type, levels=c( 
  "preds_lm_diff",
  "preds_gwr_diff",
  "preds_sar_diff",
  "preds_point_diff", 
  "preds_RF_diff"))
diff_plots$type 

colour_scale <- as.vector(material_colors()[seq(1,18,3)])
colour_scale <- c(colour_scale[1:5], "grey")

pred_plotting <- ggplot(prediction_plots, aes(y = prediction, x = type, fill = type)) +
  geom_violindot(dots_size = 5, binwidth = 10, lwd = 0.75) +
  scale_fill_manual(values = colour_scale) +
  #geom_violindot(data = richness, aes(y = prediction, x = type, fill = "midnightblue")) +
  theme_bw() +
  coord_flip()
pred_plotting

diff_plots$type
diff_plotting <- ggplot(diff_plots, aes(x = prediction_diff, y= type, fill = type)) +
  #geom_density_ridges(alpha=0.8, panel_scaling = F, scale = 0.75) +
  geom_vline(xintercept = 0, col = "black", alpha = 1, lwd = 0.5, linetype = 1) +
  geom_violin(alpha = 0.9, lwd = NA, show.legend = F) + 
  scale_fill_manual(values = colour_scale[1:5]) +
  xlim(-3500,3500) +
  #scale_fill_viridis_d() +
  theme_bw() 
  #coord_flip()
diff_plotting
pred_plotting + diff_plotting + plot_layout(guides = "collect")

max(diff_plots$prediction_diff)
min(diff_plots$prediction_diff)


ggsave("model_comparison_linear.svg", width = 15, height = 8)


a <-  ggplot(prediction_plots, aes(fill = prediction)) +
  geom_sf(size = 2.5) +
  scale_fill_viridis_c() +
  theme_modern(base_size = 12) +
  facet_wrap(~type)

ggsave("prediction_maps.svg", width = 15, height = 10)

b <-  ggplot(diff_plots, aes(fill = prediction_diff)) +
  geom_sf(size = 2.5) +
  scale_fill_continuous_diverging("Blue-Red") +
  theme_modern(base_size = 12) +
  facet_wrap(~type)

ggsave("difference_mapsUnlog.svg", width = 15, height = 10)

