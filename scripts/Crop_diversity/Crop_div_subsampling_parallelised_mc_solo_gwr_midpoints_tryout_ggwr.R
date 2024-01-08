
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

# Defining functions -----------------------------------------------------------
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
    mutate(sp = length(cumulative_namelist)) %>% 
    left_join(tdwg_3, by = c("LEVEL_COD" ="LEVEL3_COD")) %>% 
    st_sf()
  
  

  spar <- spautor(richness ~ richness_sample, data = rich_rel, spcov_type = "car")
  spar2 <- spautor(richness ~ 1, data = rich_rel, spcov_type = "car")
  
  glance(spar)
  
  rich_rel2 <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist)) %>% 
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
      mutate(sp = length(cumulative_namelist),
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
      mutate(sp = length(cumulative_namelist),
             richness_sample = richness_sample,
             richness = richness) %>% 
      st_sf() %>% 
      dplyr::select(geometry, richness, richness_sample) %>% 
      as("Spatial")
    
  
  rich_rel3 <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist)) %>% 
    left_join(tdwg_3, by = c("LEVEL_COD" ="LEVEL3_COD")) %>%  # tdwg_3_red
    #distinct(LEVEL_COD, richness, .keep_all = T) %>%
    #dplyr::select(geometry) %>% 
    st_sf()
  pred_sar <- splm(richness ~ richness_sample, rich_rel, "exponential")
  glance(pred_sar)
  pred_lm <- lm(richness ~ richness_sample, rich_rel)
  glance(pred_lm)
  
  bw <- bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
                   approach = "AIC",
                   adaptive = T)
  
  m.gwr <- gwr.basic(richness ~ richness_sample, data=rich_rel_shp,
                         adaptive = T,
                         bw = bw, 
                         cv = T, 
                         longlat = T)
  
  rich_rel3$preds <- predict(pred_sar, newdata = rich_rel3)
  rich_rel3$preds_lm <- predict(pred_lm, newdata = rich_rel3)
  
  se_preds <- do.call(cbind,predict(pred_sar, newdata = rich_rel3, se.fit = T))
  se_preds_lm <- do.call(cbind,predict(pred_lm, newdata = rich_rel3, se.fit = T))
  


  pred_gwr <- gwr.predict(richness ~ richness_sample, rich_rel_shp,
              rich_rel_shp2, bw = bw, kernel="bisquare",adaptive=T,
              longlat=T)
  
  pred_gwr_ext <- as.data.frame(pred_gwr$SDF$prediction)
  names(pred_gwr_ext) <- "preds_gwr"
  pred_gwr_ext$LEVEL_COD <- rich_rel_shp$LEVEL_COD
  

  rich_rel_pred <- rich_rel3 %>% 
    left_join(pred_gwr_ext, by = "LEVEL_COD") %>% 
    mutate(preds_diff = richness - preds, 
           preds_lm_diff = richness - preds_lm, 
           preds_gwr_diff = richness - preds_gwr) %>% 
    dplyr::select(richness, richness_sample, preds, preds_lm, preds_gwr, 
                  preds_diff, preds_lm_diff, preds_gwr_diff, 
                  LEVEL3_NAM, LEVEL1_COD.x) %>% 
    st_drop_geometry()
  
summary(rich_rel_pred)
  
  
rich_rel_pred_plot <- rich_rel3 %>% 
  left_join(pred_gwr_ext, by = "LEVEL_COD") %>% 
  mutate(preds_diff = richness - preds, 
         preds_lm_diff = richness - preds_lm, 
         preds_gwr_diff = richness - preds_gwr) %>% 
  dplyr::select(richness, richness_sample, preds, preds_lm, preds_gwr, 
                preds_diff, preds_lm_diff, preds_gwr_diff, 
                LEVEL3_NAM, LEVEL1_COD.x)

a <-  ggplot(rich_rel_pred_plot, aes(fill = preds)) +
    geom_sf(size = 2.5) +
    scale_fill_viridis_c() +
    theme_gray(base_size = 18) +
  ggtitle("preds")
  
b <-  ggplot(rich_rel_pred_plot, aes(fill = preds_lm)) +
    scale_fill_viridis_c() +
    geom_sf(size = 2.5) +
    theme_gray(base_size = 18) +
  ggtitle("preds_lm")

  
c <-  ggplot(rich_rel_pred_plot, aes(fill = preds_gwr)) +
    scale_fill_viridis_c() +
    geom_sf(size = 2.5) +
    theme_gray(base_size = 18) +
  ggtitle("preds_gwr")

  
d <-  ggplot(rich_rel_pred_plot, aes(fill = richness)) +
    scale_fill_viridis_c() +
    geom_sf(size = 2.5) +
    theme_gray(base_size = 18) + 
  ggtitle("richness")

e <- a + b + c + d
plot(e)  


library(colorspace)
a <-  ggplot(rich_rel_pred_plot, aes(fill = preds_diff)) +
  geom_sf(size = 2.5) +
  scale_fill_continuous_diverging(palette = "Blue-Red") +
  theme_minimal(base_size = 18) +
  ggtitle("preds")

b <-  ggplot(rich_rel_pred_plot, aes(fill = preds_lm_diff)) +
  scale_fill_continuous_diverging(palette = "Blue-Red") +
  geom_sf(size = 2.5) +
  theme_minimal(base_size = 18) +
  ggtitle("preds_lm")


c <-  ggplot(rich_rel_pred_plot, aes(fill = preds_gwr_diff)) +
  scale_fill_continuous_diverging(palette = "Blue-Red") +
  geom_sf(size = 2.5) +
  theme_minimal(base_size = 18) +
  ggtitle("preds_gwr")


d <-  ggplot(rich_rel_pred_plot, aes(fill = richness)) +
  scale_fill_viridis_c() +
  geom_sf(size = 2.5) +
  theme_minimal(base_size = 18) + 
  ggtitle("richness")

e <- a + b + c + d
plot(e)  
  
sulfmod <- splm(sulfate ~ 1, sulfate, spcov_type = "spherical")
  predict(sulfmod, newdata = sulfate_preds, se.fit = TRUE)
  
  sulfate_preds
  plot(rich_rel2)
  hist(sqrt(rich_rel2$richness))
  hist(sqrt(rich_rel2$richness_sample))
  plot(rich_rel2$richness ~ rich_rel2$richness_sample)
  
  plot(sqrt(rich_rel2$richness) ~ sqrt(rich_rel2$richness_sample))
  xx <- lm(sqrt(richness) ~ sqrt(richness_sample), data = rich_rel2)
  yy <- lm(richness ~richness_sample, data = rich_rel2)
  
  
  summary(xx)
  summary(yy)
  
  plot(xx)
  hist(xx$residuals)
  res <- as.numeric(xx$residuals)
  aa <- fitdistrplus::fitdist(res, "norm")
  aa <- fitdistrplus::fitdist(res, "norm")
  fit <- fitDist(res, k = 2, type = "realAll", trace = T, try.gamlss = TRUE)
  summary(fit)
  

  plot(aa)
  plot(bb)
  hist(yy$residuals)
  shapiro.test(xx$residuals)
  qqnorm(xx$residuals)
  # measuring correlation
  
  # overall
  #GWRbandwidth  <- gwr.sel(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel,
   #       coords=cbind(rich_rel$lon, rich_rel$lat), adapt = T, method = "aic")
  dmat <- gw.dist(coordinates(rich_rel2), focus=0, p=2, theta=0, longlat=F)
  
  GWRbandwidth2  <- bw.ggwr(richness ~ richness_sample, data= rich_rel2,
                            adaptive =   T, dMat = dmat, longlat = T, approach = "AIC")
  
  GWRbandwidth_gwr  <- bw.gwr(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel2,
                            adaptive =   T, longlat = F, approach = "CV")

  
  gwr_boot <- gwr.bootstrap(richness ~ richness_sample, data=rich_rel2,
                        adaptive = T, longlat = T, dMat = dmat, verbose = TRUE, 
                        approach = "AIC" ) 
  
  gwr_gen <- gwr.generalised(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel2, bw= GWRbandwidth_gwr, 
                            adaptive = T)    
  
  robust_gwr <- gwr.robust(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel2, bw= GWRbandwidth_gwr, 
             adaptive = T, )
  
  ggwr.model <- ggwr.basic(richness ~ richness_sample, data=rich_rel2,
                bw=GWRbandwidth2 , adaptive = T, maxiter = 100) 
  
  stats <- as.data.frame(ggwr.model$GW.diagnostic)
  
 gwr.model.selection(DeVar = "richness", InDeVars = c("richness_sample"), data=rich_rel2, bw=GWRbandwidth_gwr, 
                                   adaptive = T,  approach = "CV")
  
  modelling[[2]]
  
  
  gwr.model <- gwr.basic(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel,
                         bw=GWRbandwidth_gwr, 
                         adaptive = T) 
  
  gwr.model <- gwr(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel2, adapt=GWRbandwidth2 , hatmatrix=TRUE, se.fit=TRUE)
  

  bands <- gwr.model$bandwidth
  results <- as.data.frame(gwr.model$SDF) %>% 
    dplyr::select(- coords.x1, - coords.x2)  %>% 
    rename(richness_sample_pred = richness_sample) 
  
  gwr_map <- results %>% 
    cbind(rich_rel2, bands) %>% 
    left_join(tdwg_3_red, by = c("LEVEL_COD" ="LEVEL3_COD")) %>% 
    st_sf() 
  
  tm_shape(gwr_map) +
    tm_fill(c("richness","richness_sample", "richness_sample_pred", "richness_sample_se", "localR2" ,"bands"), palette = "inferno", style = "kmeans") +
    tm_layout(legend.position = c("left","bottom"), frame = F) +
      

    
  x1 <- qtm(gwr_map, fill = "richness", format = "World_wide",  style = tm_layout(frame = FALSE))
  x2 <-qtm(gwr_map, fill = "localR2")
  x3 <- qtm(gwr_map, fill = "bands")
  x4 <-qtm(gwr_map, fill = "bands")
  
  tmap_arrange(x1, x2, x3, x4, ncol = 2, nrow = 2)
  
  check <- gwr_map %>%  
    filter(LEVEL3_NAM == "Antarctica")

  plot(check)
  data("DubVoter")

  rich_rel3 <- rich_rel2
  rich_rel3 <- as(rich_rel2, "Spatial")

  
  mean(res_stats$Corr_richness.richness_sample)
  
  
  coordinates(rich_rel3)
  
  
  rich_rel_shp <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    left_join(midpoints_red, by = c("LEVEL_COD" = "LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  rich_rel_shp <- st_sf(rich_rel_shp)
  rich_rel_shp_sp <-  as(rich_rel_shp, "Spatial")
  
   bw <- bw.gwr(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel_shp_sp,
               approach = "AIC",
               adaptive = T) 
   
   m.gwr <- gwr.basic(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel_shp_sp,
                      adaptive = T,
                      bw = bw, 
                      cv = T)  
   m.gwr
   
   m.gwr <- gwr.basic(richness ~ richness_sample, data=rich_rel_shp_sp,
                      adaptive = T,
                      bw = bw, 
                      cv = T) 
   res <- as.data.frame(m.gwr$SDF@data)
   mean(res$Local_R2)
   median(res$Local_R2)
   summary(res)
   
   ss_stats <- gwss(data = rich_rel2, vars = c("richness", "richness_sample"),  adaptive = T, bw = GWRbandwidth_gwr)
   res_stats <- ss_stats$SDF
   mean(res_stats$Spearman_rho_richness.richness_sample)
   
   # GWR outouts
   names(m.gwr)
   tab.gwr <- rbind(apply(m.gwr$SDF@data[ , 1:12], 2, summary), coef(m))
   aa <- m.gwr$SDF@data

   res$LEVEL3_COD <- rich_rel_shp$LEVEL_COD
   res_map <-res %>% 
     left_join(tdwg_3_mod, by = "LEVEL3_COD", multiple = "all") %>%  
     st_as_sf()
     
   
   rownames(tab.gwr)[7] <- "Global"
   t(tab.gwr)
   
   tm_shape(res_map) +
     tm_fill(c("richness", "richness_sample", "richness_sample_SE"), palette = "inferno", style = "kmeans") +
     tm_layout(legend.position = c("left","bottom"), frame = F)
 
   
   
   tval = res_map %>%dplyr::select(all_of("richness_sample")) %>% st_drop_geometry()
   signif = tval < -1.96 | tval > 1.96
   # map the counties
   tm_shape(res_map) +
     tm_fill("PctBlack",midpoint = 0) + tm_style("col_blind")+
     tm_layout(legend.position = c("right","top"))+
     # now add the tvalues layer
     tm_shape(gwr_sf[signif,]) + tm_borders()
   # gwr.model <- ggwr(sqrt(richness) ~ sqrt(richness_sample), data=rich_rel,
  #                  coords=cbind(rich_rel$lon, rich_rel$lat), adapt=GWRbandwidth, family = "gaussian", type = "deviance") 
  
  results<-as.data.frame(gwr.model$SDF)
  cor <- as.data.frame(results$localR2)
  mean(cor$`results$localR2`)
  plot(results$localR2)
  summary(as.vector(st_distance(georgia)))

  
  # per continent 
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="pearson", data = ., exact =F)[[4]])  
  
  # bind all
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = length(cumulative_namelist)) %>% 
    left_join(spear, by = "id")
  
  if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
    print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  }
  
  # cumulative pattern
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
}

# Working directory ------------------------------------------------------------

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_3_mod <- st_read(dsn = "data/wgsrpd-master/level_3_richness")
tdwg_3_red <- st_read(dsn = "data/wgsrpd-master/level3_red") %>% 
  dplyr::select(LEVEL3_COD, geometry) %>%  
  st_make_valid()

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


  
#write.table(midpoints, "midpoints_coordinates.txt")
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
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints, by =c("LEVEL_COD" = "LEVEL3_COD"))

######




tm_shape(gwr.map) + 
  tm_fill(col = "richness", style="fixed", n=20, palette="Greens",
          breaks = c(seq(0, 20000, by=2000), Inf),   legend.show = FALSE) +
  tm_legend(outside=TRUE) +
  #tm_graticules(col = "grey", alpha = 0.7, labels.show = F) +
  #tm_compass(type = "4star", size = 1, position = c("left", "top")) +
  #'tm_scale_bar(breaks = c(0, 1500, 3000), text.size = 4) +
  tm_layout(frame = FALSE)
