library(DEoptim)
library(tidyverse)
library(data.table)
library(sf)
library(foreach)
library(doSNOW)
library(GWmodel)
library(ModelMap)
library(gstat)

output <- function(species=NULL,best_corr=NULL)
{
  me <- data.frame(
    Species = species,
    best_corr = best_corr
  )
  
  return(me)
}

objective <- function(subset_indices) {
  species <- slice_sample(n = round(subset_indices),plantlist_names)
  dist <- plantlist_dist %>%
    filter(plant_name_id %in% species$plant_name_id)
  
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
 
   rich_rel_shp <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    left_join(midpoints_red, by = c("LEVEL_COD" = "LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(species)) %>% 
    st_sf() %>% 
    as("Spatial")
   
   write.table(rich_rel_shp, "rich_rel_shp.csv", col.names = TRUE, row.names = FALSE,  sep = ",")
   get.test(proportion.test = 0.5, "rich_rel_shp.csv", seed = 42, folder = getwd(),  qdata.trainfn = "gwr_calib.csv", qdata.testfn = "gwr_valid.csv")
   
   gwr_calib <- read.table("gwr_calib.csv", header = TRUE, sep = ",") %>% 
     mutate(long = as.numeric(coords.x1), 
            lat = as.numeric(coords.x2))
   
   rich_rel_shp_calib_pred <- SpatialPointsDataFrame(gwr_calib[, 12:13], data = gwr_calib) 
   
   gwr_valid <- read.table("gwr_valid.csv", header = TRUE, sep = ",") %>% 
     mutate(long = as.numeric(coords.x1), 
            lat = as.numeric(coords.x2))
   
   rich_rel_shp_valid_pred <- SpatialPointsDataFrame(gwr_valid[, 12:13], data = gwr_valid) 
   
   dm.valid <- gw.dist(dp.locat = coordinates(rich_rel_shp_calib_pred),  rp.locat = coordinates(rich_rel_shp_valid_pred))
   dm.calib <- gw.dist(dp.locat = coordinates(rich_rel_shp_calib_pred))
   
   gwr.bw.cv <- bw.gwr(richness ~ richness_sample,data = rich_rel_shp_calib_pred,
                       approach = "CV", kernel = "bisquare", adaptive = TRUE, dMat = dm.calib)
   
   gwr.pred <- gwr.predict(richness ~ richness_sample,data = rich_rel_shp_calib_pred, predictdata = rich_rel_shp_calib_pred ,
                           bw = gwr.bw.cv, kernel = "bisquare", adaptive = TRUE, dMat1 = dm.valid, dMat2 = dm.calib)
   
   ols.pred.gwmodel <- gwr.predict(richness ~ richness_sample, data = rich_rel_shp_calib_pred, predictdata = rich_rel_shp_valid_pred,
                                   bw = 368,  kernel = "bisquare", adaptive = TRUE, dMat1 = dm.valid, dMat2 = dm.calib)
   
   
   RMSPE.gwr <- (mean((rich_rel_shp_valid_pred$richness - gwr.pred$SDF$prediction)^2))^0.5 
   MAPE.gwr <- mean(abs(rich_rel_shp_valid_pred$richness - gwr.pred$SDF$prediction)) 
   zscore.gwr <- (rich_rel_shp_valid_pred$richness - gwr.pred$SDF$prediction)/ + (gwr.pred$SDF$prediction_var)^0.5 
   MeanZ.gwr <- mean(zscore.gwr) 
   SDZ.gwr <- (var(zscore.gwr))^0.5
   
   bw <- try(bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
                    approach = "AIC",
                    adaptive = T)) # the problem is that it gives me an NA for species n < 5
   # can fix with capture output or use global bandwith
   #sink()
   if (!is.numeric(bw)) {
     bw <- bw_global
   }
   
   m.gwr <- try(gwr.basic(richness ~ richness_sample, data=rich_rel_shp,
                          adaptive = T,
                          bw = bw, 
                          cv = T)) 
   rich_rel_sum <- rich_rel_shp
   sum_stats<- gwss(rich_rel_shp,rich_rel_sum, c("richness", "richness_sample"), adaptive = T, bw = bw )
   stats_sum <- as.data.frame(sum_stats$SDF)
   
   if (class(m.gwr)=="try-error") {
     error_data <- 10000
     return(error_data)
   }
   
   else {
     # m.gwr <- gwr(richness ~ richness_sample, data=rich_rel_shp,
     #             adapt = bw_spgwr, 
     #                  gweight = gwr.bisquare) 
     
    res_extract <- as.data.frame(m.gwr$SDF@data)
    cor.gwr <- mean(res_extract$Local_R2)
  

  return(-cor.gwr)  # Return negative correlation (since DEoptim minimizes)
  #return(ifelse(correlation >= 0.95, correlation, 0))
   }
}


# Working directory ------------------------------------------------------------
# Import data ------------------------------------------------------------------
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

richness_patterns <- fread("data/richness_patterns.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")

midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints") %>% 
  dplyr::select(LEVEL3_COD,geometry)


# creating objects for the function

plantlist_names <- plants_full %>% 
 dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native %>% 
 dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)


richness_patterns_con <- plantlist_dist %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- plantlist_dist %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- plantlist_dist %>% 
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

rich_overall_bru_bw <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  mutate(richness2 = richness) %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  st_sf() %>% 
  as("Spatial") 


bw_global <- bw.gwr(richness ~ richness2, data=rich_overall_bru_bw,
                    approach = "AIC",
                    adaptive = T)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>%  
  dplyr::select(-geometry)

# Define the number of cores you want to use
num_cores <- 3  # Adjust this based on your system

# Create a cluster
cl <- makeCluster(num_cores)


# Register the cluster with foreach
registerDoSNOW(cl)



solutions <- foreach(i = 100:900, 
                     .packages = c("tidyverse", "DEoptim", "GWmodel", "sf"), 
                     .combine = "rbind",  
                     .verbose = T) %dopar% 
  {
                       
    de_result <- DEoptim(objective, lower = i, upper = i,
                         DEoptim.control(VTR = -Inf, strategy = 2,
                                         bs = FALSE, itermax = 10, CR = 0.5, F = 0.8, 
                                         storepopfreq = 1, trace = T, parallelType = "foreach"))
    output(de_result$optim$bestmem, de_result$optim$bestval*-1)
    }


write.table(solutions, "data/optimised/DEoptim_results.txt")
stopCluster(cl)



