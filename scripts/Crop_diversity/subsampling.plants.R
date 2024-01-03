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
    mutate(sp = length(cumulative_namelist))
  # measuring correlation
  
  # overall
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = as.numeric(cor(richness_sample, richness, 
                                      method="pearson"))) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
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
  
  #if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
  # print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  #}
  
  # cumulative pattern
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
}
parallel.subsampling <- function (iteration,sample_n) {
  samples <- list()
  for (i in 1:iteration)  {
    samples[[i]] <- parLapply(clust, seq(0,nrow(plantlist_names),sample_n), subsampling.plants)
  }
  return(samples)
}


aa <- parallel.subsampling(1,10000)
###### evaluation- -------------------------------------------------------------
clust <- makeCluster(7)
clusterExport(clust, c("dist_native","plantlist_dist","plantlist_names","rich_overall_bru"))
clusterEvalQ(clust, library(tidyverse))
doParallel::registerDoParallel(clust)

samples_rel <- parallel.subsampling(10,10000)

stopCluster(clust)