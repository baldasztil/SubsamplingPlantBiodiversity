library(data.table)
library(tidyverse)# fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
library(cowplot)
library(beepr)
library(caper)
subsampling.plants <- function(spec_n) {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
  # accumulative subsampling 
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
  
  dist <- plantlist_dist %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # phy_div patterns across tdwg3 areas
  area_codes <- unique(dist$area)
  sample_rich_bru <- data.frame(matrix(nrow = length(area_codes), ncol = 2))
  names(sample_rich_bru) <- c("area_code_l3", "phy_div_sample")
  
  
  # could also have phylogenies per country
  drop_tree_sample <-drop.tip(output_tree,no_tip_needed$label)
  
  for (i in 1: length(area_codes)) {
    temp1 <- dist %>% 
      filter(area == area_codes[i])
    
    bru_tips <- temp1 %>% 
      filter(label %in% drop_tree_sample$tip.label)
    
    sample_rich_bru$area_code_l3[i] <- temp1$area_code_l3[1]
    sample_rich_bru$phy_div_sample[i] <- pd.calc(megamatrix, tip.subset = bru_tips, method = "TBL", root.edge=FALSE)
    
    #print(paste0("Currently located in ",area_codes[i], " (",i,"/",length(area_codes),")", "|", time.laps))
  }
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = c("area_code_l3")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  # measuring correlation
  
  #correlation to overall phy_div 
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor.test(phy_div_sample, phy_div,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(phy_div_sample, phy_div,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  
  # per continent correlation
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$phy_div_sample, .$phy_div, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$phy_div_sample, .$phy_div, 
                      method="pearson", data = ., exact =F)[[4]])  
  
  # combining the results
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = length(cumulative_namelist)) %>% 
    left_join(spear, by = "id")
  
  if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),1000)){
    print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  }
  
  #extracting cumulative pattern
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
}

# Import data ------------------------------------------------------------------
output_tree <- read.tree("data/gbmb_matched_wcvp2022.tre")
output_tree <- read.tree("data/phylomaker/output_tree.tre")
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")




#megamatrix <- clade.matrix(output_tree)


plants_full_accepted <- fread("data/wcvp_accepted_merged.txt") %>% 
  mutate(label = paste0(genus,"_",species))

search <- plants_full_accepted %>%  
  filter (plant_name_id %in% output_tree$tip.label) %>% 
  reframe(genus_n = n_distinct(genus))


plant_labels <- dplyr::select(plants_full_accepted, plant_name_id, label)

dist_native <- fread("data/dist_native.txt") %>% 
  left_join(plant_labels, by = "plant_name_id") %>% 
  filter(label %in% output_tree$tip.label)
# Import data ------------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3, label, area)

plantlist_names <- plants_full_accepted %>% 
  dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name, label)


dist.mat <- long2sparse(dist_native, grids = "area_code_l3", species = "label")
aa <- as.data.frame(PD(dist.mat, output_tree))
bb <- phylo_endemism(dist.mat, output_tree)
cc <- cbind(aa,bb)
cc$area_code_l3 <- rownames(aa)

