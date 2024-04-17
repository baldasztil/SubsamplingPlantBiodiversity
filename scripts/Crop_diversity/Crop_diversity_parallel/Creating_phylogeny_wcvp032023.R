# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------
# This is a script to create random subsamples from the WCVP. It is using the
# taxonomic and geographic information available through the data 
# to relate species richness in the subsample to global patterns.  

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(ape)
library(U.PhyloMaker)
# Defining functions -----------------------------------------------------------
theme_default.lb.boxplot <- function () {
  theme(
    # add border 1)
    #panel.border = element_rect(colour = "black", fill = NA),
    # color background 2)
    panel.background = element_rect(fill = NA),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "black", linetype = 3, size = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, size = 0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 9),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "right", 
    legend.title = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 12), 
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 14),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F)
}
theme_default.lb <- function () {
  theme(
    # add border 1)
    #panel.border = element_rect(colour = "black", fill = NA),
    # color background 2)
    panel.background = element_rect(fill = NA),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "black", linetype = 3, size = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, size = 0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 9),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "right", 
    legend.title = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 12), 
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 14),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F)
}
getwd()
# Import data ------------------------------------------------------------------
megatree <- read.tree("data/phylomaker/plant_megatree.tre")
output_tree <- read.tree("data/gbmb_matched_wcvp2022.tre")

gen_list <- read.csv("data/phylomaker/plant_genus_list.csv")

plants_full <- fread("data/wcvp/wcvp_names_032023.csv")
plants_full$taxon_name <- str_replace(plants_full$taxon_name, " ", "_")




tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
richness_patterns <- fread("data/richness_patterns.txt")
dist_native <- fread("data/dist_native.txt")
plants_full_accepted <- fread("data/wcvp_accepted_merged.txt")
plants_full_accepted$taxon_name <- str_replace(plants_full_accepted$taxon_name, " ", "_")

search <- plants_full_accepted %>%  
  filter (taxon_name %in% megatree$tip.label) %>% 
  reframe(genus_n = n_distinct(genus))

gen_number <- length(unique(plants_full_accepted$genus))

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  dplyr::select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


#plants_formated_sample <- sample_n(plants_full_accepted, 5000) #nrow(plants_full_accepted)

plants_formated_sample <- plants_full_accepted


#sp_list <- plants_formated_sample %>% 
 # dplyr::select(species = taxon_name) %>% 
  #mutate(family = NA,
   #      genus = NA,
    #     species.relative = NA, 
     #    genus.relative = NA)  
  
sp_list <- plants_formated_sample %>% 
  dplyr::select(species = taxon_name, genus, family) %>% 
  mutate(
         species.relative = NA, 
         genus.relative = NA)  


gen_list_wcvp <- plants_formated_sample %>% 
  filter(taxon_rank == "Genus")


gen_list <- plants_full %>% 
  filter(taxon_rank == "Genus") %>% 
  group_by(genus, family, accepted_plant_name_id) %>% 
  summarise()


tip_labels <- as.data.frame(megatree$tip.label) 

tip_labels <- tip_labels %>% 
  rename(taxon_name = `megatree$tip.label`) %>% 
  mutate(taxon_name2 = taxon_name, 
         test_original = taxon_name == taxon_name2)


name_list_tips <- plants_full %>% 
  filter(taxon_name %in% tip_labels$taxon_name) %>%  
  filter(!taxon_status == "Accepted") %>% 
  dplyr::select(taxon_name, accepted_plant_name_id) %>% 
  filter(!is.na(accepted_plant_name_id))


accepted_names_list <- plants_full %>% 
  filter(accepted_plant_name_id %in% name_list_tips$accepted_plant_name_id) %>% 
  filter(taxon_status == "Accepted") %>% 
  mutate(accepted_name = taxon_name) %>% 
  dplyr::select(accepted_name, accepted_plant_name_id)


names_reformed <- name_list_tips %>% 
  left_join(accepted_names_list, by = "accepted_plant_name_id") %>% 
  mutate(accepted_name = word(accepted_name, 1)) %>% 
  dplyr::select(accepted_name, taxon_name)

names_good <- plants_full_accepted %>% 
  filter(taxon_name %in% tip_labels$taxon_name) %>%  
  filter(taxon_status == "Accepted")  %>% 
  mutate(accepted_name = taxon_name) %>% 
  dplyr::select(accepted_name, taxon_name) %>% 
  rbind(names_reformed) %>% 
  filter(!duplicated(.)) %>% 
  mutate(test = taxon_name == accepted_name)


names_comb <- names_good %>% 
  right_join(tip_labels, by = "taxon_name") %>% 
  filter(!duplicated(taxon_name)) %>% 
  mutate(accepted_name = ifelse(is.na(accepted_name), taxon_name, accepted_name)) %>% 
  dplyr::select(accepted_name, taxon_name, test)

tip_labels_new <- tip_labels %>% 
  left_join(names_comb, by = "taxon_name") %>% 
  dplyr::select(accepted_name)

megatree$tip.label <- tip_labels_new$accepted_name


genus_megatree <- unique(str_split(megatree$tip.label, "_", simplify = T)[,1])

gen_list_wcvp_accepted <- plants_formated_sample %>% 
  group_by(genus, family) %>% 
  summarise()

gen_list_match <- gen_list_wcvp_accepted %>% 
  filter(genus %in% genus_megatree)

gen_list_nomatch <- gen_list_wcvp_accepted %>% 
  filter(!genus %in% genus_megatree)


#sp_list2 <- plants_formated_sample %>% 
 # dplyr::select(species = taxon_name, genus, family) 

#xx <- build.nodes.2(megatree, sp_list2)
#table(xx$genus)

result <- phylo.maker(sp_list, megatree, gen_list_wcvp_accepted, nodes.type = 2, scenario = 3, output.tree	= T)
save.image("output_tree_wcvp")
output_tree <- result$phylo

check <- result$sp.list
unique(result$sp.list$output.note)


check2 <- check %>% 
  filter(output.note == "no insertion")
check2$species <- gsub(" ", "_", check2$species)

library(phytools)
xx <- add.random(result$phylo, tips = check2$species)
write.tree(result$phylo, "data/output_tree_20012024.tre")
write.table(result$sp.list, "data/output_splist_20012024.txt")
