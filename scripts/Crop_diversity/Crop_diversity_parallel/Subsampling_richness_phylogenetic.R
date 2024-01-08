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
gen_list <- read.csv("data/phylomaker/plant_genus_list.csv")

plants_full <- fread("data/wcvp/wcvp_names.txt")


tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
richness_patterns <- fread("data/richness_patterns.txt")
dist_native <- fread("data/dist_native.txt")
plants_full_accepted <- fread("data/wcvp_accepted_merged.txt")




tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


plants_formated_sample <- sample_n(plants_full_accepted, 500)



sp_list <- plants_full_accepted %>% 
  select(species = taxon_name) %>% 
  mutate(family = NA,
         genus = NA,
         species.relative = NA, 
         genus.relative = NA)  
  


gen_list_wcvp <- plants_full %>% 
  filter(taxon_rank == "Genus")

gen_list_wcvp_accepted <- plants_full %>% 
  filter(taxon_rank == "Genus") %>% 
  filter(taxon_status == "Accepted") %>% 
  filter(genus_hybrid == "")
  

genus_match <- gen_list %>% 
  left_join(gen_list_wcvp, by = "genus", multiple = "all") %>% 
  filter(taxon_status == "Accepted") %>% 
  filter(genus %in% plants_formated_sample$genus) %>% 
  select(genus, family = family.x)

genus_match <- gen_list %>% 
  left_join(gen_list_wcvp, by = "genus", multiple = "all") %>% 
  filter(taxon_status == "Accepted" ) %>% 
  filter(genus_hybrid == "") %>% 
  select(genus, family = family.x) 


genus_missing <- gen_list_wcvp_accepted %>%  
  filter(!genus %in% gen_list$genus)


genus_nomatch <- gen_list %>% 
  filter(!genus %in% genus_match$genus) %>% 
  select(genus)


result <- phylo.maker(sp_list, megatree, genus_match, nodes.type = 1, scenario = 3, output.tree	= T)


write.tree(result$phylo, "data/phylogeny/output_tree.tre")
write.table(result$sp.list, "data/phylogeny/output_splist.txt")
