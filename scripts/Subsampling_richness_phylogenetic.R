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
library(phytools)
library(caper)
library(dplyr)
library(ggtree)

options(ragg.max_dim = 10000000000)
# Defining functions -----------------------------------------------------------
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

# Import data ------------------------------------------------------------------
megatree <- read.tree("data/phylomaker/plant_megatree.tre")
gen_list <- read.csv("data/phylomaker/plant_genus_list.csv")

output_tree <- read.tree("data/phylomaker/output_tree.tre")
out_splist <- fread("data/phylomaker/output_splist.txt")


# Import data ------------------------------------------------------------------

dtips<-sample(tree$tip.label,10)
dt<-drop.tip(tree,dtips)
plotTree(dt)

a <- ggtree(output_tree, layout="circular")

ggsave(file = paste0("fulltree.png"),  a, 
       width = 40, height = 40, dpi = 600, limitsize = F, )

missing_sp <- out_splist %>% 
  filter(output.note == "no insertion" ) 


missing_plants <- plants_full %>% 
filter(taxon_name %in% missing_sp$species)

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
richness_patterns <- fread("data/richness_patterns.txt")
dist_native <- fread("data/dist_native.txt")
plants_full_accepted <- fread("data/wcvp_accepted_merged.txt")



tdwg_3_names <- as.data.frame(tdwg_3) %>% 
  left_join(tdwg_1, by = "LEVEL1_COD") %>% 
  dplyr::select(LEVEL3_NAM, LEVEL3_COD, LEVEL1_COD, LEVEL1_NAM)

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  dplyr::select(LEVEL1_NAM, LEVEL1_COD)

species_per_cont <- dist_native %>% 
  group_by(continent_code_l1) %>% 
  reframe(plant_name_id = unique(plant_name_id)) %>% 
  left_join(tdwg_1_names, by = c("continent_code_l1" = "LEVEL1_COD"))
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")
  

species_filt <- species_per_cont %>% 
  filter(LEVEL1_NAM == "SOUTHERN AMERICA") %>% 
  left_join(plants_full_accepted, by = "plant_name_id")


#plants_formated_sample <- sample_n(plants_full_accepted, 500)




sp_list <- as.data.frame(species_filt) %>% 
  dplyr::select(species = taxon_name) %>% 
  mutate(genus = NA,
         family = NA,
         species.relative = NA, 
         genus.relative = NA)  
  


gen_list_wcvp <- plants_full_accepted %>% 
  group_by(family) %>% 
  reframe(genus = unique(genus))
  

  

genus_match <- as.data.frame(gen_list) %>% 
  left_join(gen_list_wcvp, by = "genus", multiple = "all") %>% 
  filter(genus %in% species_filt$genus) %>% 
  dplyr::select(genus, family = family.x)


genus_missing <- gen_list_wcvp %>%  
  filter(!genus %in% gen_list$genus)


genus_nomatch <- gen_list %>% 
  filter(!genus %in% gen_list_wcvp$genus) %>% 
  dplyr::select(genus)




phylo.maker(sp_list, megatree, genus_match,  scenario = 3, output.tree	= T, r = 1)

splist <- result$sp.list
output_tree <- result$phylo
#write.tree(result$phylo, "/data/phylogeny/output_tree_500sp.tre")
#write.table(result$sp.list, "/data/phylogeny/output_splist_500sp.txt")

xx <- clade.matrix(output_tree) 

a <- pd.calc(output_tree, tip.subset = NULL, method = "TBL", root.edge=FALSE)

ed_paci <- ed.calc(xx, polytomy.cf=c("isaac","mooers","none"))

ed_paci <- ed_paci$spp

clade_mat <- as.data.frame(xx$clade.matrix)
