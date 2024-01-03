# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------
# This is a script to create randomom subsamples from the WCVP. It is using the
# taxonomic and geographic information available through the data 
# to relate species richness in the subsample to global patterns.  

# Libraries --------------------------------------------------------------------
library(tidyr)
library(purrr)
library(stringr)
library(data.table)
library(sf)
library(ggplot2)
library(ggpubr)
library(ggokabeito)
library(inlmisc)
library(wesanderson) #darjeeling 1, rushmore 1
library(ggsci)
library(tmap)
library(dplyr)
library(gridExtra)
library(forcats)
library(viridis)
library(tmaptools)
library(RColorBrewer)
tmap_options(check.and.fix = TRUE)


theme_default.lb.boxplot <- function () {
  theme(
    # add border 1)
    panel.border = element_rect(colour = "black", fill = NA),
    # color background 2)
    panel.background = element_rect(fill = NA),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "black", linetype = 3, linewidth =  0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, linewidth =  0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 26),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 26),
    axis.ticks = element_line(colour = "black"),
    #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "none", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 21), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 26), 
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.8),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 26),
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
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, linewidth =  0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 26),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 26),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black", linetype = 1, linewidth =  0.5),
    # legend at the bottom 6)
    legend.position = "bottom", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 20), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 18), 
    
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.8),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 21),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F,
    legend.key.size = unit(1, "lines")
  )
}
objective <- function(subset_indices, plantlist_names) {
  list <- list()
  for (i in 1:1000) {
    species <- slice_sample(n = round(subset_indices),plantlist_names)
    dist <- plantlist_dist %>%
      filter(plant_name_id %in% species$plant_name_id)
    sample_rich_bru <- dist %>%
      group_by(area_code_l3) %>%
      summarise(richness_sample = n_distinct(plant_name_id)) %>%
      rename(LEVEL3_COD = area_code_l3)
    
    list[[i]] <- rich_overall_bru %>%
      left_join(sample_rich_bru, by = "LEVEL3_COD") %>%
      replace(is.na(.), 0) %>%
      mutate(sp = length(species))
  }
  data_out <- as.data.frame(do.call(rbind, list))
  output_file <- data_out %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(mean_richness_samp = mean(richness_sample),
            sd_richness_samp = sd(richness_sample), 
            LEVEL1_COD = unique(LEVEL1_COD), 
            LEVEL1_NAM = unique(LEVEL1_NAM),
            sp = subset_indices)
            
#  aa <- aggregate(data=data_out,richness_sample~.,FUN = mean)
  
  return(output_file)  # Return negative correlation (since DEoptim minimizes)
  #return(ifelse(correlation >= 0.95, correlation, 0))
}

tdwg_3_mod <- st_read(dsn = "data/wgsrpd-master/level_3_richness")


tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
rich_overall_bru <- tdwg_3_mod %>% 
  dplyr::select(richness, LEVEL3_COD, LEVEL3_NAM, LEVEL1_COD, LEVEL1_NAM) %>% 
  st_drop_geometry()

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  dplyr::select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

plants_full_all <- fread("data/wcvp_accepted_merged.txt")
srli <- read.csv("data/red/srli/srli_full.csv")
redlist <- read.csv("data/red/redlist_full_syn.csv")

dist_native_all <- fread("data/dist_native.txt")

plantlist_names <- plants_full_all %>%  
  dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native_all %>% 
  filter(plant_name_id %in% plants_full_all$plant_name_id)

samples_cumulative_rel <- fread("output/cumulative_rel/fullsamples_cumulative_rel.txt") %>% 
  mutate(method = "random",
         dataset = "wcvp") %>% 
  dplyr::select(id, method, sp, cor = mean, dataset)

samples_cumulative_rel_red <- fread("output/cumulative_rel/redlist_cumulative_rel.txt") %>%  
  mutate(method = "random", 
         dataset = "redlist") %>% 
  dplyr::select(id, method, sp, cor = mean, dataset)

samples_cumulative_rel_srli <- fread("output/cumulative_rel/srli_cumulative_rel.txt") %>% 
  mutate(method = "random", 
         dataset = "srli") %>% 
  dplyr::select(id, method, sp, cor = mean, dataset)

samples_rand <- rbind(samples_cumulative_rel, samples_cumulative_rel_red, samples_cumulative_rel_srli)


max_vals <- samples_rand %>%  
  group_by(dataset,id) %>% 
  reframe(max.cor = max(cor))

max_vals$cor.missing <- max_vals$max.cor - 1 

red_max_map <- max_vals %>% 
  filter(dataset %in% c("redlist", "srli") & !id == "overall") %>% 
  mutate(LEVEL1_COD = as.numeric(id))


mapping_red <- tdwg_3_mod %>%  
  left_join(red_max_map, by = "LEVEL1_COD", multiple = "all")
unique(st_is_valid(mapping_red))


tm_shape(mapping_red) + 
  tm_fill(col = "cor.missing", n=5, palette="Reds", legend.is.portrait = FALSE) +
  tm_legend(outside=T,legend.outside.position = "bottom", legend.outside.size = 0.1)  + 
  tm_facets(by = "dataset", nrow = 1, ncol = 2)  +
  tm_layout(frame = FALSE, asp = 2, panel.show = F, legend.position=c(0.25,0), 
            panel.labels = c("Red List", "SRLI"), 
            panel.label.bg.color = "black") 
#tm_graticules(col = "grey", alpha = 0.7, labels.show = F) +
#tm_compass(type = "4star", size = 1, position = c("left", "top")) +
#tm_scale_bar(breaks = c(0, 1500, 3000), text.size = 4, position = c("right", "top")) 

mapping_redlist <- mapping_red %>%  filter(dataset == "redlist") %>% 
  rename(Bias= cor.missing) %>% 
  tm_shape() + 
  tm_fill(col = "Bias", n=5, palette="Reds")  + tm_format("World") +
  tm_layout(frame = FALSE, legend.outside = F, legend.position = c(0.04,0.1)) 

mapping_srli <- mapping_red %>%  filter(dataset == "srli") %>% 
  rename(Bias= cor.missing) %>% 
  tm_shape() + 
  tm_fill(col = "Bias", n=5, palette="Reds")  + tm_format("World") +
  tm_layout(frame = FALSE, legend.outside = F, legend.position = c(0.04,0.1), 
            main.title = "SRLI", main.title.position = "center") 


#tmap_save(mapping_redlist, "Redlist_Bias.png",  dpi = 600)
tmap_arrange(mapping_redlist, mapping_srli, nrow = 2, ncol = 1, asp = NA, outer.margins = 0)


wcvp_means_red <- objective(nrow(redlist), plants_full_all)
wcvp_means_red$dataset <- "wcvp_red"

wcvp_means_srli <- objective(nrow(srli), plants_full_all)
wcvp_means_srli$dataset <- "wcvp_srli"

redlist_means <- objective(nrow(redlist), redlist)
redlist_means$dataset <- "redlist"

srli_means <- objective(nrow(srli), srli)
srli_means$dataset <- "srli"

means_sd <- rbind(wcvp_means_red,wcvp_means_srli, redlist_means, srli_means)


a <- wcvp_means_red %>% 
  dplyr::select(LEVEL3_COD, wcvp_rich_red = mean_richness_samp, wcvp_sd = sd_richness_samp) 

b <- wcvp_means_srli %>% 
  dplyr::select(LEVEL3_COD, wcvp_rich_srli = mean_richness_samp, wcvp_sd = sd_richness_samp) 

c <- redlist_means %>% 
  dplyr::select(LEVEL3_COD, redlist_rich = mean_richness_samp, redlist_sd = sd_richness_samp) 

d <- srli_means %>% 
  dplyr::select(LEVEL3_COD, srli_rich = mean_richness_samp, srli_sd = sd_richness_samp) 

tdwg_3_means <- tdwg_3_mod %>% 
  left_join(a, by = "LEVEL3_COD") %>% 
  left_join(b, by = "LEVEL3_COD") %>% 
  left_join(c, by = "LEVEL3_COD") %>% 
  left_join(d, by = "LEVEL3_COD")

tdwg_3_means$Mean_Difference <- tdwg_3_means$redlist_rich / tdwg_3_means$wcvp_rich_red - 1
tdwg_3_means$Difference_SRLI <- tdwg_3_means$srli_rich / tdwg_3_means$wcvp_rich_srli - 1

#RColorBrewer::brewer.pal.info
RColorBrewer::brewer.pal(7, "PuOr")
pal <- c("#B35806", "#F1A340", "#FEE0B6", "#F7F7F7", "#D8DAEB", "#998EC3", "#542788")

#pal <- c("#E66101", "#FDB863", "#F7F7F7", "#B2ABD2", "#5E3C99")
#pal <- c("#D6604D" ,"#F4A582" , "#fef7cd" ,"#92C5DE" ,"#4393C3")
mapping_redlist <- 
  tm_shape(tdwg_3_means) + 
  tm_fill(col = "Mean_Difference", n=7, palette= pal,breaks = c(-1,-0.75,-0.5,-0.25,0.25, 0.5,0.75, 1), midpoint = 0, legend.is.portrait = F)  +
  tm_legend(outside=T,legend.outside.position = "bottom", legend.outside.size = 0.1, scale = 10) + tm_format("World") +
  tm_layout(frame = FALSE, legend.outside = T, legend.position = c(0.35,0)) 

mapping_redlist
tmap_save(mapping_redlist, "Redlist_Bias_meanspecies_full.png",  dpi = 1200)


mapping_srli <- 
  tm_shape(tdwg_3_means) + 
  tm_fill(col = "Difference_SRLI", n=5, palette= pal, breaks = c(-1,-0.75,-0.5,-0.25,0.25, 0.5,0.75, 1), midpoint = 0, legend.show = F)  + tm_format("World") +
  tm_layout(frame = FALSE, legend.outside = F, legend.position = c(0.04,0.1)) 

mapping_srli
tmap_save(mapping_srli, "SRLI_Bias_meanspecies_full.png",  dpi = 1200)

tmap_arrange(mapping_redlist, mapping_srli, nrow = 2, ncol = 1, asp = NA, outer.margins = 0)


