# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------
# This is a script to create random subsamples from the WCVP. It is using the
# taxonomic and geographic information available through the data 
# to relate species richness in the subsample to global patterns.  

# Libraries --------------------------------------------------------------------
library(tidyr)
library(dplyr)
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

# Defining functions and objects -----------------------------------------------
theme_default.lb <- function () {
  theme(
  # add border 1)
  panel.border = element_rect(colour = "black", fill = NA, linetype = 1, linewidth = 0.6),
  # color background 2)
  panel.background = element_rect(fill = NA),
  # modify grid 3)
  panel.grid.major.x = element_line(colour = "black", linetype = 3, linewidth = 0.5),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, linewidth = 0.5),
  panel.grid.minor.y = element_blank(),
  # modify text, axis and colour 4) and 5)
  axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 10),
  axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 10),
  axis.ticks = element_line(colour = "black"),
  #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
  # legend at the bottom 6)
  legend.position = "right", 
  legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 12), 
  legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 12), 
  legend.key.size = unit(1.5, "lines"),
  strip.background = element_rect(colour = "black", fill = alpha("#8685C4",0.5), linetype = 1, linewidth = 0.6),
  strip.text = element_text(colour = "black",face = "bold", family = "sans", size = 10),
  
  plot.title = element_text(hjust = 0.5, colour = "black", face = "plain", family = "sans", size = 14),
  panel.spacing.y = unit(1, "lines"),
  panel.ontop = F,
  plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
}
theme_default.lb.lines <- function () {
  theme(
    # add border 1)
    panel.border = element_rect(colour = "black", fill = NA, linetype = 1, linewidth = 0.6),
    # color background 2)
    panel.background = element_rect(fill = NA),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "black", linetype = 3, linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 12),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 12),
    axis.ticks = element_line(colour = "black"),
    #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "right", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 14), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
    legend.key.size = unit(1.5, "lines"),
    strip.background = element_rect(colour = "black", fill = alpha("#C2C6FFF8",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "bold", family = "sans", size = 16),
    
    plot.title = element_text(hjust = 0.5, colour = "black", face = "plain", family = "sans", size = 14),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F,
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
}
theme_default.lb.lines.big <- function () {
  theme(
    # add border 1)
    panel.border = element_rect(colour = "black", fill = NA, linetype = 1, linewidth = 0.6),
    # color background 2)
    panel.background = element_rect(fill = NA),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "black", linetype = 3, linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 14),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 16),
    axis.ticks = element_line(colour = "black"),
    #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "right", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 16), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 16), 
    legend.key.size = unit(1.5, "lines"),
    strip.background = element_rect(colour = "black", fill = alpha("#C2C6FFF8",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "bold", family = "sans", size = 16),
    
    plot.title = element_text(hjust = 0.5, colour = "black", face = "plain", family = "sans", size = 14),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F,
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
}

plotting.lines <- function (x) {

  
  cont_top5_filt <- cont_top5_mean %>% 
    filter(Continent == x) 
    
  cont_families_top5 <- threshold_cont %>%
    filter(Continent == x) %>%
    filter(family %in% cont_top5_filt$family) %>% 
    ungroup() %>%
    slice_max(ratio, by = c(Continent,family)) %>%
    arrange(desc(ratio)) %>%
    slice(1:5)
  
  plot <- samples_cumulative_rel %>% 
    filter(family %in% cont_families_top5$family) %>%
    ggplot(aes(x= log10(sp), y=mean, col = family)) +
    geom_line(lwd = 1)  +
    xlab("Sample size species (log10)") +
    ylab("Correlation coefficient") +
    facet_wrap(~Continent, nrow = 2, ncol = 5) +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
                linetype = 0, alpha = 0.1, ) +
    ylim(0,1) +
    scale_colour_manual("Family", values = pallete) +
    theme_default.lb.lines.big()
  plot
  ggsave(paste0("family_curves_top_5", x, 
                "_increasingsteps.png"),  
         width = 27, height = 15, dpi = 600)
  plot
}
plotting.lines.percont <- function (x) {

  cont_top5_filt <- cont_top5_mean %>% 
    filter(Continent == x) 
  
  cont_families_top5 <- threshold_cont %>%
    filter(Continent == x) %>%
    filter(family %in% cont_top5_filt$family) %>% 
    ungroup() %>%
    slice_max(ratio, by = c(Continent,family)) %>%
    arrange(desc(ratio)) %>%
    slice(1:5)
  
  plot <- samples_cumulative_rel %>% 
    filter(family %in% cont_families_top5$family & Continent == x[1]) %>%
    ggplot(aes(x= log10(sp), y=mean, col = family)) +
    geom_line(lwd = 1)  +
    xlab("Sample size species (log10)") +
    ylab("Correlation coefficient") +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
                linetype = 0, alpha = 0.1, ) +
    ylim(0,1) +
    facet_wrap(~Continent) +
    scale_colour_manual("Family", values = pallete) +
    theme_default.lb.lines()
ggsave(paste0("family_curves_top_5", x, 
              "_increasingsteps_solo.png"),  
       width = 10, height = 8, dpi = 600)
plot
}

plotting.maps <- function (x) {
  
  mapping <- rich_rel_out %>%
    filter(family == x) %>% 
    select(family, LEVEL3_COD = LEVEL_COD, richness_sample, n_sample)

  tdwg_left <- tdwg_3 %>% 
    left_join(mapping, "LEVEL3_COD") %>% 
    left_join(tdwg_1_names, by = c("LEVEL1_COD"), multiple = "all")
  
  assign("tdwg_left", tdwg_left ,envir = .GlobalEnv)
  plot <-  ggplot(data = tdwg_left)  +
    geom_sf(aes(fill = richness_sample)) +
    scale_fill_gradient("Species (n)", low = "grey90", high = "grey1", na.value = "white") +
    ggtitle(paste0(x, " (n = ", unique(mapping$n_sample), ")")) +  
    theme_default.lb()
  plot
  
 # ggsave(paste0("family_map_", x,".jpeg"),  
  #       width = 25, height = 15, dpi = 600)
}
plotting.maps.perc <- function (x) {
  
  mapping <- rich_rel_out %>%
    filter(family == x) %>% 
    select(family, LEVEL3_COD = LEVEL_COD, richness_sample, n_sample)
  
  tdwg_left <- tdwg_3 %>% 
    left_join(mapping, "LEVEL3_COD") %>% 
    left_join(tdwg_1_names, by = c("LEVEL1_COD"), multiple = "all") %>% 
    filter(LEVEL1_NAM == continent)
  
  assign("tdwg_left", tdwg_left ,envir = .GlobalEnv)
  plot <-  ggplot(data = tdwg_left)  +
    geom_sf(aes(fill = richness_sample)) +
    scale_fill_gradient("Species (n)", low = "grey90", high = "grey1", na.value = "white") +
    ggtitle(paste0(x, " (n = ", unique(mapping$n_sample), ")")) +  
    theme_default.lb()
  plot
  
  # ggsave(paste0("family_map_", x,".jpeg"),  
  #       width = 25, height = 15, dpi = 600)
}

plotting.families <- function (continent) {
  cumulative_sample <- list()
  assign("continent", continent, envir = .GlobalEnv)
  
  cont_top5_filt <- cont_top5_mean %>% 
    filter(Continent == continent) 
  
  families <- threshold_cont %>%
    filter(Continent == continent) %>%
    filter(family %in% cont_top5_filt$family) %>% 
    ungroup() %>%
    slice_max(ratio, by = c(Continent,family)) %>%
    arrange(desc(ratio)) %>%
    slice(1:5) 
  
  for (i in 1: nrow(families)) {
    plantlist_names <- plants_full %>% 
      filter(family %in% families$family[i])
    samplesize <- families %>% 
      filter(family %in% families$family[i])
    sample <- plantlist_names %>% 
      slice_sample(n = samplesize$sp) %>% 
      mutate(sp = samplesize$sp)
    cumulative_sample[[i]] <- sample   
    
  }
  
  cumulative_sample <- as.data.frame(do.call(bind_rows, cumulative_sample))
  index <- cumulative_sample %>%  
    select(family,plant_name_id,sp)
  
  dist <- dist_native %>% 
    filter(plant_name_id %in% cumulative_sample$plant_name_id) %>% 
    left_join(index, by = "plant_name_id")
  
  sample_rich_bru <- dist %>% 
    group_by(family, area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id),
              n_sample = unique(sp)) %>% 
    rename(LEVEL_COD = area_code_l3) 
  #pivot_wider(values_from = richness_sample names_from = family)
  
  rich_rel <- richness_patterns %>% 
    filter(ID == "bru") %>% 
    full_join(sample_rich_bru, by = "LEVEL_COD", multiple = "all") %>% 
    replace(is.na(.), 0) 
  
  assign("rich_rel_out", rich_rel, envir = .GlobalEnv)
  
  tdwg_3_mod <- tdwg_3 %>% 
    full_join(rich_rel, by = c("LEVEL3_COD" = "LEVEL_COD"), multiple = "all") %>% 
    full_join(tdwg_1_names, by = "LEVEL1_COD")
  
  assign("x", unique(index$family),  envir = .GlobalEnv)
  plotlist <- lapply(x, plotting.maps)
  
  pattern <- tdwg_3_mod %>%
    ggplot() +
    geom_sf(aes(fill = richness)) +
    scale_fill_gradient("Species (n)", low = "grey90", high = "grey1", na.value = "white") +
    ggtitle(paste0("Species richness observed")) +  
    theme_default.lb()
  
  plotlist[[6]] <- pattern
  
  assign(paste0("plotlist_",continent), plotlist, envir = .GlobalEnv)
  ggarrange(plotlist = plotlist, ncol =2, nrow =3)
  ggsave(paste0("family_maps_top5_global_",continent,".jpeg"),  
         width = 15, height = 13, dpi = 600)
  
}
plotting.families.perc <- function (continent) {
  cumulative_sample <- list()
  assign("continent", continent, envir = .GlobalEnv)
  
  cont_top5_filt <- cont_top5_mean %>% 
    filter(Continent == continent) 
  
  families <- threshold_cont %>%
    filter(Continent == continent) %>%
    filter(family %in% cont_top5_filt$family) %>% 
    ungroup() %>%
    slice_max(ratio, by = c(Continent,family)) %>%
    arrange(desc(ratio)) %>%
    slice(1:5)
  
  for (i in 1: nrow(families)) {
    plantlist_names <- plants_full %>% 
      filter(family %in% families$family[i])
    samplesize <- families %>% 
      filter(family %in% families$family[i])
    sample <- plantlist_names %>% 
      slice_sample(n = samplesize$sp) %>% 
      mutate(sp = samplesize$sp)
    cumulative_sample[[i]] <- sample   
    
  }
  
  cumulative_sample <- as.data.frame(do.call(bind_rows, cumulative_sample))
  index <- cumulative_sample %>%  
    select(family,plant_name_id,sp)
  
  dist <- dist_native %>% 
    filter(plant_name_id %in% cumulative_sample$plant_name_id) %>% 
    left_join(index, by = "plant_name_id")
  
  sample_rich_bru <- dist %>% 
    group_by(family, area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id),
              n_sample = unique(sp)) %>% 
    rename(LEVEL_COD = area_code_l3) 
  #pivot_wider(values_from = richness_sample names_from = family)
  
  rich_rel <- richness_patterns %>% 
    filter(ID == "bru") %>% 
    full_join(sample_rich_bru, by = "LEVEL_COD", multiple = "all") %>% 
    replace(is.na(.), 0) 
  
  assign("rich_rel_out", rich_rel, envir = .GlobalEnv)
  
  
  tdwg_3_mod <- tdwg_3 %>% 
    full_join(rich_rel, by = c("LEVEL3_COD" = "LEVEL_COD"), multiple = "all") %>% 
    full_join(tdwg_1_names, by = "LEVEL1_COD")
  
  assign("x", unique(index$family),  envir = .GlobalEnv)
  plotlist <- lapply(x, plotting.maps.perc)
  
  pattern <- tdwg_3_mod %>%
    filter(LEVEL1_NAM == continent) %>% 
    ggplot() +
    geom_sf(aes(fill = richness)) +
    scale_fill_gradient("Species (n)", low = "grey90", high = "grey1", na.value = "white") +
    ggtitle(paste0("Species richness observed")) +  
    theme_default.lb()
  
  plotlist[[6]] <- pattern
  
  assign(paste0("plotlist_",continent), plotlist, envir = .GlobalEnv)
  ggarrange(plotlist = plotlist, ncol =2, nrow =3)
  ggsave(paste0("family_maps_top5_",continent,".jpeg"),  
         width = 10, height = 10, dpi = 600)
  
}

# Import data ------------------------------------------------------------------
richness_patterns <- read.csv("output/Overall_richness.csv")
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/output/families/", full.names = TRUE)
#Read them in a list
samplelist <- lapply(list_data, read.table)

#make big dataframe
data_out_raw <- as.data.frame(do.call(rbind, samplelist))

data_out <- data_out_raw %>% 
  filter(!family == "Acanthaceae" | !sp== 5422)

# Reframe ----------------------------------------------------------------------
samples_cumulative_rel <- data_out %>%
  group_by(id,sp,family) %>% 
  mutate( cor.sp = ifelse(is.na(cor.sp), 0, cor.sp)) %>% 
  summarise(mean = mean(cor.sp, na.rm = TRUE),
            sd = sd(cor.sp, na.rm = TRUE),
            n.sample = n()) %>%
  mutate(se = sd / sqrt(n.sample),
         lower.ci = mean - qt(1 - (0.05 / 2), n.sample - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.sample - 1) * se) %>% 
  left_join(tdwg_1_names, by = c("id" = "LEVEL1_COD")) %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "OVERALL", Continent))
  
#Output analysis----------------------------------------------------------------

threshold <- samples_cumulative_rel %>% 
  mutate(threshold = case_when(mean > 0.8  ~ 'above', 
                               mean < 0.8 ~ 'below')) %>% 
  group_by(family) %>% 
  filter(min(sp) & threshold == "above") %>%
  slice_min(sp) %>% 
  arrange(sp) %>% 
  head(n= 5)

threshold_cont <- samples_cumulative_rel %>%
  mutate(threshold = case_when(mean > 0.8 & sp > 100 ~ 'above', 
                               mean < 0.8 | sp < 100 ~ 'below')) %>% 
  group_by(family) %>% 
  filter(min(sp) & threshold == "above") %>% 
  mutate(ratio = mean/sp)

  
overall_top5_mean <- threshold_cont %>%
  ungroup() %>%
  filter(Continent == "OVERALL") %>% 
  slice_max(mean, by = c(Continent,family)) %>%
  arrange(desc(mean)) %>%
  slice(1:5) 

cont_top5_mean <- threshold_cont %>%
  ungroup() %>% 
  slice_max(mean, by = c(Continent,family)) %>%
  group_by(Continent) %>% 
  arrange(desc(mean)) %>%
  slice(1:5) 

cont_top5_ratio <- threshold_cont %>%
  ungroup() %>% 
  slice_max(ratio, by = c(Continent,family)) %>%
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:5) 

overall_top5_ratio <- threshold_cont %>%
  ungroup() %>%
  filter(Continent == "OVERALL") %>% 
  slice_max(ratio, by = c(Continent,family)) %>%
  arrange(desc(ratio)) %>%
  slice(1:5)

cont_top5 <- threshold_cont %>%
  filter(family %in% cont_top5_mean$family) %>% 
  ungroup() %>% 
  slice_max(ratio, by = c(Continent,family)) %>%
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:5) 

overall_top5 <- threshold_cont %>%
  filter(family %in% cont_top5_mean$family) %>% 
  ungroup() %>%
  filter(Continent == "OVERALL") %>% 
  slice_max(ratio, by = c(Continent,family)) %>%
  arrange(desc(ratio)) %>%
  slice(1:5)
  
tdwg_1_names$LEVEL1_COD <-  as.numeric(tdwg_1_names$LEVEL1_COD)


samples_cumulative_rel$Continent <- as.factor(samples_cumulative_rel$Continent)
samples_cumulative_rel$Continent<-factor(samples_cumulative_rel$Continent, levels=c("AFRICA","ANTARCTIC","ASIA-TEMPERATE","ASIA-TROPICAL",
                                                                                    "AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC","SOUTHERN AMERICA","OVERALL"))



#Plotting-----------------------------------------------------------------------
#pallete <- wes_palette(name  = "BottleRocket2", 5, type = c("discrete"))

#correlation curves
pallete <- as.vector(GetColors(5, scheme = "mute"))
continent_names <- unique(threshold_cont$Continent)

#correlation plots for each continent global
lapply(continent_names, plotting.lines)

#by continent single  
continent_corr_top5 <- lapply(continent_names, plotting.lines.percont)

# combining the plots of each continent
ggarrange(plotlist = continent_corr_top5, ncol = 5, nrow = 2, labels = "AUTO")
ggsave(paste0("family_curves_top_5", 
              "_increasingsteps_solo_combined.png"),  
       width = 35, height = 10, dpi = 600)



# richness map

#richness per continent global
lapply(continent_names,plotting.families)

#by continent, mapping continent single 
continent_names <- continent_names[!(continent_names == "OVERALL")]
lapply(continent_names,plotting.families.perc)




# richness overall 
tdwg_3_mod <- tdwg_3 %>% 
  left_join(richness_patterns, by = c("LEVEL3_COD" = "LEVEL_COD"), multiple = "all") 

pattern <- tdwg_3_mod %>%
  ggplot() +
  geom_sf(aes(fill = richness)) +
  scale_fill_gradient("Species (n)", low = "grey90", high = "grey1", na.value = "white") +
  ggtitle(paste0("Species richness observed")) +  
  theme_pubclean() +
  theme(panel.background = element_rect(color = NA))

pattern
