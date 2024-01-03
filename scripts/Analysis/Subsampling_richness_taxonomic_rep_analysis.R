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
library(RColorBrewer)
library(ggokabeito)
library(khroma)
library(tmaptools)
library(reactablefmtr)
library(gridExtra)

# Defining functions -----------------------------------------------------------
theme_default.lb.boxplot <- function () {
  theme(
    # add border 1)
    #panel.border = element_rect(colour = "black", fill = NA),
    # color background 2)
    panel.background = element_rect(fill = NA),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "black", linetype = 3, linewidth =  0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, linewidth =  0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 28),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 28),
    axis.ticks = element_line(colour = "black"),
    #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "right", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 28), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 28), 
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 28),
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
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 28),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 28),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black", linetype = 1, linewidth =  0.5),
    # legend at the bottom 6)
    legend.position = "bottom", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 28), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 28), 
    
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 28),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F,
    legend.key.size = unit(1, "lines")
  )
}

plotting.lines <- function (x) {
  overall_top5 <- threshold_cont %>%
    ungroup() %>%
    filter(Continent == x) %>% 
    slice_max(ratio, by = c(Continent,family)) %>%
    arrange(desc(ratio)) %>%
    slice(1:5) 
  
  plot <- samples_cumulative_rel %>% 
    filter(family %in% overall_top5$family) %>%
    ggplot(aes(x= log10(sp), y=mean, col = family)) +
    geom_line(lwd = 1)  +
    xlab("Sample size species (log10)") +
    ylab("Correlation coefficient") +
    facet_wrap(~Continent, nrow = 2, ncol = 5) +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
                linetype = 0, alpha = 0.1, ) +
    ylim(0,1) +
    scale_color_nejm() +
    theme_default.lb()
  plot
  ggsave(paste0("family_curves_top_", x, 
                "_increasingsteps_logsp_newtheme.jpeg"),  
         width = 25, height = 15, dpi = 600)
  plot
}

options(digits = 2)  
# Import data ------------------------------------------------------------------

richness_patterns <- read.csv("output/Overall_richness.csv")
plants_full <- fread("data/wcvp_accepted_merged.txt") %>% 
  select(-continent)
dist_native <- fread("data/dist_native.txt")

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_3_names <- as.data.frame(tdwg_3) %>% 
  select(LEVEL3_COD, LEVEL2_COD, LEVEL1_COD) %>%  
  left_join(tdwg_1_names, by ="LEVEL1_COD")

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/taxonomic_rep", full.names = TRUE)
#Read them in a list
samplelist <- lapply(list_data, fread)

#make big dataframe
data_out <- as.data.frame(do.call(rbind, samplelist))

# Reframe ----------------------------------------------------------------------
samples_cumulative_rel <- data_out %>%
  group_by(id,sp) %>% 
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
  mutate(threshold = case_when(mean > 0.95  ~ 'above', 
                               mean < 0.95 ~ 'below')) %>% 
  group_by(Continent) %>% 
  filter(threshold == "above") %>%
  mutate(ratio = mean/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 



richness_patterns <- dist_native %>%
  left_join(plants_full, by = "plant_name_id") %>% 
  group_by(area_code_l3) %>% 
  summarise(rich_family = n_distinct(family)) %>% 
  left_join(tdwg_3_names, by =c("area_code_l3" = "LEVEL3_COD")) 





#Summary statistics analysis----------------------------------------------------

#Plotting-----------------------------------------------------------------------




#a <- colour("muted")
colouring <- c( '#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', "#CC3311", "#00183E")

viri <- viridisLite::viridis(100)


samples_cumulative_rel$Continent <- as.factor(samples_cumulative_rel$Continent)
samples_cumulative_rel$Continent<-factor(samples_cumulative_rel$Continent, levels=c("AFRICA","ANTARCTIC","ASIA-TEMPERATE","ASIA-TROPICAL",
                                            "AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC","SOUTHERN AMERICA","OVERALL"))



curves <- ggplot(data = samples_cumulative_rel, aes(x=log10(sp), y=mean, col = Continent)) +
  geom_line(lwd = 1, alpha = 0.5)  +
  geom_line(data = subset(samples_cumulative_rel, Continent == 'OVERALL'), 
            lwd = 2, color = "midnightblue", alpha = 0.7) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  scale_colour_manual(values = c(colouring)) +
  theme_default.lb() + 
  xlim(c(0,6)) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 4, nrow = 3))
  
ggsave(file = paste0("curves_taxrep", 
                     "_cropdiv_it100x.jpeg"),curves,   
       width = 20, height = 15, dpi = 600)

richness_map <- tdwg_3 %>% 
  left_join(richness_patterns, by =c("LEVEL3_COD" = "area_code_l3"))

pattern <- ggplot(data = richness_map) +
  geom_sf(aes(fill = rich_family)) +
  scale_fill_gradient("Family (n)", low = "grey90", high = "grey1", na.value = "white") +
  ggtitle("Family richness (observed)") +  
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
pattern

# richness maps + curves
g <- arrangeGrob(curves, pattern, ncol = 2)
ggsave(file = paste0("maps_curves", 
                     "_cropdiv_taxonomic_rep.jpeg"),g,   
       width = 25, height = 10, dpi = 600)


#ggarrange(curves, pattern, sample, ncol =2, nrow =2, align = c("hv"), heights = c(2, 1))


#grid.arrange(curves, arrangeGrob(pattern, sample), ncol=2,top="Correlation Subsample")






tdwg_2_mod <- tdwg_2 %>% 
  filter(LEVEL2_COD %in%richness_patterns_reg$LEVEL2_COD) %>% 
  left_join(richness_patterns_reg, by = "LEVEL2_COD") 
reg <- ggplot(data = tdwg_2_mod) +
  geom_sf(aes(fill = richness)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick2", na.value = "white")

ggarrange(bru, reg, con,ncol =3, nrow =1)




sample_ratios <- samples_cumulative_rel %>% 
  filter(sp >100) %>% 
  mutate(ratio = mean/sp) %>% 
  ggplot(aes(x=log10(sp), y=ratio, col = Continent)) +
  geom_point() +
  scale_colour_manual(values = c(colouring))
sample_ratios




