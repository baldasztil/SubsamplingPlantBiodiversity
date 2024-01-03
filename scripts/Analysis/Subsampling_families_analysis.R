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
# Defining functions -----------------------------------------------------------
theme_default.lb <- function () {
  theme(
  # add border 1)
  panel.border = element_rect(colour = "black", fill = NA, linetype = 1, linewidth = 0.6),
  # color background 2)
  panel.background = element_rect(fill = NA),
  # modify grid 3)
  panel.grid.major.x = element_line(colour = "black", linetype = 3, size = 0.5),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, size = 0.5),
  panel.grid.minor.y = element_blank(),
  # modify text, axis and colour 4) and 5)
  axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 12),
  axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 14),
  axis.ticks = element_line(colour = "black"),
  #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
  # legend at the bottom 6)
  legend.position = "bottom", 
  legend.title = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
  legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 12), 
  
  strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
  strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 14),
  panel.spacing.y = unit(1, "lines"),
  panel.ontop = F)
  }

# Import data ------------------------------------------------------------------
richness_patterns <- read.csv("output/Overall_richness.csv")
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/output/families/", full.names = TRUE)
#Read them in a list
samplelist <- lapply(list_data, read.table)

#make big dataframe
data_out_raw <- as.data.frame(do.call(rbind, samplelist))

?
data_out <- data_out_raw %>% 
  filter(!family == "Acanthaceae" | !sp== 5422)

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

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

cont_top5 <- threshold_cont %>%
  ungroup() %>% 
  slice_max(ratio, by = c(Continent,family)) %>%
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:5) 
  
overall_top5 <- threshold_cont %>%
  ungroup() %>%
  filter(Continent == "OVERALL") %>% 
  slice_max(ratio, by = c(Continent,family)) %>%
  arrange(desc(ratio)) %>%
  slice(1:5) 
  
  


#Plotting-----------------------------------------------------------------------
pallete <- as.vector(GetColors(5, scheme = "light"))
pallete <- wes_palette(name  = "BottleRocket2", 5, type = c("discrete"))

samples_threshold <- samples_cumulative_rel %>% 
  filter(family %in% overall_top5$family) %>%
  ggplot(aes(x= log10(sp), y=mean, col = family)) +
  geom_line(lwd = 1)  +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  facet_wrap(~Continent) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1, ) +
  scale_color_discrete("Family") +
  scale_color_nejm() +
  theme_bw() +
  theme(panel.border = element_blank()) 
samples_threshold

  
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
plotting.lines("OVERALL")

x <- unique(threshold_cont$Continent)

lapply(x, plotting.lines)


plotting.lines.percont <- function (x) {
  overall_top5 <- threshold_cont %>%
    ungroup() %>%
    filter(Continent == x) %>% 
    slice_max(ratio, by = c(Continent,family)) %>%
    arrange(desc(ratio)) %>%
    slice(1:5) 
  
  plot <- samples_cumulative_rel %>% 
    filter(family %in% overall_top5$family & Continent == x[1]) %>%
    ggplot(aes(x= log10(sp), y=mean, col = family)) +
    geom_line(lwd = 1)  +
    xlab("Sample size species (log10)") +
    ylab("Correlation coefficient") +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
                linetype = 0, alpha = 0.1, ) +
    ylim(0,1) +
    facet_wrap(~Continent) +
    scale_color_discrete("Family", type = pallete) +
    theme(
      # add border 1)
      panel.border = element_rect(colour = "black", fill = NA, linetype = 1, linewidth = 0.6),
      # color background 2)
      panel.background = element_rect(fill = NA),
      # modify grid 3)
      panel.grid.major.x = element_line(colour = "black", linetype = 3, size = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, size = 0.5),
      panel.grid.minor.y = element_blank(),
      # modify text, axis and colour 4) and 5)
      axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 12),
      axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 14),
      axis.ticks = element_line(colour = "black"),
      #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
      # legend at the bottom 6)
      legend.position = "right", 
      legend.title = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
      legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 12), 
      
      strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
      plot.title = element_text(colour = "black",face = "plain", family = "sans", size = 16, hjust = 0.5),
      panel.spacing.y = unit(1, "lines"),
      panel.ontop = F
    )
  ggsave(paste0("family_curves_top_", x, 
                "_increasingsteps_logsp_newtheme_solo.jpeg"),  
         width = 25, height = 15, dpi = 600)
}

x <- unique(threshold_cont$Continent)


plotlist <- lapply(x, plotting.lines.percont)


ggarrange(plotlist = plotlist,
          ncol = 5,
          nrow = 2, 
          labels = "AUTO",
          align = "hv"
          )

ggsave(paste0("family_curves_top_", "individual", 
              "_increasingsteps_logsp_newtheme.jpeg"),  
       width = 40, height = 15, dpi = 600)


samples_threshold_cont <- samples_cumulative_rel %>% 
  filter(family %in% overall_top5$family) %>%
  ggplot(aes(x= log10(sp), y=mean, col = family)) +
  geom_line(lwd = 1)  +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  facet_wrap(~Continent) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1, ) +
  theme_bw() +
  theme(panel.border = element_blank())
samples_threshold_cont

ggsave(paste0("family_curves_top_", nrow(threshold), 
              "_increasingsteps_logsp.jpeg"),  
       width = 25, height = 15, dpi = 600)




# richness maps 
rich_overall_con <- richness_patterns %>% 
  filter(ID =="con") %>% 
  left_join(tdwg_1_names, by =c("LEVEL_COD" = "LEVEL1_COD"))

rich_species_3001 <- samples_cumulative_rich %>% 
  filter(sp == 3001) %>% 
  group_by(LEVEL_COD) %>% 
  summarise(mean_richness_sample = mean(richness_sample),
            richness_overall = mean(richness), 
            LEVEL1_COD = unique(LEVEL1_COD))


tdwg_3_mod <- tdwg_3 %>% 
  filter(LEVEL3_COD %in% rich_species_3001$LEVEL_COD) %>% 
  left_join(rich_species_3001, by = c("LEVEL3_COD" = "LEVEL_COD")) 

pattern <- ggplot(data = tdwg_3_mod) +
  geom_sf(aes(fill = richness_overall)) +
  scale_fill_gradient(low = "grey90", high = "grey1", na.value = "white") +
  theme_bw()

sample <- ggplot(data = tdwg_3_mod) +
  geom_sf(aes(fill = mean_richness_sample)) +
  scale_fill_gradient(low = "grey90", high = "grey1", na.value = "white")+ 
  theme_bw()

ggarrange(pattern, sample, ncol =2, nrow =1)
ggsave(paste0("maps_", dataset, 
              "_increasingsteps",diff(rich_rel_cumulative_names$sp)[1],".jpeg"),  
       width = 15, height = 8, dpi = 600)

tdwg_1_mod <- tdwg_1 %>% 
  filter(LEVEL1_COD %in%richness_patterns_con$LEVEL1_COD) %>% 
  left_join(richness_patterns_con, by = "LEVEL1_COD") 
  

con <- ggplot(data = tdwg_1_mod) +
  geom_sf(aes(fill = richness)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick1", na.value = "white")

tdwg_2_mod <- tdwg_2 %>% 
  filter(LEVEL2_COD %in%richness_patterns_reg$LEVEL2_COD) %>% 
  left_join(richness_patterns_reg, by = "LEVEL2_COD") 
reg <- ggplot(data = tdwg_2_mod) +
  geom_sf(aes(fill = richness)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick2", na.value = "white")

ggarrange(bru, reg, con,ncol =3, nrow =1)



rich_species_3001 <- samples_cumulative_rich %>% 
  filter(sp == 3001) %>% 
  ggplot( aes(x =Continent, y =richness_sample, fill = Continent ))+
  geom_jitter(color="black", size=0.4, alpha=0.09) +
  geom_boxplot() +
  theme_bw()

rich_species_full <- samples_cumulative_rich %>% 
  filter(sp == 3001) %>% 
  ggplot(aes(x =Continent, y =richness, fill = Continent))+
  geom_boxplot() +
  theme_bw()

ggarrange (rich_species_full,rich_species_3001, ncol = 2, nrow =1)

ggsave(paste0("speciesboxplot_",dataset, "_steps",
              diff(rich_rel_cumulative_names$sp)[1],".jpeg"),  
       width = 25, height = 8, dpi = 600)
