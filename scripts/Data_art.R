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
library(khroma)

# Defining functions and objects -----------------------------------------------



# Import data ------------------------------------------------------------------
richness_patterns <- read.csv("output/Overall_richness.csv")
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")
data_out_fullsamp <- read.table("output/full_samples_rel_steps_increasing_niter100.txt")


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
  filter(id == "overall") 

samples_families_plot <- data_out %>% 
  #filter(id == "overall") %>% 
  group_by(family, sp) %>% 
  mutate(n = 1:1000) %>%  
  filter(n == 1:1000) %>% 
  mutate(n = as.character(n), 
         cor.sp = ifelse(is.na(cor.sp), 0, cor.sp)) #%>% 
  #summarise(mean =  mean(cor.sp))


viri <- viridis::inferno(451)

ggplot() + 
  geom_line(data = samples_families_plot, 
            aes(x=log10(sp), y= cor.sp, col = family)
            , lwd = 2, alpha = 0.4, show.legend = F, 
            lineend = "round",
            linejoin = "round") + 
  scale_color_manual(values = viri) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  theme_void()


ggsave(paste0("family_curves_artsy_allsamples", 
              "_increasingsteps_solo_combined.png"),  
       width = 27, height = 10, dpi = 600)

# Import data ------------------------------------------------------------------


samples_total <- data_out_fullsamp %>%  
  mutate(LEVEL1_COD = as.character(id)) %>% 
  left_join(tdwg_1_names, by ="LEVEL1_COD") %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "OVERALL", Continent))




samples_cumulative_rel <- data_out_fullsamp %>%
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

ribbon <- subset(samples_cumulative_rel, Continent == 'OVERALL')

samples_total_plot <- samples_total %>% 
  #filter(Continent == "OVERALL") %>% 
  group_by(sp) %>% 
  mutate(n = 1:1000) %>% 
  filter(n == 1:1000) %>%
  mutate(log_sp = log10(sp)) %>%
  mutate(n = as.character(n), 
         cor.sp = ifelse(is.na(cor.sp), 0, cor.sp)) %>% 
  ungroup()

viri <- viridisLite::viridis(1000)

smooth_rainbow <- colour("sunset")
## Start at purple instead of off-white
plot(smooth_rainbow(1000, range = c(0.25, 1)))

viri500 <- viri[500]

ggplot() + 
  geom_line(data = samples_total_plot, 
            aes(x=log10(sp), y= cor.sp, col = n), lwd = 0.5, alpha = 0.3, show.legend = F) +
  geom_ribbon(data = ribbon, 
              aes(x = log10(sp), ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.2, col = viri500) +
  #scale_color_manual(values = c(viri)) +
  scale_colour_acton(discrete = T, range = c(0.25, 1)) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  theme_void()


sample_points <- samples_total_plot %>% 
  filter(!cor.sp == 0)

ggplot() + 
  geom_point(data = sample_points, 
            aes(x=n, y= cor.sp, col = n), size = 2, alpha = 1, show.legend = F) +
  scale_colour_smoothrainbow(discrete = T, range = c(0.25, 1)) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  theme_void()

ggsave(paste0("sample_points_artsy_all_void_rainbow.png"),  
       width = 35, height = 17, dpi = 600)
