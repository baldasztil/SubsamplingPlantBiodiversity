# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------
# This is a script to create random subsamples from the WCVP. It is using the
# taxonomic and geographic information available through the data 
# to relate species richness in the subsample to global patterns.  


# Packages ---------------------------------------------------------------------


library(tidyverse)
library(data.table)
library(sf)
library(tmap)
library(spdep)
library(tmaptools)


# Data -------------------------------------------------------------------------

richness_patterns <- fread("data/richness_patterns.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
tdwg_3_adj <- st_read(dsn = "data/wgsrpd-master/level3_adj")
#adjusting the plotting file 
tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)
tdwg_3_mod <- tdwg_3_adj %>% 
  left_join(richness_patterns, by = c("LEVEL3_COD" = "LEVEL_COD"), multiple = "all") %>% 
  left_join(tdwg_1_names, by = c("LEVEL1_COD"), multiple = "all") %>% 
  filter(!is.na(richness))


# Options ----------------------------------------------------------------------
sf_use_s2(FALSE)
tmap_options(check.and.fix = TRUE)


# mapping projections available in r
proj_db <- system.file("proj/proj.db", package = "sf")
crs_table <- sf::read_sf(proj_db, "crs_view") 

# choosing a coordinate system with epsg codes 
crs <- st_crs("ESRI:53042") # winkel 3


# Analysis ---------------------------------------------------------------------




# richness patterns 
tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1_names$LEVEL1_COD)
rich_overall_con <- richness_patterns %>% 
  filter(ID =="con") %>% 
  left_join(tdwg_1_names, by =c("LEVEL_COD" = "LEVEL1_COD")) %>% 
  mutate(sum = nrow(plants_full),
  freq = richness / sum)

# projecting and validating spatial files 
tdwg_3_wintri <- st_transform(tdwg_3_mod, crs = "ESRI:54012")
tdwg_winktri <- st_make_valid(tdwg_3_wintri)
tdwg_wgs <- st_make_valid(tdwg_3_mod) %>% 
  distinct(LEVEL3_NAM, richness, .keep_all= TRUE)

# checking geometries
table(st_is_valid(tdwg_3_mod))
tdwg_needadj <- tdwg_3_mod[!st_is_valid(tdwg_3_mod)== TRUE,]

tmap_mode("plot") # can also specify "view"
# mapping normal species map 
png(filename = "world.png",width = 12, height = 8, units = "in", res = 600)

tm_shape(tdwg_wgs) + 
  tm_fill(col = "richness", style="fixed", n=20, palette="Greens",
                             breaks = c(seq(0, 20000, by=2000), Inf),   legend.show = FALSE) +
  tm_legend(outside=TRUE) +
  #tm_graticules(col = "grey", alpha = 0.7, labels.show = F) +
  #tm_compass(type = "4star", size = 1, position = c("left", "top")) +
  #'tm_scale_bar(breaks = c(0, 1500, 3000), text.size = 4) +
  tm_layout(frame = FALSE)

dev.off()


# a bubble map of species richness per continent 
tmaptools::palette_explorer()
cols <- get_brewer_pal("Paired", n = 9, stretch = T)
tm_shape(tdwg_winktri) +
  tm_polygons(col = "LEVEL1_NAM", palette = cols, alpha = 1
              , contrast=.7, id="LEVEL1_NAM", title="Continent", 
              border.col = "grey1", border.alpha = 0.7) +
  tm_bubbles("richness", col = "richness", alpha = 0.9,
             border.col = NA,border.lwd = 0, border.alpha = .5, 
             style="fixed", breaks=c(seq(0, 10000, by=2000), Inf),
             palette="YlOrRd", contrast=1, 
             title.size="Sp richness", 
             title.col="Sp richness", id="name", 
             legend.shape.is.portrait = F) +
  tm_graticules(col = "grey", alpha = 0.7) +
  tm_compass(type = "4star", size = 2, position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 1500, 3000), text.size = 0.7) +
  tm_layout(frame = FALSE) 






# testing for spatial autocorrelation

tdwg_3_mod$sqrtrich <- sqrt(tdwg_3_mod$richness)

hist(sqrt(tdwg_3_mod$richness))


boxplot(tdwg_3_mod$sqrtrich, horizontal = TRUE)


tdwg_3_nb <- poly2nb(tdwg_3_mod, queen=TRUE)


tdwg_3_lw <- nb2listw(tdwg_3_nb, style="W", zero.policy=TRUE)


moran.test(tdwg_3_mod$richness,tdwg_3_lw, alternative="greater", zero.policy = TRUE)


MC<- moran.mc(tdwg_3_mod$richness,tdwg_3_lw, nsim=9099, alternative="greater", zero.policy = TRUE)
plot(MC)

 


a<- is.na(tdwg_3_mod$richness)

 table(a)


 tdwg_3_mod
 
 set.seed(131)
 tdwg_3_mod$rand1 <- sample(tdwg_3_mod$richness, length(tdwg_3_mod$richness), replace = FALSE)
 tdwg_3_mod$rand2 <- sample(tdwg_3_mod$richness, length(tdwg_3_mod$richness), replace = FALSE)
 tdwg_3_mod$rand3 <- sample(tdwg_3_mod$richness, length(tdwg_3_mod$richness), replace = FALSE)
 
 tm_shape(tdwg_3_mod) + tm_fill(col=c("richness", "rand1", "rand2", "rand3"),
                       style="quantile", n=8, palette="Greens", legend.show = FALSE) +
   tm_facets( nrow=1)

 
 