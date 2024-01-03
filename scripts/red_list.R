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
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 14),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    #axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "right", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 14), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
    
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
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, linewidth =  0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 14),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black", linetype = 1, linewidth =  0.5),
    # legend at the bottom 6)
    legend.position = "bottom", 
    legend.title = element_text(colour = "black", face = "bold", family = "sans", size = 14), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
    
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 14),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F,
    legend.key.size = unit(1, "lines")
  )
}

redlist_analysis <- read.csv("data/redlist/assessments.csv")
richness_patterns_allplants <- read.csv("output/Overall_richness.csv")


plants_full_all <- fread("data/wcvp_accepted_merged.txt")
dist_native_all <- fread("data/dist_native.txt")


tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")




tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

plants_full <- plants_full_all %>% 
  filter(taxon_name %in% redlist_analysis$scientificName) %>% 
  left_join(redlist_analysis, by = c("taxon_name" = "scientificName")) 

a <- as.data.frame(table(plants_full$family))

dist_native <- dist_native_all %>% 
  filter(plant_name_id %in% plants_full$plant_name_id)



richness_patterns_con <- dist_native %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- dist_native %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- dist_native %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

richness_patterns1 <- richness_patterns_con %>% 
  rename( LEVEL_COD = LEVEL1_COD)
richness_patterns2 <- richness_patterns_reg %>% 
  rename( LEVEL_COD = LEVEL2_COD)
richness_patterns3 <- richness_patterns_bru %>% 
  rename( LEVEL_COD = LEVEL3_COD)

richness_patterns <- rbind(richness_patterns1,richness_patterns2, 
                           richness_patterns3)


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


tdwg_3_mod_richcomp <- tdwg_3_mod %>% 
  left_join(richness_patterns_allplants, by = c("LEVEL3_COD" = "LEVEL_COD"), multiple = "all") 


family_rich_red <- plants_full %>%  
  group_by(genus) %>% 
  summarise(n_red= n_distinct(plant_name_id))


family_rich_full <- plants_full_all %>%  
  group_by(genus) %>% 
  summarise(n_full= n_distinct(plant_name_id))

family_rich_comp <- family_rich_full %>%  
  left_join(family_rich_red, by = "genus") %>% 
  replace(is.na(.), 0)  
  
not_assessed <- family_rich_comp  %>% 
  filter(n_red == 0)

assessed <- family_rich_comp  %>% 
  filter(n_red == 0)

hist(sqrt(family_rich_comp$n_red))



p1 <- lm(family_rich_comp$n_full ~ family_rich_comp$n_red)



cor.test(log10(family_rich_comp$n_full), log10(family_rich_comp$n_red), family = "pearson")

cor.test(tdwg_3_mod_richcomp$richness.x, tdwg_3_mod_richcomp$richness.y, family = "spearman")


plot(log10(family_rich_comp$n_red) ~ log10(family_rich_comp$n_full))


richness_patterns <- read.csv("output/Overall_richness.csv")

# names 
plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

# distribution
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

# Region names 
tdwg_codes <- as.data.frame(tdwg_3) %>% 
  select(LEVEL3_COD,LEVEL2_COD,LEVEL1_COD)

# Overall patterns at botanical recording unit scale  
rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))


# Subsampling process
start.time <- Sys.time() # time keeping 

# empty objects to fill 
cumulative_namelist <- data.frame()
rich_rel_cumulative <- data.frame()
rich_cumulative <- data.frame()
samplelist <- list()
samples <- list()

for (x in 1:100) {
  repeat{
    # randomly sampling the data for species ids 
    plantlist_names_left <- plantlist_names %>% 
      filter(!plant_name_id %in%  cumulative_namelist$plant_name_id)
    # stepwise increase of sample size (cuts down running time)
    if (nrow(cumulative_namelist) < 10) {
      species_sample <- sample_n(plantlist_names_left, 1)
    }
    if (nrow(cumulative_namelist) > 10 & nrow(cumulative_namelist) < 100) {
      species_sample <- sample_n(plantlist_names_left, 10)
    }
    if (nrow(cumulative_namelist) > 100 & nrow(cumulative_namelist) < 1000){
      species_sample <- sample_n(plantlist_names_left, 100)
    } 
    if (nrow(cumulative_namelist) > 1000 & nrow(plantlist_names_left) > 1000){
      species_sample <- sample_n(plantlist_names_left, 1000)
    } 
    if (nrow(plantlist_names_left) < 1000){
      species_sample <- plantlist_names_left
    }
    if (nrow(plantlist_names_left) == 0){
      break
    }
    # create cumulative dataframe with all names that have been sample
    cumulative_namelist <- rbind(cumulative_namelist, species_sample)
    species <- cumulative_namelist
    
    # extracting distribution of sample 
    dist <- dist_native %>% 
      filter(plant_name_id %in% species$plant_name_id)
    
    # richness patterns across botanical recording units 
    sample_rich_bru <- dist %>% 
      group_by(area_code_l3) %>% 
      summarise(richness_sample = n_distinct(plant_name_id)) %>% 
      rename(LEVEL_COD = area_code_l3)
    
    # combing with datasets with overall patterns 
    rich_rel <- rich_overall_bru %>% 
      left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
      replace(is.na(.), 0) %>% # replacing NA with 0 (no correlation)
      mutate(sp = nrow(cumulative_namelist))
    
    rich_cumulative <- rbind(rich_cumulative, rich_rel)
    
    # measuring correlation between richness patterns 
    # overall
    rich_rel_bru <- rich_rel %>% 
      summarise(cor.sp = cor.test(richness_sample, richness,
                                  method="spearman", exact =F)[[4]],
                cor.pe = as.numeric(cor(richness_sample, richness, 
                                        method="pearson"))) %>% 
      mutate(id = "overall",
             sp = nrow(cumulative_namelist))
    
    # per continent 
    cor.sp <- rich_rel %>% 
      split(.$LEVEL1_COD) %>% 
      map_dbl(~cor.test(.$richness_sample, .$richness, 
                        method="spearman", data = ., exact =F)[[4]])
    cor.pe <- rich_rel %>% 
      split(.$LEVEL1_COD) %>% 
      map_dbl(~cor.test(.$richness_sample, .$richness, 
                        method="pearson", data = ., exact =F)[[4]])  
    
    # combining it all
    spear <- as.data.frame(cor.sp) %>% 
      mutate(id = rownames(.))
    rich_rel_con <- as.data.frame(cor.pe) %>% 
      mutate(id = rownames(.),
             sp = nrow(cumulative_namelist)) %>% 
      left_join(spear, by = "id")
    
    # cumulative patterns
    rich_rel_cumulative <- rbind(rich_rel_cumulative,rich_rel_bru, 
                                 rich_rel_con)
    if (nrow(cumulative_namelist) == nrow(plantlist_names)){
      break # stop when all names have been sampled 
    }
    if (nrow(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
      print(paste0("There are ", nrow(cumulative_namelist) ,
                   " species in the subsample")) # print subsample size 
    }
  }
  # extracting results to objects saved outside
  samplelist[[x]] <- rich_rel_cumulative 
  samples[[x]] <- rich_cumulative
  cumulative_namelist <- data.frame()
  rich_rel_cumulative <- data.frame()
  rich_cumulative <- data.frame()
  print(paste0("This is iteration ", x))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}



data_out_red <- as.data.frame(do.call(rbind, samplelist))

samples_cumulative_rel_red <- data_out_red %>%
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

threshold_cont_red <- samples_cumulative_rel_red %>%
  mutate(threshold = case_when(mean > 0.8 & sp > 100 ~ 'above', 
                               mean < 0.8 | sp < 100 ~ 'below')) %>% 
  filter(min(sp) & threshold == "above") %>% 
  mutate(ratio = mean/sp) %>% 
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:1) 

threshold_red <- samples_cumulative_rel_red %>% 
  group_by(Continent) %>% 
  slice_max(mean) 
  filter(min(sp) & threshold == "above")

colouring <- c( '#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', "#CC3311", "#00183E")
  
curves_red <- ggplot(data = samples_cumulative_rel_red, aes(x=log10(sp), y=mean, col = Continent)) +
  geom_line(lwd = 1, alpha = 0.5)  +
  geom_line(data = subset(samples_cumulative_rel_red, Continent == 'OVERALL'), 
            lwd = 2, color = "midnightblue", alpha = 0.7) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  xlim(c(0,6)) +
  scale_colour_manual(values = c(colouring)) +
  theme_default.lb()

curves_red



ggsave(file = paste0("curves_redlist_cropdiv.png"),  curves_red, 
       width = 12, height = 7, dpi = 600)



red_box <- data_out_red %>% 
  filter(id == "overall") %>% 
  slice_max(sp)


sample_extract <- ggplot(red_box, aes(x = id, y =cor.sp)) +
  geom_jitter(color = "black", alpha = 0.05) +
  labs(x = "Continent", y = "Correlation coefficient") +
  stat_boxplot(geom = "errorbar", width = 0.2, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.95, width =0.5, color = "grey10")+
  theme_default.lb.boxplot()
# base -------------------------------------------------------------------------

data_out <- read.table("output/full_samples_rel_steps_increasing_niter100.txt")


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

threshold_cont <- samples_cumulative_rel %>%
  mutate(threshold = case_when(mean > 0.95 & sp > 100 ~ 'above', 
                               mean < 0.95 | sp < 100 ~ 'below')) %>% 
  filter(min(sp) & threshold == "above") %>% 
  mutate(ratio = mean/sp) %>% 
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:1) 

threshold <- samples_cumulative_rel %>% 
  mutate(threshold = case_when(mean > 0.95  ~ 'above', 
                               mean < 0.95 ~ 'below')) %>% 
  group_by(Continent) %>% 
  filter(min(sp) & threshold == "above") %>%
  slice_min(sp) 

threshold_red$diff <- threshold_red$sp - threshold$sp

steps <- as.data.frame(unique(data_out_red$sp))

steps_10 <- steps %>% 
  slice_tail(n = 10) %>% 
  rename(sp = `unique(data_out_red$sp)`)

red_box_data <- data_out_red %>% 
  filter(id == "overall") %>% 
  mutate(id = "overall redlist") %>% 
  filter(sp %in% steps_10$sp) 



all_box_data <- data_out %>% 
  filter(id == "overall") %>% 
  filter(sp %in% steps_10$sp) 
  


all_box <- ggplot(all_box_data, aes(x = id, y =cor.sp, fill = "orange")) +
  geom_jitter(color = "black", alpha = 0.05) +
  labs(x = "Dataset", y = "Correlation coefficient") +
  stat_boxplot(geom = "errorbar", width = 0.2, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.95, width =0.2, color = "grey10") +
  scale_fill_manual(values = "grey") +
  theme_default.lb.boxplot()


red_box <- ggplot(red_box_data, aes(x = id, y =cor.sp, fill = "orange")) +
  geom_jitter(color = "black", alpha = 0.05) +
  labs(x = "Dataset", y = "Correlation coefficient") +
  stat_boxplot(geom = "errorbar", width = 0.2, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.95, width =0.2, color = "grey10")+
  scale_fill_manual(values = "darkred") +
  theme_default.lb.boxplot()


ggarrange(all_box,red_box, ncol = 2, nrow =1, legend = "none", common.legend = T)


ggsave(file = paste0("boxes_redlist_cropdiv_last10.png"),
       width = 12, height = 7, dpi = 600)
