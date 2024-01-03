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

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  dplyr::select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

plants_full_all <- fread("data/wcvp_accepted_merged.txt")
dist_native_all <- fread("data/dist_native.txt")

plantlist_names <- plants_full_all %>%  
  dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native_all %>% 
  filter(plant_name_id %in% plants_full_all$plant_name_id)

samples_cumulative_rel <- fread("output/cumulative_rel/fullsamples_cumulative_rel.txt") %>% 
  mutate(method = "random",
         dataset = "World Checklist") %>% 
  dplyr::select(id, method, sp, cor = mean, dataset)

samples_cumulative_rel_red <- fread("output/cumulative_rel/redlist_cumulative_rel.txt") %>%  
  mutate(method = "random", 
         dataset = "IUCN Red List") %>% 
  dplyr::select(id, method, sp, cor = mean, dataset)

samples_cumulative_rel_srli <- fread("output/cumulative_rel/srli_cumulative_rel.txt") %>% 
  mutate(method = "random", 
         dataset = "Sampled Red List Index") %>% 
  dplyr::select(id, method, sp, cor = mean, dataset)

samples_cumulative_rel_red_optimisedim  <- fread("output/optimised/DEoptim_results_redlist.txt") %>% 
  mutate(method = "optimised", 
         id = "overall",
         dataset = "IUCN Red List") %>% 
  dplyr::select(id, method, sp = Species, cor = best_corr, dataset)

samples_cumulative_rel_srli_optimisedim <- fread("output/optimised/DEoptim_results_srli.txt") %>% 
  mutate(method = "optimised", 
         id = "overall", 
         dataset = "Sampled Red List Index") %>% 
  dplyr::select(id, method, sp = Species, cor = best_corr, dataset)

samples_cumulative_rel_optimisedim <- fread("output/optimised/DEoptim_results_1000.txt") %>% 
  mutate(method = "optimised", 
         id = "overall", 
         dataset = "World Checklist") %>% 
  dplyr::select(id, method, sp = Species, cor = best_corr, dataset)

samples_comb <- rbind(samples_cumulative_rel, samples_cumulative_rel_red, samples_cumulative_rel_srli, 
                      samples_cumulative_rel_optimisedim, samples_cumulative_rel_red_optimisedim, samples_cumulative_rel_srli_optimisedim)
samples_rand <- rbind(samples_cumulative_rel, samples_cumulative_rel_red, samples_cumulative_rel_srli)
unique(samples_comb$dataset)

xx <- samples_comb %>% 
  filter(id == "overall") %>% 
  filter(cor > 0.95) %>% 
  group_by(method, dataset) %>%
  reframe(minsp = min(sp))

xx <- samples_comb %>% 
  filter(id == "overall") %>% 
  group_by(method, dataset) %>%
  reframe(minsp = min(cor))

overall_patterns <-samples_comb %>% 
  filter(id == "overall") %>% 
  filter(sp > 1389 & sp < 5124)


patterns_rearranged <- overall_patterns %>%
  arrange(desc(cor)) %>% 
  group_by(dataset,method) %>% 
  group_modify(~ add_row(.x,.before=0)) %>% 
  mutate(cor = ifelse(is.na(cor), 0, cor),
         sp = ifelse(is.na(sp), 0, sp),
         id = ifelse(is.na(id), "overall", id))  %>% 
  mutate(dataset = as.factor(dataset))
patterns_rearranged$dataset <-  factor(patterns_rearranged$dataset, levels=c('IUCN Red List', 'Sampled Red List Index', 'World Checklist'))

patterns_rearranged <- overall_patterns

# Reorder data
ggplot(patterns_rearranged, aes(x=dataset, y= cor, fill= dataset)) +
  geom_violin(aes(linetype = method), width=2,  alpha = 0.7, position = position_dodge(0.9), lwd = 0.75, show.legend = F) +
  geom_boxplot(aes(col = method), fill = "white", width=0.15, outlier.alpha = 0.15, outlier.size = 0, alpha = 0.3, position = position_dodge(0.9),  show.legend = F) +
  scale_fill_manual(values = c("darkred","red","midnightblue")) +
  scale_colour_manual(values = c("black","black","black")) +
  theme_bw() + xlab("Dataset") +
  ylab("Correlation coefficient") + ylim(c(0.8,1)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.4), lwd = 1) 

#stat_summary(
#  fun.data = "mean_se",  fun.args = list(mult = 1), 
#  geom = "pointrange", color = "black", 
#  position = position_dodge(0.9)) +

#facet_wrap(~method)

ggsave(file = paste0("violins_redlist_srli__fullsamples_1390to5123_optimised_random_adjnames.png"),
       width = 6, height = 5, device = "png", dpi = 1200)

t.test(cor ~ method, data = patterns_rearranged) 
pairwise.t.test(cor,factor(patterns_rearranged$dataset), data = patterns_rearranged, p.adjust.method = "bonferroni")
patterns_rearranged %>%
  group_by(dataset) %>%
  get_summary_stats(cor, type = "mean_sd")

stats <- patterns_rearranged %>%
  group_by(dataset, method) %>%
  get_summary_stats(cor, type = "full")

anova_model <- aov(cor ~ dataset, data = patterns_rearranged)
TukeyHSD(anova_model, conf.level=.95) 


curves <- 
  ggplot(data = patterns_rearranged, aes(x=log10(sp), y= cor, col = dataset, linetype = method)) +
  geom_line(lwd = 1.6, alpha = 0.8)  +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  # geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
  #            linetype = 0, alpha = 0.1) +
  scale_color_manual(values = c("darkred","red","midnightblue")) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.4), lwd = 1.5) +
  xlim(c(0,4)) +
  theme_default.lb() 

curves

ggsave(file = paste0("curves_redlist_srli__fullsamples_5123_optimised_random.png"),
       width = 15, height = 12, dpi = 600)






threshold_cont <- samples_cumulative_rel %>%
  filter(mean > 0.95) %>% 
  mutate(ratio = mean/sp) %>% 
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:1) 


threshold_cont_red <- samples_cumulative_rel_red %>%
  filter(mean > 0.95) %>% 
  mutate(ratio = mean/sp) %>% 
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:1) 

threshold_cont_srli <- samples_cumulative_rel_srli %>%
  filter(mean > 0.95) %>% 
  mutate(ratio = mean/sp) %>% 
  group_by(Continent) %>% 
  arrange(desc(ratio)) %>%
  slice(1:1) 





