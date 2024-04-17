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
library(rWCVP)
library(MKinfer)
library(arrow)
library(plotrix)                              


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
  GLOBAL_top5 <- threshold_cont %>%
    ungroup() %>%
    filter(Continent == x) %>% 
    slice_max(ratio, by = c(Continent,family)) %>%
    arrange(desc(ratio)) %>%
    slice(1:5) 
  
  plot <- samples_cumulative_rel %>% 
    filter(family %in% GLOBAL_top5$family) %>%
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

options(digits = 5)  
# Import data ------------------------------------------------------------------

#data_out2 <- read.table("output/full_samples_steps_increasing_niter100.txt")
richness_patterns <- read.csv("output/overall_richness.csv")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")

tdwg_3 <- wgsrpd3
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/fullsamples_test", full.names = TRUE)
#Read them in a list
samplelist <- lapply(list_data, read_parquet)

#make big dataframe
data_out <- rbindlist(samplelist) %>% 
  filter(sp %in% c(1:300)) 
  #as.data.frame(do.call(rbind, samplelist)) 
gc()
setDT(data_out)
# Reframe ----------------------------------------------------------------------


xx <- data_out[, list("mean_corgw" = mean(cor.sp_gwr_mean), 
                      "sd_corgw_mean" = mean(cor.sp_gwr_sd, na.rm =T), 
                      "se_corgw" = std.error(cor.sp_gwr_mean, na.rm =T),
                      
                      "mean_pseudor" = mean(pseudo.r.squared),
                      "sd_pseudor" = sd(pseudo.r.squared, na.rm =T),
                      "se_pseudor" = std.error(pseudo.r.squared, na.rm =T),
                      
                      "n" = length(cor.sp_gwr_mean),
                      
                      "lower_ci_corgw" = cor.sp_gwr_mean - error(cor.sp_gwr_mean), 
                      "upper_ci_corgw" = cor.sp_gwr_mean + error(cor.sp_gwr_mean),
               
                      "lower_ci_pseudor" = cor.sp_gwr_mean - error(pseudo.r.squared), 
                      "upper_ci_pseudor" = cor.sp_gwr_mean + error(pseudo.r.squared)),
      by = c("continent", "sp")] 
  
samples_cumulative_rel <- data_out %>%
  group_by(continent,sp) %>%  
  mutate(se = sd / sqrt(n.sample),
         lower.ci = mean - qt(1 - (0.05 / 2), 
                              n.sample - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.sample - 1) * se)  
  
  reframe(across(
    .cols = c(starts_with(c("pseudo", "cor")), 
                          ends_with("res")), 
    .fns = list(mean = mean, sd = sd(na.rm = T)), 
    .names = "{col}_{fn}"
  ), 
  n.sample = n()) 

error <- function(x) {
  qnorm(0.975)*sd(x)/sqrt(length(x))
} 

upper.ci <- function(x) {
  mean(x) + qt(1 - (0.05 / 2), 
               x$n.sample - 1) * std.error(x)
} 

head(samples_cumulative_rel)

 


 unique(samples_cumulative_rel$id)
samples_cumulative_rich <- data_out2 %>% 
  mutate(LEVEL1_COD = as.character(LEVEL1_COD)) %>% 
  left_join(tdwg_1_names, by = "LEVEL1_COD") %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "GLOBAL", Continent))



results <- read.csv("fullsamples_test.csv") %>% 
  filter(id == "overall")


samples_total <- data_out %>%  
  mutate(LEVEL1_COD = as.character(id)) %>% 
  left_join(tdwg_1_names, by ="LEVEL1_COD") %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "GLOBAL", Continent))
  
  
#Output analysis----------------------------------------------------------------
samples_cumulative_rel <- fread("output/cumulative_rel/fullsamples_cumulative_rel.txt")

threshold <- samples_cumulative_rel %>% 
 filter(mean > 0.95) %>% 
  group_by(Continent) %>% 
  mutate(ratio = mean/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 

#Summary statistics analysis----------------------------------------------------


richness_patterns_con <- dist_native %>%
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  mutate(LEVEL_COD = as.character(continent_code_l1)) %>%  
  mutate(ID = "con")


index_endemic <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n == 1) # 56%

richness_patterns_con_en <- dist_native %>%
  filter(plant_name_id %in% index_endemic$plant_name_id) %>% 
  group_by(continent_code_l1) %>% 
  summarise(rich_en = n_distinct(plant_name_id))%>% 
  mutate(LEVEL_COD = as.character(continent_code_l1)) %>%  
  mutate(ID = "con")

index_widespread <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) # 43.9%

richness_patterns_con_nonen <- dist_native %>%
  filter(plant_name_id %in% index_widespread$plant_name_id) %>% 
  group_by(continent_code_l1) %>% 
  summarise(rich_wide = n_distinct(plant_name_id))%>% 
  mutate(LEVEL_COD = as.character(continent_code_l1)) %>% 
  mutate(ID = "con")


richvsamplesize <- threshold %>%
  filter(!Continent == "GLOBAL") %>% 
  left_join(richness_patterns, by = join_by("id" == "LEVEL_COD")) %>%
  left_join(richness_patterns_con_en, by = join_by("id" == "LEVEL_COD")) %>%
  left_join(richness_patterns_con_nonen, by = join_by("id" == "LEVEL_COD")) %>% 
  select(Continent, sp, richness, rich_en, rich_wide) %>% 
  mutate(prop_en = rich_en / (rich_en + rich_wide), 
         prop_wide = rich_wide / (rich_en + rich_wide)) %>% 
  ungroup() %>%
  mutate(sp = sp,
         ratio =   log10(sp / richness),
         rank_rich = rank(richness),
         rank_en = rank(rich_en),
         rank_wide = rank(rich_wide),
         rank_threshold =  rank(sp))



GLOBAL_patterns <- dist_native %>% 
  summarise(Continent = "GLOBAL", 
            sp = as.numeric(threshold[8,2]),
            richness = n_distinct(plant_name_id),
            rich_en = n_distinct(index_endemic$plant_name_id),
            rich_wide = n_distinct(index_widespread$plant_name_id),
            prop_en = rich_en / (rich_en + rich_wide), 
            prop_wide = rich_wide / (rich_en + rich_wide),
            ratio =  log10(sp / richness)) %>%
  bind_rows(richvsamplesize) %>% 
  replace(is.na(.), 0) 

rich_summary <- GLOBAL_patterns %>% 
  rename(threshold = sp, ratio_rich_tresh = ratio) 

rich_table <- reactable(
  rich_summary,
  defaultSorted = "threshold",
  defaultSortOrder = "desc",
  defaultColDef = colDef(
    cell = data_bars(rich_summary, 
                     #number_fmt = scales::number_format(accuracy = 0.01),
                     fill_color = rev(viridis::viridis(100)), 
                     #background = "grey",
                     text_position = "above", 
                     bar_height = 6, 
                     #fill_gradient = TRUE
                     box_shadow = T), align = "center", width = 150
    
  ),
    columns = list(
      Continent = colDef(show = T),
      threshold = colDef( show = T, format = colFormat(digits = 0)),
      richness = colDef( show = T, format = colFormat(digits = 0)), 
      rich_en = colDef( show = T, format = colFormat(digits = 0)),
      rich_wide = colDef( show = T, format = colFormat(digits = 0)),
      prop_en = colDef( show = T, format = colFormat(digits = 2)),
      prop_wide = colDef( show = T, format = colFormat(digits = 2)),
      ratio_rich_thresh = colDef( show = T, format = colFormat(digits = 1)),
      rank_rich = colDef( show = T, format = colFormat(digits = 0)),
      rank_en = colDef( show = T, format = colFormat(digits = 0)),
      rank_wide = colDef( show = T, format = colFormat(digits = 0)),
      rank_threshold = colDef( show = T, format = colFormat(digits = 0))
    )
)

rich_table
#save_reactable(rich_table, "rich__summary_table.png", vwidth = 1900)

#Plotting-----------------------------------------------------------------------

a <- palette_okabe_ito()






sample_rich_threshold <- as.data.frame(do.call(rbind, xx)) # from extraction script 

samples_cumulative_rich <- sample_rich_threshold %>% 
  mutate(LEVEL1_COD = as.character(LEVEL1_COD)) %>% 
  left_join(tdwg_1_names, by = "LEVEL1_COD") %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "GLOBAL", Continent))

rich_species_1390 <- samples_cumulative_rich %>% 
  group_by(Continent) %>% 
  mutate(mean = mean(richness_sample))

rich_species_1390_box <- 
  ggplot(rich_species_1390, aes(x = reorder(Continent,mean), y =richness_sample, fill = reorder(Continent,mean))) +
  geom_jitter(color = "black", alpha = 0.05) +
  labs(x = "Continent", y = "Number of species (sample)") +
  scale_fill_brewer(palette = "Reds" ) +
  stat_boxplot(geom = "errorbar", width = 0.2, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.95, width =0.5, color = "grey10")+
  theme_default.lb.boxplot()
rich_species_1390_box

rich_species_full <- samples_cumulative_rich %>% 
  group_by(Continent) %>% 
  mutate(mean = mean(richness)) 

rich_species_full_box <-  
  ggplot(rich_species_full, aes(x = reorder(Continent,mean), y =richness, fill = reorder(Continent,mean))) +
  labs(x = "Continent", y = "Number of species (observed)") +
  scale_fill_brewer(palette = "Reds" ) +
  stat_boxplot(geom = "errorbar", width = 0.2, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.95, width =0.5, color = "grey10")+
  theme_default.lb.boxplot()
rich_species_full_box




ggarrange(rich_species_full_box,rich_species_1390_box, ncol = 2, nrow =1, legend = "none", common.legend = T)
ggsave(file = paste0("boxplots", 
                     "cropdiv","@",1390,".jpeg"),  
       width = 35, height = 10, dpi = 600)



# richness maps 
rich_GLOBAL_con <- richness_patterns %>% 
  filter(ID =="con") %>% 
  left_join(tdwg_1_names, by =c("LEVEL_COD" = "LEVEL1_COD"))

rich_species_1390 <- samples_cumulative_rich %>% 
  group_by(LEVEL_COD) %>% 
  summarise(mean_richness_sample = mean(richness_sample),
            richness_GLOBAL = mean(richness), 
            LEVEL1_COD = unique(LEVEL1_COD))


tdwg_3_mod <- tdwg_3 %>% 
  filter(LEVEL3_COD %in% rich_species_1390$LEVEL_COD) %>% 
  left_join(rich_species_1390, by = c("LEVEL3_COD" = "LEVEL_COD")) 

pattern <- ggplot(data = tdwg_3_mod) +
  geom_sf(aes(fill = richness_GLOBAL)) +
  scale_fill_gradient("Species (n)", low = "grey90", high = "grey1", na.value = "white") +
  ggtitle("Species richness (observed)") +  
  theme_bw(base_size = 28) +
  theme(plot.title = element_text(hjust = 0.5))

sample <- ggplot(data = tdwg_3_mod) +
  geom_sf(aes(fill = mean_richness_sample)) +
  scale_fill_gradient("Species (n)", low = "grey90", high = "grey1", na.value = "white")+ 
  ggtitle("Species richness (sample)") +
  theme_bw(base_size = 28) +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange( pattern, sample, ncol = 2, nrow =1, common.legend = F)

ggsave(file = paste0("maps", 
                     "cropdiv","@",1390,".jpeg"),  
       width = 35, height = 10, dpi = 600)


#a <- colour("muted")
colouring <- c( '#CC6677', '#722280', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', "#CC3311", "#00183E")

viri <- viridisLite::viridis(100)


samples_cumulative_rel$Continent <- as.factor(samples_cumulative_rel$Continent)
samples_cumulative_rel$Continent<-factor(samples_cumulative_rel$Continent, levels=c("AFRICA","ANTARCTIC","ASIA-TEMPERATE","ASIA-TROPICAL",
                                            "AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC","SOUTHERN AMERICA","GLOBAL"))



curves <- ggplot(data = samples_cumulative_rel, aes(x=log10(sp), y=mean, col = Continent)) +
  geom_line(lwd = 1.5, alpha = 0.5)  +
  geom_line(data = subset(samples_cumulative_rel, Continent == 'GLOBAL'), 
            lwd = 8, color = "midnightblue", alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 1.5) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.25) +
  scale_colour_manual(values = c(colouring)) +
  xlim(c(0,6)) +
  theme_default.lb() + 
  theme(legend.background=element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 4, nrow = 3))
  
ggsave(file = paste0("curves_fullsample", 
                     "_cropdiv_it100_poster.png"),curves,   
       width = 18, height = 15, dpi = 600)


# richness maps + curves
g <- arrangeGrob(curves, arrangeGrob(pattern, sample), ncol = 2)
ggsave(file = paste0("maps_curves", 
                     "_increasingsteps","@",1390,".jpeg"),g,   
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




