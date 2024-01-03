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

#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/fullsamples", full.names = TRUE)
#Read them in a list
samplelist <- lapply(list_data, fread)

#make big dataframe
data_out <- as.data.frame(do.call(rbind, samplelist))

#Get the path of filenames
list_data_red <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/redlist/samples/red", full.names = TRUE)
#Read them in a list
samplelist_red <- lapply(list_data_red, fread)

data_out_red <- as.data.frame(do.call(rbind, samplelist_red))

#Get the path of filenames
list_data_srli <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/redlist/samples/srli", full.names = TRUE)
#Read them in a list
samplelist_srli <- lapply(list_data_srli, fread)

data_out_srli <- as.data.frame(do.call(rbind, samplelist_srli))


plants_full_all <- fread("data/wcvp_accepted_merged.txt")
dist_native_all <- fread("data/dist_native.txt")


tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  dplyr::select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


samples_cumulative_rel <- fread("output/cumulative_rel/fullsamples_cumulative_rel.txt")
samples_cumulative_rel_red <- fread("output/cumulative_rel/redlist_cumulative_rel.txt")
samples_cumulative_rel_srli <- fread("output/cumulative_rel/srli_cumulative_rel.txt")

samples_cumulative_rel_optim  <- fread("output/optimised/DEoptim_results_redlist.txt") %>% 
  mutate(method = "opt")
samples_cumulative_rel_red_optim <- fread("output/optimised/DEoptim_results_srli.txt") %>% 
  mutate(method = "opt")
samples_cumulative_rel_srli_optim <- fread("output/optimised/DEoptim_results_1000.txt") %>% 
  mutate(method = "opt")

aa <- samples_cumulative_rel_red %>% 
  filter(id == "overall") #never

bb <- samples_cumulative_rel_srli %>% 
  filter(id == "overall") %>% 
  filter(mean > 0.95) #3707


cc <- samples_cumulative_rel %>% 
  filter(id == "overall") %>% 
  filter(mean > 0.95) #1390

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



steps <- as.data.frame(unique(samples_cumulative_rel_srli$sp))


steps_10 <- steps %>% 
  slice_tail(n = 1000) %>% 
  rename(sp = `unique(samples_cumulative_rel_srli$sp)`) 

srli_box_data <- data_out_srli %>% 
  filter(id == "overall") %>% 
  mutate(id = "srli") %>% 
  filter(sp %in% steps_10$sp) 

red_box_data <- data_out_red %>% 
  filter(id == "overall") %>% 
  mutate(id = "redlist") %>% 
  filter(sp %in% steps_10$sp) 

all_box_data <- data_out %>% 
  filter(id == "overall") %>% 
  filter(sp %in% steps_10$sp) 

box <- bind_rows(srli_box_data, red_box_data, all_box_data)
comp <- bind_rows(data_out_srli, data_out_red, data_out)

stats <- box %>% 
  group_by(id,sp) %>% 
  summarise(mean = mean(cor.sp))

stats_comp <- spread(stats,
                     key = id, 
                     value = mean)

#library(pgirmess)
#library(FSA)
#a <- aov(cor.sp ~ id, data = box)
#krusk <- kruskalmc(cor.sp  ~  id, data = box)
#TukeyHSD(a)

all_box <- ggplot(all_box_data, aes(x = id, y =cor.sp, fill = "")) +
  geom_jitter(color = "grey80", alpha = 0.01) +
  labs(x = "Dataset", y = "Correlation coefficient") +
  stat_boxplot(geom = "errorbar", width = 0.1, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.99, width =0.2, color = "grey10") +
  scale_fill_manual(values = "grey80") +
  ylim(c(0.8,1)) +
  theme_default.lb.boxplot()


red_box <- ggplot(red_box_data, aes(x = id, y =cor.sp, fill = "")) +
  geom_jitter(color = "red", alpha = 0.01) +
  labs(x = "Dataset", y = "Correlation coefficient") +
  stat_boxplot(geom = "errorbar", width = 0.1, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.99, width =0.2, color = "grey10")+
  scale_fill_manual(values = "red") +
  ylim(c(0.8,1)) +
  theme_default.lb.boxplot()


srli_box <- ggplot(srli_box_data, aes(x = id, y =cor.sp, fill = "")) +
  geom_jitter(color = "darkred", alpha = 0.01) +
  labs(x = "Dataset", y = "Correlation coefficient") +
  stat_boxplot(geom = "errorbar", width = 0.1, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.99, width =0.2, color = "grey10")+
  ylim(c(0.8,1)) +
  scale_fill_manual(values = "darkred") +
  theme_default.lb.boxplot()



ggsave(file = paste0("boxes_redlist_srli_cropdiv_last10_fullsamples.png"),
       width = 20, height = 8, dpi = 600)


bb <- samples_cumulative_rel_srli %>% 
  filter(id == "overall") %>% 
  filter (sp <= 5423)

aa <- samples_cumulative_rel_red %>% 
  filter(id == "overall") %>% #never
  filter (sp <= 5423)


cc <- samples_cumulative_rel %>% 
  filter(id == "overall") %>% 
  filter (sp > 0 & sp <= 5423) %>% 
  mutate(ID = "baseline")


threshold_comb <- bind_rows(threshold_cont, threshold_cont_red, threshold_cont_srli)
combo_raw <- bind_rows(aa, bb, cc)

#write.csv(threshold_comb, "threshold_redlist_slri_wcvp.csv")

combo <- combo_raw %>% 
  filter(sp > 0) %>% 
  select(Dataset = ID, mean, sp, lower.ci, upper.ci)


combo_red <- combo %>% 
  arrange(desc(mean)) %>% 
  group_by(Dataset) %>% 
  slice_head(n = 1000)

# insert box instead of combo for different plots 
violins <- ggplot(combo, aes(x = mean, y= Dataset, fill = Dataset)) +
  #geom_jitter(color = "darkred", alpha = 0.01) +
  labs(x = "Correlation coefficient", y = "Dataset") +
  geom_violin(lwd =0.8, scale = "width", trim = T) +
  #stat_summary(fun=mean, colour="black", geom="point", 
   #            shape=18, size=3, show.legend=FALSE) +
  #stat_boxplot(geom = "errorbar", width = 0.1, alpha = 0.95) +
  geom_boxplot(notch=TRUE,alpha=0.99, width =0.1, color = "grey1", fill = NA, outlier.shape = NA)+
  scale_fill_manual(values = c("midnightblue","red","darkred")) +
  xlim(0.8, 1) +
  geom_vline(xintercept=0.95, linetype="dashed", color = alpha("black", 0.4), lwd = 1.5) +
  theme_default.lb.boxplot()
violins

ggsave(file = paste0("violins_redlist_srli_cropdiv_all_fullsamples_poster.png"),
       width = 12, height = 10, dpi = 800)



curves <- ggplot(data = combo, aes(x=log10(sp), y=mean, col = Dataset)) +
  geom_line(lwd = 1.6, alpha = 0.8)  +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  scale_color_manual(values = c("midnightblue","red","darkred")) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("black", 0.4), lwd = 1.5) +
  xlim(c(0,6)) +
  theme_default.lb() 
curves

ggsave(file = paste0("comp_curves_red_srli_wcvp.png"),curves,
       width = 12, height = 10, dpi = 600)

# curves 
colouring <- c( '#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', "#CC3311", "#00183E")

viri <- viridisLite::viridis(100)


samples_cumulative_rel$Continent <- as.factor(samples_cumulative_rel$Continent)
samples_cumulative_rel$Continent<-factor(samples_cumulative_rel$Continent, levels=c("AFRICA","ANTARCTIC","ASIA-TEMPERATE","ASIA-TROPICAL",
                                                                                    "AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC","SOUTHERN AMERICA","OVERALL"))
samples_cumulative_rel_red$Continent <- as.factor(samples_cumulative_rel_red$Continent)
samples_cumulative_rel_red$Continent<-factor(samples_cumulative_rel_red$Continent, levels=c("AFRICA","ANTARCTIC","ASIA-TEMPERATE","ASIA-TROPICAL",
                                                                                    "AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC","SOUTHERN AMERICA","OVERALL"))
samples_cumulative_rel_srli$Continent <- as.factor(samples_cumulative_rel_srli$Continent)
samples_cumulative_rel_srli$Continent<-factor(samples_cumulative_rel_srli$Continent, levels=c("AFRICA","ANTARCTIC","ASIA-TEMPERATE","ASIA-TROPICAL",
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
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.4), lwd = 1.5) +
  xlim(c(0,6)) +
  theme_default.lb() + 
  ggtitle("World Checklist of Vascular Plants") 


curves_red <- ggplot(data = samples_cumulative_rel_red, aes(x=log10(sp), y=mean, col = Continent)) +
  geom_line(lwd = 1, alpha = 0.5)  +
  geom_line(data = subset(samples_cumulative_rel_red, Continent == 'OVERALL'), 
            lwd = 2, color = "midnightblue", alpha = 0.7) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  scale_colour_manual(values = c(colouring)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.8), lwd = 1.5) +
  xlim(c(0,6)) +
  theme_default.lb() + 
  ggtitle("IUCN Red List") 


curves_srli <- ggplot(data = samples_cumulative_rel_srli, aes(x=log10(sp), y=mean, col = Continent)) +
  geom_line(lwd = 1, alpha = 0.5)  +
  geom_line(data = subset(samples_cumulative_rel_srli, Continent == 'OVERALL'), 
            lwd = 2, color = "midnightblue", alpha = 0.7) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  scale_colour_manual(values = c(colouring)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.8), lwd = 1.5) +
  xlim(c(0,6)) +
  theme_default.lb() + 
  ggtitle("Sampled Red List Index") 


g <- arrangeGrob(curves_srli,curves_red,curves, ncol = 1, nrow = 3)



ggarrange(curves_red, curves_srli, curves, ncol = 1, nrow = 3, common.legend = T, labels = "AUTO")
ggsave(file = paste0("curves_red_srli_wcvp.png"),
       width = 15, height = 30, dpi = 600)









#  First look ------------------------------------------------------------------




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


redlist_full <- fread("data/red/redlist_full_syn.csv")

redlist_sub <- redlist_full %>% 
  subset(., select=which(!duplicated(names(.)))) 


redlist_sub$value <- 1

#redlist_red <-  redlist_sub %>% 
 # separate(redlistCriteria, c("Crit1", "Crit2"), sep = " + ")
 
#redlist_red$redlistCriteria <-
 # gsub("\\(.*","",redlist_red$redlistCriteria) 
  
xx <- as.data.frame(table(redlist_sub$redlistCriteria))
names(xx) <- c("criteria", "criteria_count")


redlist_table <- xx %>% 
  filter(!criteria == "") %>% 
  arrange(desc(criteria_count)) %>% 
  slice_head(n = 10)

threatened <- redlist_sub %>% 
  filter(!redlistCriteria == "")

threat_criterion <- as.data.frame(table(threatened$redlistCriteria))


redlist_table <- xx %>% 
  filter(!criteria == "") %>% 
  arrange(desc(criteria_count)) %>% 
  slice_head(n = 20)

red_summary <- threatened %>%
  mutate(short_crit = case_when(
    grepl("A", redlistCriteria) & !grepl("B", redlistCriteria) ~ "Population",
    !grepl("A", redlistCriteria) & grepl("B", redlistCriteria) ~ "Range",
    grepl("A", redlistCriteria) & grepl("B", redlistCriteria) ~ "Population & Range",
    grepl("D", redlistCriteria) ~ "Range",
    grepl("C", redlistCriteria) ~ "Population",
    TRUE ~ "Other"
  ))

applied_crit <- as.data.frame(table(red_summary$short_crit))


ext <- red_summary %>% 
  filter(short_crit == "Other")

redlist_table <- xx %>% 
  filter(!criteria == "") %>% 
  arrange(desc(criteria_count)) %>% 
  slice_head(n = 20)


sum(redlist_table$criteria_count)

redlist_comp <- redlist_sub %>% 
  dplyr::select(yearPublished, value)



ggplot(data = redlist_table, aes(x=as.factor(criteria), y = value)) +
  geom_bar(stat = "identity", fill = "black", col = "black") + 
  labs(x = "Criteria applied", y = "Count") +
  theme_bw(base_size = 26)

ggsave(file = paste0("redlist_criteria_freq.png"),
       width = 20, height = 5, dpi = 600)


redlist_table
time <- as.data.frame(table(redlist_sub$yearPublished))
names(time) <- c("year", "count")

time %>% 
  arrange(desc(year)) %>% 
  slice_head(n = 5) %>% 
  summarise(m = mean(count))

ggplot(data = redlist_comp, aes(x = value, group=value)) +
  geom_density(fill = "black", col = "black") + 
  labs(x = "Year published", y = "Number of assessments") +
  theme_bw(base_size = 26)


redlist_sub$yearPublished