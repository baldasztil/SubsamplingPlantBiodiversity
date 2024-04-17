library(tidyverse)
library(data.table)
library(ggsci)
library(patchwork)
library(hrbrthemes)
#library(extrafont)
library(RColorBrewer)
library(forcats)
library(viridis)
library(ggdist)
library(directlabels)
library(arrow)


#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/samples_aggregated/parquet", full.names = TRUE)


sample_data_con <- open_dataset(
  sources = "data/samples_aggregated/parquet/", 
  format = "parquet"
)

sample_data_con


samples_cumulative_rel_red <- sample_data_con %>% 
  filter(sp <= 150000) %>% 
  collect() %>% 
  setDT()

threshold <- sample_data_con %>% 
  filter(mean_corgw > 0.95 ) %>% 
  group_by(continent, index, source) %>% 
  collect() %>% 
  slice_min(sp) 


unique(samples_cumulative_rel_red)
summary(samples_cumulative_rel_red)

unique(samples_cumulative_rel_red$index)
#a <- colour("muted")
colouring <- c( '#CC6677', '#722280', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', "#CC3311", "#00183E")

colouring <- c( '#CC6677',  '#DDCC77', '#117733', '#88CCEE', '#882255')
viri <- viridisLite::viridis(100)

colouring[1:2]
samples_cumulative_rel$continent <- as.factor(samples_cumulative_rel$continent)


samples_cumulative_rel$continent<-factor(samples_cumulative_rel$continent, levels=
                                           c("AFRICA","ANTARCTICA","ASIA-TEMPERATE","ASIA-TROPICAL"
                                             ,"AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC"
                                             ,"SOUTHERN AMERICA","GLOBAL"))



samples_cumulative_rel$continent2 <- samples_cumulative_rel$continent

unique(samples_cumulative_rel$source)
samples_notglobal <- samples_cumulative_rel %>% 
  filter(!continent == "GLOBAL" & source == "gbif.parquet") 


samples_global <- samples_cumulative_rel %>% 
  filter(continent == "GLOBAL" & source == "gbif.parquet") 

curves_notglobal <-  ggplot(data =samples_notglobal, aes(x=log10(sp), y=mean_corgw, 
                                                         lty = source)) +
  geom_line(data = samples_cumulative_rel %>%  
              filter(source == "wcvp") %>% 
              dplyr::select(-continent), 
            aes(group = continent2), col= "grey", lwd = 0.3, alpha = 0.5)  +
  geom_line(data = samples_notglobal, aes(color = continent), 
            lwd = 0.75,  alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_corgw, ymax = upper_ci_corgw), outline.type = "both", 
              linetype = 0, alpha = 0.25) +
  #scale_colour_manual(values = c(colouring)) +
  scale_color_frontiers() +
  #scale_color_npg() +
  #scale_color_brewer(palette = "Set3") +  
  facet_wrap(~continent, ncol = 3, nrow = 3) +
  theme_ipsum(axis_title_size = 11.5, base_family = "Helvetica") +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  )


curves_global <-  ggplot(data = samples_global, aes(x=log10(sp), y=mean_corgw)) +
  geom_line(data = samples_cumulative_rel %>%  
              filter(source == "wcvp") %>% 
              dplyr::select(-continent), 
            aes(group = continent2), col= "grey", lwd = 0.75, alpha = 0.5)  +
  geom_line(data = samples_global, aes(color = continent), 
            lwd = 1.5, color = "black", alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 1) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) 
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_corgw, ymax = upper_ci_corgw), outline.type = "both", 
              linetype = 0, alpha = 0.25) +
  scale_colour_manual(values = c(colouring)) +
  facet_wrap(~continent, ncol = 3, nrow = 3) +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  )

gc()
curves <- curves_global + curves_notglobal 
ggsave("curves_corgw_gbif.png", curves, width = 15, height = 10, dpi = 300)


  
library(ggrepel)

curves_global_split <-  ggplot(data = samples_global, aes(x=log10(sp), y=mean_corgw)) +
  geom_line(data = samples_global, aes(color = index, lty = source), 
            alpha = 0.8, lwd = 0.5) +
  geom_hline(yintercept=0.9, linetype="dashed", color = alpha("red", 0.6), lwd = 0.5) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_corgw, ymax = upper_ci_corgw, fill = index), 
              outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  geom_text_repel(data = samples_global %>%
                    filter(sp == last(sp)),
                  aes(label = index, 
                      group = continent, 
                      x = 5.542966 + 0.2, 
                      y = mean_corgw, 
                      color = index), 
                  force = 0.01,
                  force_pull = 1, 
                  verbose = T, 
                  box.padding = 0, 
                  xlim = c(5.542966 + 0.2, 5.542966 + 0.2)) + 
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~continent, ncol = 3, nrow = 4) +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  coord_cartesian(clip = "off") +
  theme(
    legend.position="bottomleft", 
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  ) 

curves_global_split

curves_notglobal_split <-  ggplot(data = samples_notglobal, aes(x=log10(sp), y=mean_corgw)) +
  geom_line(data = samples_notglobal, aes(color = index, lty = source), 
            alpha = 0.8, lwd = 0.5) +
  geom_hline(yintercept=0.9, linetype="dashed", color = alpha("red", 0.6), lwd = 0.5) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_corgw, ymax = upper_ci_corgw, fill = index), outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~continent, ncol = 3, nrow = 4) +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  theme(
    legend.position="bottomleft", 
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  )


curves_split <- curves_global_split + curves_notglobal_split 
ggsave("curves_split_corgw_gbif.png", curves_split, width = 17, height = 10, dpi = 300)


unique(samples_cumulative_rel$continent)

threshold <- samples_cumulative_rel %>% 
  filter(mean_corgw > 0.95 ) %>% 
  group_by(continent, index) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 
threshold
threshold$sp


whereat <- samples_cumulative_rel %>% 
  group_by(continent) %>% 
  filter(sp == 2264)

threshold <- samples_cumulative_rel %>% 
  filter(mean_corgw > 0.9) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 


threshold %>%
  ungroup() %>% 
  mutate(continent = fct_reorder(continent, sp))  %>% 
  ggplot() +
  geom_col(aes(sp, continent), fill = "black", width = 0.6) +
  theme_bw()



threshold_gbif <- samples_cumulative_rel %>% 
  filter(source == "gbif") %>% 
  filter(mean_corgw > 0.9) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 

threshold_wcvp <- samples_cumulative_rel %>% 
  filter(source == "wcvp") %>% 
  filter(mean_corgw > 0.9) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(continent)) %>%
  slice_min(sp) 

threshold_wcvp_splm <- samples_cumulative_rel %>% 
  filter(source == "wcvp") %>% 
  filter(mean_pseudor > 0.9) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(continent)) %>%
  slice_min(sp) 

threshold_wcvp$diff_sp <- (threshold_wcvp$sp - threshold_wcvp_splm$sp) 
threshold_wcvp$diff_mae <- (threshold_wcvp$mean_mae - threshold_wcvp_splm$mean_mae)


threshold_both <- samples_cumulative_rel %>% 
  filter(source == "gbif") %>% 
  right_join(threshold_wcvp, by = c("sp", "continent"))

threshold_both_sp <- threshold_wcvp %>% 
  left_join(threshold_gbif, by = c("continent"))




dumbell <- ggplot() +
  geom_segment(data = threshold_both, aes(x=continent, xend=continent, y=mean_corgw.x, yend=mean_corgw.y), 
               color="grey83", linewidth = 1, show.legend = F) +
  geom_point(data = threshold_both, aes(x=continent, y=mean_corgw.x), color="darkred", size=3, show.legend = T) +
  geom_point(data = threshold_both,aes(x=continent, y=mean_corgw.y), color="black", size=3, show.legend = T) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "topleft", 
        axis.text.x = element_text(size = 14)) +
  xlab("") +
  ylab("Correlation coefficient") 
  #ylim(c(0,1))

dumbell
  

dumbell_sp <- ggplot() +
  geom_segment(data = threshold_both_sp, aes(x=continent, xend=continent, y=sp.y, yend=sp.x), 
               color="grey83", linewidth = 1, show.legend = F) +
  geom_point(data = threshold_both_sp, aes(x=continent, y=sp.y), color="darkred", size=3, show.legend = T) +
  geom_point(data = threshold_both_sp,aes(x=continent, y=sp.x), color="black", size=3, show.legend = T) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "topleft", 
        axis.text.x = element_text(size = 14)) +
  xlab("") +
  ylab("Species") 
#ylim(c(0,1))

dumbell_sp

samples_full <- rbind(samples_cumulative_rel, samples_cumulative_rel_gbif)
samples_full_glob <- samples_full %>% 
  filter(continent == "GLOBAL") 


samples_notglobal_start <- samples_notglobal %>% 
  filter(sp %in% 1:max(threshold_both_sp$sp.x))

violin_not_global <-  ggplot(data = samples_notglobal_start, aes(x=log10(sp), y=mean_corgw, fill = source, lty = index)) +
  geom_violin() +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 0.5) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  scale_colour_manual(values = c("darkred", "black")) +
  scale_fill_manual(values = c("darkred", "black")) +
  facet_wrap(~continent, ncol = 3, nrow = 4) +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  theme(
    legend.position="bottomleft", 
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  )
violin_not_global


samples_global_start <- samples_global %>% 
  filter(sp %in% 1:max(threshold_both_sp$sp.x))

violin_global <-  ggplot(data = samples_global_start, aes(x=log10(sp), y=mean_corgw, fill = source, lty = index)) +
  geom_violin() +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 0.5) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  scale_colour_manual(values = c("darkred", "black")) +
  scale_fill_manual(values = c("darkred", "black")) +
  facet_wrap(~continent, ncol = 3, nrow = 4) +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  theme(
    legend.position="bottomleft", 
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  )
violin_global

violins <- violin_global + violin_not_global 

ggsave("violins_global.png", violins, width = 15, height = 10, dpi = 300)


samples_global_hit95 <- samples_global %>% 
  group_by(index) %>% 
  filter(mean_corgw > 0.95) %>% 
  slice_min(sp)

curves_notglobal_split <-  ggplot(data = samples_global_hit95, aes(y= log10(sp), x = 1)) +
  geom_text(aes(label = index))  +
  theme_void() +
curves_notglobal_split
