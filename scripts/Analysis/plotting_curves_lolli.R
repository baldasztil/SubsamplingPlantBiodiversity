library(tidyverse)
library(data.table)
library(ggsci)
library(patchwork)
library(hrbrthemes)
#library(extrafont)
library(RColorBrewer)
library(forcats)


samples_cumulative_rel <- fread("fullsamples_splm_agg_wcvp_rich.txt")
#fwrite(samples_cumulative_rel, "fullsamples_splm_agg.txt", sep = ",")
#a <- colour("muted")
colouring <- c( '#CC6677', '#722280', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', "#CC3311", "#00183E")

viri <- viridisLite::viridis(100)


samples_cumulative_rel$continent <- as.factor(samples_cumulative_rel$continent)


samples_cumulative_rel$continent<-factor(samples_cumulative_rel$continent, levels=
                                           c("AFRICA","ANTARCTICA","ASIA-TEMPERATE","ASIA-TROPICAL"
                                             ,"AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC"
                                             ,"SOUTHERN AMERICA","GLOBAL"))



samples_cumulative_rel$continent2 <- samples_cumulative_rel$continent


samples_notglobal <- samples_cumulative_rel %>% 
  filter(!continent == "GLOBAL")

curves_notglobal <-  ggplot(data =samples_notglobal, aes(x=log10(sp), y=mean_pseudor)) +
  geom_line(data = samples_cumulative_rel %>% dplyr::select(-continent), 
            aes(group = continent2), col= "grey", lwd = 0.3, alpha = 0.5)  +
  geom_line(data = samples_notglobal, aes(color = continent), 
            lwd = 0.75,  alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_pseudor, ymax = upper_ci_pseudor), outline.type = "both", 
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

samples_global <- samples_cumulative_rel %>% 
  filter(continent == "GLOBAL")

curves_global <-  ggplot(data = samples_global, aes(x=log10(sp), y=mean_pseudor)) +
  geom_line(data = samples_cumulative_rel %>% dplyr::select(-continent), 
            aes(group = continent2), col= "grey", lwd = 0.75, alpha = 0.5)  +
  geom_line(data = samples_global, aes(color = continent), 
            lwd = 1.5, color = "black", alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 1) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_pseudor, ymax = upper_ci_pseudor), outline.type = "both", 
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
ggsave("curves_splm.png", curves, width = 15, height = 10, dpi = 300)


diff <- samples_cumulative_rel %>% 
  group_by(continent, sp) %>% 
  summarise (diff = mean_corgw - mean_pseudor)

diff %>% 
  group_by(continent) %>% 
  summarise( mean = mean(diff), 
             sd = sd(diff))

boxplot(diff ~ continent, col = colouring, data = diff)
abline(h = 0, lty = 2)

plot(mean_pseudor ~ mean_corgw,  data = samples_global)




library(viridis)
ggplot(samples_cumulative_rel, aes(x=log10(sp), y=mean, group=Continent, fill=Continent)) +
  geom_area(alpha = 0.4) +
  geom_line(data = samples_cumulative_rel, aes(x=log10(sp), y=mean, color = Continent), linewidth=0.75)+
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(legend.position="none") +
  #ggtitle("Popularity of American names in the previous 30 years") +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
    plot.title = element_text(size=14)) +
  facet_wrap(~Continent)


ggplot(samples_cumulative_rel, aes(x=log10(sp), y=mean, group=Continent, fill=Continent)) +
  geom_area() +
  scale_fill_manual(values = c(colouring)) +
  scale_fill_viridis_d() +
  theme(legend.position="none") +
  #ggtitle("Popularity of American names in the previous 30 years") +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
    plot.title = element_text(size=14)) 


unique(samples_cumulative_rel$continent)
threshold <- samples_cumulative_rel %>% 
  filter(mean_pseudor > 0.95) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 

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

samples_cumulative_rel_gbif <- fread("fullsamples_splm_agg_gbif_rich.txt")


threshold_gbif <- samples_cumulative_rel_gbif %>% 
  filter(mean_pseudor > 0.95) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 

threshold_both <- threshold %>% 
  left_join(samples_cumulative_rel_gbif, by = c("sp", "continent"))

threshold_both_sp <- threshold %>% 
  left_join(threshold_gbif, by = c("continent"))




dumbell <- ggplot() +
  geom_segment(data = threshold_both, aes(x=continent, xend=continent, y=mean_pseudor.y, yend=mean_pseudor.x), 
               color="grey83", linewidth = 1, show.legend = F) +
  geom_point(data = threshold_both, aes(x=continent, y=mean_pseudor.y), color="darkred", size=3, show.legend = T) +
  geom_point(data = threshold_both,aes(x=continent, y=mean_pseudor.x), color="black", size=3, show.legend = T) +
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

curves_global <-  ggplot(data = samples_full_glob, aes(x=log10(sp), y=mean_corgw)) +
  geom_line(data = samples_full_glob, aes(color = source, lty = index), 
           alpha = 0.7) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 0.5) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_corgw, ymax = upper_ci_corgw, fill = source), outline.type = "both", 
              linetype = 0, alpha = 0.1) +
  scale_colour_manual(values = c("darkred", "black")) +
  scale_fill_manual(values = c("darkred", "black")) +
  facet_wrap(~continent, ncol = 3, nrow = 4) +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  theme(
    legend.position="bottomleft", 
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  )
curves_global

library(ggdist)

violin_global <-  ggplot(data = samples_full_glob, aes(x=log10(sp), y=mean_corgw, fill = source, lty = index)) +
  stat_halfeye() +
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

ggsave("curves_global.png", curves_global, width = 15, height = 10, dpi = 300)
