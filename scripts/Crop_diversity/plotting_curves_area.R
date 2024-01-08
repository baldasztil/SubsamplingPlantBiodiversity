library(tidyverse)
library(data.table)
library(ggsci)
library(patchwork)
library(hrbrthemes)
#library(extrafont)
library(RColorBrewer)

samples_cumulative_rel <- fread("output/cumulative_rel/fullsamples_cumulative_rel.txt") %>% 
  filter(sp < 2500) %>% 
  mutate(Continent = ifelse(!Continent == "OVERALL", Continent, "GLOBAL"))

results <- read.csv("fullsamples_test.csv") 

brewer.pal(9, "Set3")
#a <- colour("muted")
colouring <- c( '#CC6677', '#722280', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', "#CC3311", "#00183E")



samples_cumulative_rel$Continent<-factor(samples_cumulative_rel$Continent, levels=
                                           c("AFRICA","ANTARCTIC","ASIA-TEMPERATE","ASIA-TROPICAL"
                                             ,"AUSTRALASIA","EUROPE","NORTHERN AMERICA","PACIFIC"
                                             ,"SOUTHERN AMERICA","GLOBAL"))



samples_cumulative_rel$Continent2 <- samples_cumulative_rel$Continent


samples_notglobal <- samples_cumulative_rel %>% 
  filter(!Continent == "GLOBAL")
  
(curves_notglobal <-  ggplot(data =samples_notglobal, aes(x=log10(sp), y=mean)) +
  geom_line(data = samples_cumulative_rel %>% dplyr::select(-Continent), 
            aes(group = Continent2), col= "grey", lwd = 0.3, alpha = 0.5)  +
  geom_line(data = samples_notglobal, aes(color = Continent), 
            lwd = 0.75,  alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.25) +
  #scale_colour_manual(values = c(colouring)) +
  scale_color_frontiers() +
  #scale_color_npg() +
  #scale_color_brewer(palette = "Set3") +  
  facet_wrap(~Continent, ncol = 3, nrow = 3) +
  theme_ipsum(axis_title_size = 11.5, base_family = "Helvetica") +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  ))

samples_global <- samples_cumulative_rel %>% 
  filter(Continent == "GLOBAL")

(curves_global <-  ggplot(data = samples_global, aes(x=log10(sp), y=mean)) +
  geom_line(data = samples_cumulative_rel %>% dplyr::select(-Continent), 
            aes(group = Continent2), col= "grey", lwd = 0.75, alpha = 0.5)  +
  geom_line(data = samples_global, aes(color = Continent), 
            lwd = 1.5, color = "black", alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 1) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +

  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.25) +
  scale_colour_manual(values = c(colouring)) +
  facet_wrap(~Continent, ncol = 3, nrow = 3) +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  ))


curves_global + curves_notglobal 




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
