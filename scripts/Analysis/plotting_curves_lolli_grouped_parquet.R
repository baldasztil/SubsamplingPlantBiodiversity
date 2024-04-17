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
library(ggrepel)



sample_data <- open_dataset(
  sources = "data/samples_aggregated/parquet/", 
  format = "parquet"
)

threshold_corgw <- sample_data %>% 
  filter(mean_corgw > 0.95) %>% 
  group_by(index, source, continent) %>% 
  collect() %>% 
  slice_min(sp)


samples_cumulative_rel <- sample_data %>% 
  filter(source == "wcvp") %>% 
  collect() %>% 
  setDT()
  

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


samples_notglobal <- samples_cumulative_rel %>% 
  filter(!continent == "GLOBAL") 


samples_global <- samples_cumulative_rel %>% 
  filter(continent == "GLOBAL")




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
  #facet_wrap(~index, ncol = 3, nrow = 4) +
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
ggsave("curves_split_corgw.png", curves_split, width = 17, height = 10, dpi = 300)






curves_global_mae <-  ggplot(data = samples_global, aes(x=log10(sp), y=mean_mae)) +
  geom_line(data = samples_cumulative_rel %>%  
              filter(source == "wcvp") %>% 
              dplyr::select(-continent), 
            aes(group = continent2), col= "grey", lwd = 0.75, alpha = 0.5)  +
  geom_line(data = samples_global, aes(color = continent), 
            lwd = 0.25, color = "darkorange", alpha = 0.8) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) 
  scale_x_continuous(name = "Sample size species (log10)", limits = c(0,6)) +
  ylab("Mean absolute error") +
  scale_colour_manual(values = c(colouring)) +
  facet_wrap(~index, ncol = 2, nrow = 2, scales = "free") +
  theme_ipsum(axis_title_size = 12, base_family = "Helvetica") +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    panel.grid = element_blank()
  )
#curves_global_mae

ggsave("curves_global_mae.png", curves_global_mae, width = 15, height = 5, dpi = 600)





curves_global_range <-  ggplot(data = samples_global, aes(x=log10(sp), y=resid_range)) +
  geom_line(data = samples_cumulative_rel %>%  
              filter(source == "wcvp") %>% 
              dplyr::select(-continent), 
            aes(group = continent2), col= "grey", lwd = 0.75, alpha = 0.5)  +
  geom_line(data = samples_global, aes(color = continent), 
            lwd = 0.25, color = "darkorange", alpha = 0.8) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) 
  scale_x_continuous(name = "Sample size species (log10)", limits = c(0,6)) +
  ylab("Residual range") +
  scale_colour_manual(values = c(colouring)) +
  facet_wrap(~index, ncol = 2, nrow = 2, scales = "free") +
  theme_bw( base_family = "Helvetica") +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
  )
#curves_global_mae

ggsave("curves_global_resid_range.png", curves_global_range, width = 15, height = 5, dpi = 600)


ggplot()
threshold_corgw



xx <- samples_global %>% 
  filter(index == "fun") %>% 
  filter(continent == "GLOBAL")


threshold_corgw2 <- sample_data %>% 
  filter(mean_corgw > 0.9) %>% 
  group_by(index, continent) %>% 
  collect() %>% 
  slice_min(sp)

sums <- samples_global %>% 
  filter(index == "rich") %>% 
  filter(resid_range < 349113/ 100 ) %>% 
  slice_min(sp, n = 1 )


threshold_corgw2 <- samples_global %>% 
  filter(mean_corgw > 0.95) %>% 
  group_by(index, continent) %>% 
  collect() %>% 
  slice_min(sp)

threshold_corgw2 <- samples_global %>% 
  filter(mae > 0.95) %>% 
  group_by(index, continent) %>% 
  collect() %>% 
  slice_min(sp)

unique_combinations <- unique(threshold_corgw2 %>% select(continent, sp, index))


threshold_global_wide <- sample_data %>% 
  filter(continent %in% unique_combinations$continent &
         sp %in% unique_combinations$sp &
         index %in% unique_combinations$index) %>% 
    collect() %>% 
  filter(
    paste(continent, sp, index) %in% paste(unique_combinations$continent, 
                                           unique_combinations$sp, 
                                           unique_combinations$index)) %>% 
  dplyr::select(mean_corgw, continent, sp, source, index) %>% 
  pivot_wider(
    names_from = source,
    names_glue = "{source}_{.value}",
    values_from = c(mean_corgw)
  ) %>% 
  arrange(desc(gbif_mean_corgw))

threshold_global_wide_red <-threshold_global_wide %>% 
  filter(continent == "GLOBAL")
  

dumbell <- ggplot() +
  geom_segment(data = threshold_global_wide, aes(x=continent, xend=continent, y=wcvp_mean_corgw, yend=gbif_mean_corgw), 
               color="grey83", linewidth = 1, show.legend = F) +
  geom_point(data = threshold_global_wide, aes(x=continent, y=gbif_mean_corgw), color="darkred", size=3, show.legend = T) +
  geom_point(data = threshold_global_wide,aes(x=continent, y=wcvp_mean_corgw), color="black", size=3, show.legend = T) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "topleft", 
        axis.text.x = element_text(size = 14)) +
  xlab("") +
  ylab("Correlation coefficient") +
  facet_wrap(~index)
  #ylim(c(0,1))

dumbell
  

threshold_global_wide_sp <- sample_data %>% 
  filter(mean_corgw > 0.9) %>% 
  group_by(index, continent) %>% 
  collect() %>% 
  slice_min(sp) %>% 
  dplyr::select(mean_corgw, continent, sp, index) %>% 
  mutate(log10_sp = log10(sp)) %>% 
  pivot_wider(
    names_from = source,
    names_glue = "{source}_{.value}",
    values_from = c(mean_corgw, sp, log10_sp)
  ) 

names(threshold_global_wide)


dumbell_sp <- ggplot() +
  geom_segment(data = threshold_global_wide_sp, aes(x=continent, xend=continent, y= wcvp_log10_sp, yend=gbif_log10_sp), 
               color="grey83", linewidth = 1, show.legend = F) +
  geom_point(data = threshold_global_wide_sp, aes(x=continent, y=gbif_log10_sp), color="darkred", size=3, show.legend = T) +
  geom_point(data = threshold_global_wide_sp,aes(x=continent, y=wcvp_log10_sp), color="black", size=3, show.legend = T) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "topleft", 
        axis.text.x = element_text(size = 14)) +
  xlab("") +
  ylab("Species") +
  facet_wrap(~index, scales = "fixed")#ylim(c(0,1))
#ylim(c(0,1))

dumbell_sp

samples_full <- rbind(samples_cumulative_rel, samples_cumulative_rel_gbif)
samples_full_glob <- samples_full %>% 
  filter(continent == "GLOBAL") 


samples_notglobal_start <- samples_notglobal %>% 
  filter(sp %in% 1:max(threshold_global_wide_sp$wcvp_sp))

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
  filter(sp %in% 1000:max(threshold_global_wide_sp$wcvp_sp)) 

library(ggridges)

samples_global <- sample_data %>% 
  filter(continent == "GLOBAL") %>% 
  collect()

violin_global <-  ggplot(data = samples_global, aes(x = mean_mae, y = index,  fill = source)) +
 #geom_point() +
  geom_boxplot() +
  #geom_density_ridges(alpha = 0.7) +
  
  geom_hline(yintercept=0.95, linetype="dashed", color = alpha("red", 0.6), lwd = 0.5) +
  #geom_vline(xintercept= log10(701), linetype="dashed", color = alpha("red", 0.6), lwd = 0.8) +
  #geom_density() +
  xlab("Sample size species (log10)") +
  ylab("Correlation coefficient") +
  #scale_colour_manual(values = c( '#CC6677', '#722280', '#117733')) +
  scale_colour_manual(values = c("darkred", "black")) +
  scale_fill_manual(values = c("darkred", "black")) +
  facet_wrap(~index, ncol = 3, nrow = 4, scales = "free") +
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
  group_by(index, source) %>% 
  filter(mean_corgw > 0.9) %>% 
  slice_min(sp)


samples_maxgbif <- samples_global %>% 
  filter(source == "gbif") %>% 
  group_by(index) %>% 
  slice_max(mean_corgw)



curves_notglobal_split <-  ggplot(data = samples_global_hit95, aes(y= log10(sp), x = 1)) +
  geom_text(aes(label = index))  +
  theme_void() +
curves_notglobal_split





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

