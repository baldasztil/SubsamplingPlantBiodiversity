# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------

# Libraries --------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(kewr)
library(purrr)
library(progress)
library(stringr)
library(kewr)
library(taxize)
library(expowo)
library(data.table)
library(skimr)
library(sf)
library(ggplot2)
library(ggpubr)
# Defining functions -----------------------------------------------------------

# Working directory ------------------------------------------------------------
# Import data ------------------------------------------------------------------
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
richness_patterns <- fread("data/richness_patterns.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")

# manipulate data --------------------------------------------------------------
# extracting names from shapefile  
tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

# making it a character for later
tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

tdwg_codes <- as.data.frame(tdwg_3) %>% 
  select(LEVEL3_COD,LEVEL2_COD,LEVEL1_COD)

tdwg_codes$LEVEL1_COD <- as.character(tdwg_codes$LEVEL1_COD)

index_endemic <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n == 1) # 56%

index_widespread <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) # 43.9%

# endemic ----------------------------------------------------------------------
dist_endemic <- dist_native %>% 
  filter(plant_name_id %in% index_endemic$plant_name_id)

plants_endemic <- plants_full %>% 
  filter(plant_name_id %in% index_endemic$plant_name_id)

plantlist_dist <- dist_endemic %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_endemic %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

######
start.time <- Sys.time()
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
      # too include last sample (smaller than 500)
      if (nrow(cumulative_namelist) < 10) {
        species_sample <- sample_n(plantlist_names_left, 1)
        }
        if (nrow(cumulative_namelist) >= 10 & nrow(cumulative_namelist) < 100) {
          species_sample <- sample_n(plantlist_names_left, 10)
          }
          if (nrow(cumulative_namelist) >= 100 & nrow(cumulative_namelist) < 1000){
            species_sample <- sample_n(plantlist_names_left, 100)
           } 
           if (nrow(cumulative_namelist) >= 1000 & nrow(plantlist_names_left) >= 1000){
             species_sample <- sample_n(plantlist_names_left, 1000)
              } 
              if (nrow(plantlist_names_left) < 1000){
                species_sample <- plantlist_names_left
              }
                if (nrow(plantlist_names_left) == 0){
                  break
      }
      cumulative_namelist <- rbind(cumulative_namelist, species_sample)
      species <- cumulative_namelist
      
      # caluclating sample species richness
      dist <- dist_native %>% 
        filter(plant_name_id %in% species$plant_name_id)
      
      # richness patterns across brus
      sample_rich_bru <- dist %>% 
        group_by(area_code_l3) %>% 
        summarise(richness_sample = n_distinct(plant_name_id)) %>% 
        rename(LEVEL_COD = area_code_l3)
      
      rich_rel <- rich_overall_bru %>% 
        left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
        replace(is.na(.), 0) %>% 
        mutate(sp = nrow(cumulative_namelist))
      
      rich_cumulative <- rbind(rich_cumulative, rich_rel)
      # measuring correlation
      
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
        
        # bind all
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
        break
      }
    if (nrow(cumulative_namelist) %in% seq(0,nrow(plantlist_names),10000)){
    print(paste0("There are ", nrow(cumulative_namelist) ,
                 " species in the subsample")) 
    }
  }
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

data_out <- as.data.frame(do.call(rbind, samplelist))
data_out2 <- as.data.frame(do.call(rbind, samples))
write.table(data_out, "output/full_samples_rel_steps_increasing_niter100_endemic.txt")
write.table(data_out2, "output/full_samples_steps_increasing_niter100_endemic.txt")

# Formatting output ------------------------------------------------------------

rich_average <- data_out %>%
  group_by(id,sp) %>% 
  mutate( cor.sp = ifelse(is.na(cor.sp), 0, cor.sp)) %>% 
  summarise(mean = mean(cor.sp, na.rm = TRUE),
            sd = sd(cor.sp, na.rm = TRUE),
            n.sample = n()) %>%
  mutate(se = sd / sqrt(n.sample),
         lower.ci = mean - qt(1 - (0.05 / 2), n.sample - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.sample - 1) * se)

rich_rel_cumulative_names <- rich_average %>% 
  left_join(tdwg_1_names, by = c("id" = "LEVEL1_COD")) %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "OVERALL", Continent))


rich_species <- data_out2
rich_species$LEVEL1_COD <- as.character(rich_species$LEVEL1_COD)

rich_cumulative_species <- rich_species %>% 
  left_join(tdwg_1_names, by = "LEVEL1_COD") %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "OVERALL", Continent))

rich_rel_cumulative_names %>% 
  summarise_all(max())

#Output analysis----------------------------------------------------------------

threshold <- rich_rel_cumulative_names %>% 
  mutate(threshold = case_when(mean > 0.9  ~ 'above', 
                               mean < 0.9 ~ 'below')) %>% 
  group_by(Continent) %>% 
  filter(min(sp) & threshold == "above") %>%
  slice_min(sp)



#Plotting-----------------------------------------------------------------------

rich_species_3000 <- rich_cumulative_species %>% 
  filter(sp == 3001) %>% 
  ggplot( aes(x =Continent, y =richness_sample, fill = Continent ))+
  geom_jitter(color="black", size=0.4, alpha=0.09) +
  geom_boxplot() +
  theme_bw()

rich_species_full <- rich_cumulative_species %>% 
  filter(sp == 3001) %>% 
  ggplot(aes(x =Continent, y =richness, fill = Continent))+
  geom_boxplot() +
  theme_bw()

ggarrange (rich_species_full,rich_species_3000, ncol = 2, nrow =1)
ggsave(paste0("speciesboxplot", "n",max(rich_rel_cumulative_names$sp),"_steps",
       diff(rich_rel_cumulative_names$sp)[1],"niter_",x,".jpeg"),  
       width = 25, height = 8, dpi = 600)

rich_rel_cumulative_names %>% 
  ggplot(aes(x=log10(sp), y=mean, col = Continent)) +
  geom_line(lwd = 0.7)  +
  xlab("Sample size species") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1, ) +
  theme_bw() +
  theme(panel.border = element_blank())

ggsave(paste0("endemic_curvesfirst1000", "n",max(rich_rel_cumulative_names$sp),
              "_increasingsteps",diff(rich_rel_cumulative_names$sp)[1],"niter_",x,".jpeg"),  
       width = 12, height = 8, dpi = 600)


rich_overall_con <- richness_patterns %>% 
  filter(ID =="con") %>% 
  left_join(tdwg_1_names, by =c("LEVEL_COD" = "LEVEL1_COD"))


tdwg_3_mod <- tdwg_3 %>% 
  filter(LEVEL3_COD %in% rich_rel$LEVEL_COD) %>% 
  left_join(rich_rel, by = c("LEVEL3_COD" = "LEVEL_COD")) 
pattern <- ggplot(data = tdwg_3_mod) +
  geom_sf(aes(fill = richness)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick3", na.value = "white")
sample <- ggplot(data = tdwg_3_mod) +
  geom_sf(aes(fill = richness_sample)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick3", na.value = "white")

ggarrange(pattern, sample, ncol =2, nrow =1)


tdwg_1_mod <- tdwg_1 %>% 
  filter(LEVEL1_COD %in%richness_patterns_con$LEVEL1_COD) %>% 
  left_join(richness_patterns_con, by = "LEVEL1_COD") 

con <- ggplot(data = tdwg_1_mod) +
  geom_sf(aes(fill = richness)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick1", na.value = "white")

tdwg_2_mod <- tdwg_2 %>% 
  filter(LEVEL2_COD %in%richness_patterns_reg$LEVEL2_COD) %>% 
  left_join(richness_patterns_reg, by = "LEVEL2_COD") 
reg <- ggplot(data = tdwg_2_mod) +
  geom_sf(aes(fill = richness)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick2", na.value = "white")

ggarrange(bru, reg, con,ncol =3, nrow =1)

# non-endemic ------------------------------------------------------------------
dist_widespread <- dist_native %>% 
  filter(plant_name_id %in% index_widespread$plant_name_id)

plants_widespread <- plants_full %>% 
  filter(plant_name_id %in% index_widespread$plant_name_id)


plantlist_dist <- dist_widespread %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_widespread %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

######
start.time <- Sys.time()
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
    # too include last sample (smaller than 500)
    if (nrow(cumulative_namelist) < 10) {
      species_sample <- sample_n(plantlist_names_left, 1)
    }
    if (nrow(cumulative_namelist) >= 10 & nrow(cumulative_namelist) < 100) {
      species_sample <- sample_n(plantlist_names_left, 10)
    }
    if (nrow(cumulative_namelist) >= 100 & nrow(cumulative_namelist) < 1000){
      species_sample <- sample_n(plantlist_names_left, 100)
    } 
    if (nrow(cumulative_namelist) >= 1000 & nrow(plantlist_names_left) >= 1000){
      species_sample <- sample_n(plantlist_names_left, 1000)
    } 
    if (nrow(plantlist_names_left) < 1000){
      species_sample <- plantlist_names_left
    }
    if (nrow(plantlist_names_left) == 0){
      break
    }
    cumulative_namelist <- rbind(cumulative_namelist, species_sample)
    species <- cumulative_namelist
    
    # caluclating sample species richness
    dist <- dist_native %>% 
      filter(plant_name_id %in% species$plant_name_id)
    
    # richness patterns across brus
    sample_rich_bru <- dist %>% 
      group_by(area_code_l3) %>% 
      summarise(richness_sample = n_distinct(plant_name_id)) %>% 
      rename(LEVEL_COD = area_code_l3)
    
    rich_rel <- rich_overall_bru %>% 
      left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
      replace(is.na(.), 0) %>% 
      mutate(sp = nrow(cumulative_namelist))
    
    rich_cumulative <- rbind(rich_cumulative, rich_rel)
    # measuring correlation
    
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
    
    # bind all
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
      break
    }
    if (nrow(cumulative_namelist) %in% seq(0,nrow(plantlist_names),10000)){
      print(paste0("There are ", nrow(cumulative_namelist) ,
                   " species in the subsample")) 
    }
  }
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

data_out <- as.data.frame(do.call(rbind, samplelist))
data_out2 <- as.data.frame(do.call(rbind, samples))
write.table(data_out, "output/full_samples_rel_steps_increasing_niter100_nonendemic.txt")
write.table(data_out2, "output/full_samples_steps_increasing_niter100_nonendemic.txt")

rich_average <- data_out %>%
  group_by(id,sp) %>% 
  mutate( cor.sp = ifelse(is.na(cor.sp), 0, cor.sp)) %>% 
  summarise(mean = mean(cor.sp, na.rm = TRUE),
            sd = sd(cor.sp, na.rm = TRUE),
            n.sample = n()) %>%
  mutate(se = sd / sqrt(n.sample),
         lower.ci = mean - qt(1 - (0.05 / 2), n.sample - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.sample - 1) * se)

rich_rel_cumulative_names <- rich_average %>% 
  left_join(tdwg_1_names, by = c("id" = "LEVEL1_COD")) %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "OVERALL", Continent))


rich_species <- data_out2
rich_species$LEVEL1_COD <- as.character(rich_species$LEVEL1_COD)

rich_cumulative_species <- rich_species %>% 
  left_join(tdwg_1_names, by = "LEVEL1_COD") %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "OVERALL", Continent))

rich_rel_cumulative_names %>% 
  map_dbl(max(.))

rich_rel_cumulative_names %>% 
  ggplot(aes(x=log10(sp), y=mean, col = Continent)) +
  geom_line(lwd = 0.7)  +
  xlab("Sample size species (log 10)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", 
              linetype = 0, alpha = 0.1, ) +
  theme_bw() +
  theme(panel.border = element_blank())

ggsave(paste0("nonendemic_curves", "n",max(rich_rel_cumulative_names$sp),
              "_increasingsteps",diff(rich_rel_cumulative_names$sp)[1],"niter_",x,".jpeg"),  
       width = 12, height = 8, dpi = 600)



# go to Edit/Folding/Collapse all to collapse all sections
# save a pdf and a png of your files
# objects names -> x_y 
# function names -> x.y
# files -> x_y_z_date.R
# for inline comments -> <space><space>#<space>comment
# library("formatR") -> function tidy_source() & tidy_dir() to clean old scripts
# to replace -> gsub(".", "_", names(dataframe), fixed = TRUE)
# change case -> tolower(names(dataframe))
