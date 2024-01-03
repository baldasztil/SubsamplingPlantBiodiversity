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

wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
#-------------------------------------------------------------------------------

wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

wcvp_sp <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(!taxon_rank %in% c("Genus"))

wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) 

#-------------------------------------------------------------------------------

dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% 
  filter(!area == "")

dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))

plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")

plants_full_extinct_norange <- wcvp_accepted %>% 
  filter(!plant_name_id %in% plants_full$plant_name_id)

#-------------------------------------------------------------------------------

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

richness_patterns <- rbind(richness_patterns1,richness_patterns2, richness_patterns3)

#-------------------------------------------------------------------------------


# Adjust data ------------------------------------------------------------------

#table of number of unique species per botanical unit here at different levels 


# Analysis ---------------------------------------------------------------------
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

#plantlist_names <- plants_full %>%  
 # select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name) %>%
  #sample_n(50000)

tdwg_codes <- as.data.frame(tdwg_3) %>% 
  select(LEVEL3_COD,LEVEL2_COD,LEVEL1_COD)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

start.time <- Sys.time()
######
cumulative_namelist <- data.frame()
rich_rel_cumulative <- data.frame()
samplelist <- list()
for (x in 1:100) {
    repeat{
      
      # randomly sampling the data for species ids 
      plantlist_names_left <- plantlist_names %>% 
        filter(!plant_name_id %in%  cumulative_namelist$plant_name_id)
      if (nrow(plantlist_names_left) > 500) {
        species_sample <- sample_n(plantlist_names_left, 500)
      } else {
        species_sample <- plantlist_names_left
      }
      cumulative_namelist <- rbind(cumulative_namelist, species_sample)
      species <- cumulative_namelist
      # caluclating sample species richness
      dist <- plantlist_dist %>% 
        filter(plant_name_id %in% species$plant_name_id)
      
      # make richness patterns with brus and group by continent/region instead of ID to look at within patterns
      sample_rich_bru <- dist %>% 
        group_by(area_code_l3) %>% 
        summarise(richness_sample = n_distinct(plant_name_id)) %>% 
        rename(LEVEL_COD = area_code_l3) 
      # measuring correlation
      rich_rel <- rich_overall_bru %>% 
        left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
        replace(is.na(.), 0)
      
      riche_rel_bru <- rich_rel %>% 
        summarise(cor.sp = cor(richness_sample, richness, method="spearman"),
                  cor.pe = cor(richness_sample, richness, method="pearson")) %>% 
        mutate(sp = nrow(cumulative_namelist),
               id = "overall")
      rich_rel_cumulative <- rbind(rich_rel_cumulative,riche_rel_bru)
          for (i in 1:9) {
            continent_pattern <- rich_rel %>% 
              filter(LEVEL1_COD == i)
            rich_rel_con <- continent_pattern %>% 
              summarise(cor.sp = cor(richness_sample, richness, method="spearman"),
                        cor.pe = cor(richness_sample, richness, method="pearson")) %>% 
              mutate(sp = nrow(cumulative_namelist),
                     id = as.character(i))
            if (rich_rel_con$id > 0) {
              rich_rel_cumulative <- rbind(rich_rel_cumulative,rich_rel_con)
            }
            if (i == 9) {
              break
            }
          }
      if (nrow(cumulative_namelist) == nrow(plantlist_names)){
        break
      }
    if (nrow(cumulative_namelist) %in% seq(0,nrow(plantlist_names),1000)){
    print(paste0("There are ", nrow(cumulative_namelist) ," species in the subsample")) 
    }
  }
    samplelist[[x]] <- rich_rel_cumulative 
    cumulative_namelist <- data.frame()
    rich_rel_cumulative <- data.frame()
    print(paste0("This is iteration ", x))
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
}


data_out <- as.data.frame(do.call(rbind, samplelist))

rich_average <- data_out %>%
  group_by(id,sp) %>% 
  summarise(mean = mean(cor.sp, na.rm = TRUE),
            sd = sd(cor.sp, na.rm = TRUE),
            n.sample = n()) %>%
  mutate(se = sd / sqrt(n.sample),
         lower.ci = mean - qt(1 - (0.05 / 2), n.sample - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.sample - 1) * se)



tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


rich_rel_cumulative_names <- rich_average %>% 
  left_join(tdwg_1_names, by = c("id" = "LEVEL1_COD")) %>% 
  mutate(Continent = LEVEL1_NAM,
         Continent = ifelse(is.na(Continent), "OVERALL", Continent))

rich_rel_cumulative_names %>% 
  ggplot(aes(x=sp, y=mean, col = Continent)) +
  geom_line(lwd = 0.7)  +
  xlab("Sampe size (n species)") +
  ylab("Correlation coefficient") +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), outline.type = "both", linetype = 0, alpha = 0.1, ) +
  theme_classic()

ggsave(paste0("curves_", "n",max(rich_rel_cumulative_names$sp),"_steps",diff(rich_rel_cumulative_names$sp)[1],".jpeg"),  width = 12, height = 8, dpi = 600)

?geom_ribbon
repeat{
  
  # randomly sampling the data for species ids 
  plantlist_names_left <- plantlist_names %>% 
    filter(!plant_name_id %in%  cumulative_namelist$plant_name_id)
  species_sample <- sample_n(plantlist_names_left, 1) 
  cumulative_namelist <- rbind(cumulative_namelist, species_sample)
  species <- cumulative_namelist
  # caluclating sample species richness
  dist <- plantlist_dist %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # make richness patterns with brus and group by continent/region instead of ID to look at within patterns
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3) %>% 
    left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))
  # measuring correlation
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0)
  
  riche_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor(richness_sample, richness, method="spearman"),
              cor.pe = cor(richness_sample, richness, method="pearson")) %>% 
    mutate(sp = nrow(cumulative_namelist))
           
  rich_rel_cumulative <- rbind(rich_rel_cumulative,riche_rel_bru)
  # data entry
  # can try an make it into long dataframe with scale (con,reg,biu) as group 
  # variable and then calculate the means 
  
  if (nrow(cumulative_namelist) == nrow(plantlist_names)){
    break
  }
  print(paste0("There are ", nrow(cumulative_namelist) ," species in the subsample"))
}






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



repeat{

    # randomly sampling the data for species ids 
    plantlist_names_left <- plantlist_names %>% 
      filter(!plant_name_id %in%  cumulative_namelist$plant_name_id)
    
    species_sample <- sample_n(plantlist_names_left, 1) 
    cumulative_namelist <- rbind(cumulative_namelist, species_sample)
    species <- cumulative_namelist
    # caluclating sample species richness
    dist <- plantlist_dist %>% 
      filter(plant_name_id %in% species$plant_name_id)
    
    # make richness patterns with brus and group by continent/region instead of ID to look at within patterns
    sample_rich_con <- dist %>% 
      group_by(continent_code_l1) %>% 
      summarise(richness_sample = n_distinct(plant_name_id))%>% 
      rename(LEVEL_COD = continent_code_l1) 
    
    sample_rich_reg <- dist %>% 
      group_by(region_code_l2) %>% 
      summarise(richness_sample = n_distinct(plant_name_id)) %>% 
      rename(LEVEL_COD = region_code_l2) 
    
    sample_rich_bru <- dist %>% 
      group_by(area_code_l3) %>% 
      summarise(richness_sample = n_distinct(plant_name_id)) %>% 
      rename(LEVEL_COD = area_code_l3) 
    
    sample_rich_pat <- rbind(sample_rich_bru,sample_rich_reg, 
                               sample_rich_con)
    # measuring correlation
    rich_rel <- rich_overall %>% 
      left_join(sample_rich_pat, by = "LEVEL_COD") %>% 
      group_by(ID) %>% 
      summarise(cor.sp = as.numeric(cor.test(richness, richness_sample, method="spearman")[[4]]),
                cor.pear = as.numeric(cor.test(richness, richness_sample, method="pearson")[[4]])) %>% 
      mutate(sp = nrow(cumulative_namelist),
             ID = ID)
    
    rich_rel_cumulative <- rbind(rich_rel_cumulative,rich_rel)
    # data entry
    # can try an make it into long dataframe with scale (con,reg,biu) as group 
    # variable and then calculate the means 
    
    if (nrow(cumulative_namelist) == nrow(plantlist_names)/10){
      break
    }
  print(paste0("There are ", nrow(cumulative_namelist) ," species in the subsample"))
}



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

# useful tips ------------------------------------------------------------------
# go to Edit/Folding/Collapse all to collapse all sections
# save a pdf and a png of your files
# objects names -> x_y 
# function names -> x.y
# files -> x_y_z_date.R
# for inline comments -> <space><space>#<space>comment
# library("formatR") -> function tidy_source() & tidy_dir() to clean old scripts
# to replace -> gsub(".", "_", names(dataframe), fixed = TRUE)
# change case -> tolower(names(dataframe))
