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
# Defining functions -----------------------------------------------------------

# Working directory ------------------------------------------------------------
getwd()
# Import data ------------------------------------------------------------------

wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 

wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

wcvp_sp <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(!taxon_rank %in% c("Genus"))
  
colnames(wcvp_raw)
str(wcvp_raw)
wcvp_raw$accepted <- wcvp_raw$accepted_plant_name_id

# no additional life-forms in synonyms
life_form_check <- wcvp_placed %>%  
  group_by(accepted_plant_name_id) %>% 
  reframe(life_form = unique(lifeform_description)) %>% 
  filter(!life_form == "")

wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) 


dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 

dist_native <- dist_raw %>% 
  filter(introduced == 0)


length(unique(dist_native$plant_name_id))

table_check <- as.data.frame(table(wcvp_accepted_life_form$accepted_plant_name_id))

dist_patterns <- dist_raw %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))

length(unique(dist_patterns$plant_name_id))

wcvp_accepted_distfilt <- wcvp_accepted_dupfilt %>% 
  filter(plant_name_id %in% dist_raw$plant_name_id)


distraw_filt <- dist_raw %>% 
  filter(!plant_name_id %in% wcvp_accepted$plant_name_id)

plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")
table(plants_full$continent)                   

# need to adam plant name id but this is weird, az and wcs at then end?? 
#df <- wcvp_accepted %>% 
 # filter(grepl('wcs', plant_name_id)) 


#wcvp_accepted$plant_name_id <- as.numeric(gsub("[^0-9]", "", wcvp_accepted$plant_name_id))
#length(unique(wcvp_accepted$plant_name_id))


table_check <- as.data.frame(table(wcvp_accepted$plant_name_id))
table_check_dup <- table_check %>% 
  filter(Freq >1)

wcvp_accepted_dup <- wcvp_accepted %>% 
  filter(plant_name_id %in% table_check_dup$Var1)%>% 
  filter(grepl('az', accepted_plant_name_id))

#dist_dup <- dist_patterns %>%  
 # filter (plant_name_id %in% c(329046,329055, 329058))



wcvp_accepted_dupfilt <- wcvp_accepted %>% 
  filter(!accepted_plant_name_id %in% wcvp_accepted_dup$accepted_plant_name_id)

length(unique(wcvp_accepted$plant_name_id))

wcvp_accepted_dup <- wcvp_accepted_dupfilt %>% 
  filter(grepl('az', accepted_plant_name_id))


wcvp_accepted_distfilt <- wcvp_accepted_dupfilt %>% 
  filter(!plant_name_id %in% dist_raw$plant_name_id)


distraw_filt <- dist_raw %>% 
  filter(!plant_name_id %in% wcvp_accepted_distfilt$plant_name_id)

plants_full <- left_(wcvp_accepted, dist_patterns, by = "plant_name_id")

non_present <- filter(plants_full, )

library(splitstackshape)
aa <- cSplit(distribution_patterns, "continent", ",")
aa <- cSplit(distribution_patterns, "continent", ",")

?cSplit_e
# Adjust data ------------------------------------------------------------------

#table of number of unique species per botanical unit here at different levels 

global_speciesrichness_con <- as.data.frame(table(plantlist_raw$con))
global_speciesrichness_reg <- as.data.frame(table(plantlist_raw$reg))
global_speciesrichness_biu <- as.data.frame(table(plantlist_raw$biu))

# Analysis ---------------------------------------------------------------------

# R objects setup for null model 
list <- list()
list2 <- list()
#
frame <- data.frame(matrix(ncol = 1, nrow = floor(length(plantlist)/500)))
names(frame) <- c("corr.test")

cumulative_namelist <- data.frame()


# make progress bar for models 
n_iter <- length(frame) # Number of iterations of the loop
# Initializes the progress bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")  

# for same results 
for (i in 1:100) {
  for (x in 1:n_iter){
    setTxtProgressBar(pb,x)
    # randomly sampling the data for species names 
    species_sample <- as.data.frame(sample(plantlist, 
                                           length(plantlist$names),500))
    cumulative_namelist <- rbind(cumulative_namelist, species_sample)
    species <- cumulative_namelist
    species(names) <- names(plantlist)
    # caluclating sample species richness
    sample_speciesrichness_con <- as.data.frame(table(species$con))
    sample_speciesrichness_reg <- as.data.frame(table(species$reg))
    sample_speciesrichness_biu <- as.data.frame(table(species$biu))
    # measuring correlation 
    cor_con <- cor.test(global_speciesrichness_con$freq, 
                        sample_speciesrichness_con$freq, family="spearman")
    cor_reg <- cor.test(global_speciesrichness_reg$freq, 
                        sample_speciesrichness_reg$freq, family="spearman")
    cor_biu <- cor.test(global_speciesrichness_biu$freq, 
                        sample_speciesrichness_biu$freq, family="spearman")
    # data entry
    # can try an make it into long dataframe with scale (con,reg,biu) as group 
    # variable and then calculate the means 
    frame$cor[x] <- cor_con$estimate[1]
    frame$scale[x] <- "con"
    
    frame$cor[x+1] <- cor_reg$estimate[1]
    frame$scale[x+1] <- "reg"
    
    frame$cor_biu[x+2] <- cor_biu$estimate[1]
    frame$scale[x+2] <- "biu"
  }
stat_summary$mean[i] <- frame %>% 
  group_by(scale) %>% 
  summary(mean())

}

  setTxtProgressBar(pb, x)
  

s_spear <- sd(frame$corr.test)
ci_spear <- CI(frame$corr.test) 
b <- c(ci_spear,s_spear)

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
