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
library(rWCVP)
library(arrow)
library(plotrix)


# Defining functions -----------------------------------------------------------

error <- function(x) {
  qnorm(0.975)*sd(x)/sqrt(length(x))
} 


options(digits = 5)  
# Import data ------------------------------------------------------------------

#data_out2 <- read.table("output/full_samples_steps_increasing_niter100.txt")
richness_patterns <- read.csv("data/richness_patterns.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")

tdwg_3 <- wgsrpd3
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


#Get the path of filenames
list_data <- list.files("D:/samples/richness/gbif/parquet", full.names = TRUE)

# split between the part that comes before the numerics and the "1.img" etc.--adjust appropriately
split <- strsplit(list_data, ".*Samples_iteration_splm_") 
# strip the "1.img" etc such that only the numeric part is left
# turn the characters in numeric







split2 <- lapply(split, function(x) {
  temp <- gsub('[0-9]+.txt.gz','', x)
  extract <- unique(unlist(temp))
  pattern <- extract[!extract == ""]
  pattern_split <- str_split_1(pattern, "_")
  if (length(pattern_split) == 1) {
    pattern_split[2] <- "wcvp"
  }
  pattern_split
})

index_source <- unique(unlist(split2))[1:2]


list_data2 <- str_extract(list_data, '\\b\\w+$')

#Read them in a list
list_parquet <- list_data[which(list_data2 == "parquet")]
samplelist_parq <- lapply(list_parquet, read_parquet)

list_txt <- list_data[which(list_data2 == "gz")]
samplelist_txt <- lapply(list_txt, fread)

list_full <- c(samplelist_parq, samplelist_txt)


rm(list_parquet, list_txt)
gc()

#make big dataframe
data_out <- rbindlist(list_full) 

rm(list_full)
#filter(sp %in% c(1:3000)) 
#as.data.frame(do.call(rbind, samplelist)) 
setDT(data_out)
gc()

# Reframe ----------------------------------------------------------------------


samples_cumulative_rel <- rollup(data_out, list("mean_corgw" = mean(cor.sp_gwr_mean), 
                                                
                                                "sd_corgw_mean" = mean(cor.sp_gwr_sd, na.rm =T), 
                                                "se_corgw" = std.error(cor.sp_gwr_mean, na.rm =T),
                                                
                                                "mean_pseudor" = mean(pseudo.r.squared),
                                                "sd_pseudor" = sd(pseudo.r.squared, na.rm =T),
                                                "se_pseudor" = std.error(pseudo.r.squared, na.rm =T),
                                                
                                                "n" = length(cor.sp_gwr_mean),
                                                
                                                "lower_ci_corgw" = mean(cor.sp_gwr_mean) - error(cor.sp_gwr_mean), 
                                                "upper_ci_corgw" = mean(cor.sp_gwr_mean) + error(cor.sp_gwr_mean),
                                                
                                                "lower_ci_pseudor" = mean(pseudo.r.squared) - error(pseudo.r.squared), 
                                                "upper_ci_pseudor" = mean(pseudo.r.squared) + error(pseudo.r.squared),
                                                
                                                "resid_range" = mean(maxres - minres), 
                                                "mean_mae" = mean(mae), 
                                                "mean_sumres" = mean(sumres2)), 
                                 
                                 by = c("continent", "sp")) |> 
  drop_na() %>% 
  mutate(index = index_source[1], 
         source = index_source[2])


#Output analysis----------------------------------------------------------------

#samples_cumulative_rel <- fread("output/cumulative_rel/fullsamples_cumulative_rel.txt")

threshold <- samples_cumulative_rel %>% 
  filter(mean_pseudor > 0.95 & mean_corgw > 0.95) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 


threshold <- samples_cumulative_rel %>% 
  filter(mean_corgw > 0.95) %>% 
  group_by(continent) %>% 
  mutate(ratio = mean_pseudor/sp) %>% 
  arrange(desc(ratio)) %>%
  slice_min(sp) 
