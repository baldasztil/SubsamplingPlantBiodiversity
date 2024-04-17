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
richness_patterns <- read.csv("output/overall_richness.csv")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")

tdwg_3 <- wgsrpd3
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)


#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/samples/tax", full.names = TRUE)

# split between the part that comes before the numerics and the "1.img" etc.--adjust appropriately
split <- strsplit(list_data, ".*Samples_iteration_splm_") 
# strip the "1.img" etc such that only the numeric part is left
# turn the characters in numeric



split2 <- lapply(split, function(x) {
  temp <- gsub('[0-9]+.parquet','', x)
  extract <- unique(unlist(temp))
  pattern <- extract[!extract == ""]
  pattern_split <- str_split_1(pattern, "_")
  if (length(pattern_split) == 1) {
    pattern_split[2] <- "wcvp"
  }
  pattern_split
})

index_source <- unique(unlist(split2))

#Read them in a list
samplelist <- lapply(list_data, read_parquet)

#make big dataframe
data_out <- rbindlist(samplelist) 
  #filter(sp %in% c(1:3000)) 
  #as.data.frame(do.call(rbind, samplelist)) 
gc()
setDT(data_out)
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
  
fwrite(samples_cumulative_rel, paste0("./data/samples_aggregated/",index_source[1],"_", index_source[2],"_",
                                      "fullsamples_splm_agg.txt"))

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

