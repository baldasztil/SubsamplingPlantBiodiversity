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






rewrite_parquet<- function(x) {
  aa <- unique(x$index)
  bb <- unique(x$source)
  write_parquet(x, paste0("samples_agg_", aa, "_", bb, ".parquet"))
  gc()
}


#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/samples_aggregated", full.names = TRUE)

#Read them in a list
samples_cumulative_rel <- setDT(rbindlist(lapply(list_data, fread))) %>% 
  mutate(source = ifelse(source == "gbif.parquet", "gbif", source),
         index = ifelse(index == "rich1158059.parquet", "rich", index))


setDT(samples_cumulative_rel)

samples_cumulative_rel[, fwrite(.SD, paste0("samples_agg_",  "_", ".parquet")), by = c(index, source)]


xx <- sample_data_con %>% 
  group_by(source, index)  %>% 
  summarise() %>% 
  collect() 
  



aa <- unique(samples_cumulative_rel$index)
bb <- unique(samples_cumulative_rel$source)


list_samples <- split(samples_cumulative_rel, by = c("index", "source")) %>% 
  lapply(rewrite_parquet)
  
  
  
  split(index, source) %>% 
  rewrite_parquet()

