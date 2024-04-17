library(tidyverse)
library(arrow)
library(feather)
gbif_parq <- read_parquet("data/Samples_iteration_splm_rich1034079.parquet", as_data_frame = T)  %>% 
  mutate_if(is.numeric,
            round,
            digits = 3) %>% 
  dplyr::select(-n, -p, -npar, -value, -AIC, - range)

file <- paste0("data/Samples_iteration_splm_rich_gbif",
               sample(1:10000000, 1, replace=TRUE),".parquet")
write_parquet(gbif_parq, file)
