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



rewrite_parquet<- function(file_path) {
  x <- fread(file = file_path)
  write_parquet(x, paste0(file_path, ".parquet"))
  gc()
}
#Get the path of filenames
list_data <- list.files("D:/samples/functional/wcvp", full.names = TRUE)

# split between the part that comes before the numerics and the "1.img" etc.--adjust appropriately
split <- strsplit(list_data, ".*Samples_iteration_splm_") 
# strip the "1.img" etc such that only the numeric part is left
# turn the characters in numeric  
list_data2 <- str_extract(list_data, '\\b\\w+$')


list_txt <- list_data[which(list_data2 == "gz")]
lapply(list_txt, rewrite_parquet)

gc()


