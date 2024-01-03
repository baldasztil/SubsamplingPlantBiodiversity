library(pryr)
library(devtools)
library(lineprof)
library(profvis)

source("scripts/Crop_diversity/Crop_div_subsampling_parallelised.R")
prof <- profvis(parLapply(clust, seq(0,nrow(plantlist_names),1000), 
                                             subsampling.plants))
