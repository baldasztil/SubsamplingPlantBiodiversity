#install.packages("tidyverse", 
 #                  lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", dependencies = T, repos = "https://cran.r-project.org")
  #install.packages("data.table", 
   #                lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", dependencies = T, repos = "https://cran.r-project.org")
install.packages("spmodel", 
                   lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", dependencies = T, repos = "https://cran.r-project.org")
  install.packages("ape", 
                   lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", dependencies = T, repos = "https://cran.r-project.org")
  install.packages("phyloregion", 
                   lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", dependencies = T, repos = "https://cran.r-project.org")
  install.packages("GWmodel",
                   lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", dependencies = T, repos = "https://cran.r-project.org")
  install.packages("arrow", 
                   lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", 
                   dependencies = T, repos = "https://cran.r-project.org")

install.packages("feather", 
                   lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", 
                   dependencies = T, repos = "https://cran.r-project.org")

install.packages("DEoptim", 
                 lib = "/exports/eddie3_homes_local/s2117440/R/x86_64-pc-linux-gnu-library/4.3", 
                 dependencies = T, repos = "https://cran.r-project.org")
  
  
  #conda config --add pkgs_dirs /exports/eddie3_homes_local/s2117440/.conda/pkgs/
  #conda install arrow 
  
  #conda create --name R_clonal r-tidyverse r-sf r-terra r-arrow 
  
  
#rsync -r  /exports/eddie3_homes_local/s2117440/ s2117440@baltic10.geos.ed.ac.uk:/home/s2117440/data

#mv /exports/eddie3_homes_local/s2117440/array-jobs/first-array-job.sh /exports/eddie3_homes_local/s2117440/

# $ qsub -pe sharedmem 8 -R y jobscript.sh # parallel job 




    

aa <- feather::read_feather("data/checkpoint_phylo_10.feather")


install.packages("BIEN", dependencies = T, repos = "https://cran.r-project.org")