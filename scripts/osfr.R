install.packages("osfr")
remotes::install_github("Between-the-Fjords/dataDownloader")
library(osfr)
library(dplyr)
library(dataDownloader)

OSF_PAT <- "Xedg2ur8YUNtgYufQiE74hHE4alDaMVKOU8Yy85zJwwiBPatCLlsTKs7YnbpMW6EWYUjcp"
osf_auth("Xedg2ur8YUNtgYufQiE74hHE4alDaMVKOU8Yy85zJwwiBPatCLlsTKs7YnbpMW6EWYUjcp")
get_file(node = "fuhzx", file = "wcvp_names.csv", path = "data", 
               remote_path = "data")


my_project <- osf_retrieve_node("https://osf.io/fuhzx/")
osf_retrieve_file("fuhzx")
osf_ls_files(my_project)
osf_retrieve_user("me")
