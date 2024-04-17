plants_full <- fread("data/wcvp_accepted_merged.txt")
growth <- read.csv("data/growths_rf.csv") %>% 
  dplyr::select(growth_form = apd_plant_growth_form, plant_name_id)


aa <- plants_full %>%
  filter(family == "Musaceae")

table(aa$growth_form)
table(aa$)


aa <- plants_full %>%
  filter(family == "Orchidaceae")
table(aa$growth_form)

aa <- plants_full %>%
  filter(family == "Poaceae")
table(aa$growth_form)

aa <- plants_full %>%
  filter(family == "Solanaceae")
table(aa$growth_form)

aa <- plants_full %>%
  filter(genus == "Solanum")
table(aa$growth_form)

aa <- plants_full %>%
  filter(genus == "Poa")

table(aa$growth_form)

table(aa$growth_form)
