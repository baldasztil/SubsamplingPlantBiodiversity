library(DEoptim)
library(progress)
library(tidyverse)
library(data.table)
library(sf)
library(foreach)
library(doParallel)

# Working directory ------------------------------------------------------------
# Import data ------------------------------------------------------------------
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

richness_patterns <- fread("data/richness_patterns.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")


sample <- sample(plants_full$plant_name_id, 10000)

plants_sample <- plants_full %>%  
  filter(plant_name_id %in% sample)

dist_sample <- dist_native %>%  
  filter(plant_name_id %in% sample)

# creating objects for the function

plantlist_names <- plants_full %>% 
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)




richness_patterns_con <- plantlist_dist %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- plantlist_dist %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- plantlist_dist %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

richness_patterns1 <- richness_patterns_con %>% 
  rename( LEVEL_COD = LEVEL1_COD)
richness_patterns2 <- richness_patterns_reg %>% 
  rename( LEVEL_COD = LEVEL2_COD)
richness_patterns3 <- richness_patterns_bru %>% 
  rename( LEVEL_COD = LEVEL3_COD)

richness_patterns <- rbind(richness_patterns1,richness_patterns2, 
                           richness_patterns3)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

species_indices <-  seq(1, 1000, 10)
# Define fitness function (absolute difference from correlation target)
fitness <- function(species_indices) {
  species <- plantlist_names[species_indices, ]
  dist <- dist_native %>%
    filter(plant_name_id %in% species$plant_name_id)
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
  rich_rel <- rich_overall_bru %>%
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    replace(is.na(.), 0) %>%
    mutate(sp = length(species_indices))
  cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "spearman", exact = FALSE)[[4]]
  cor.pe <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "pearson", exact = FALSE)[[4]]
  abs(0.95 - max(cor.sp, cor.pe))  # Minimize the absolute difference from 0.95
}

# Define the optimization problem with progress bar
pb <- progress_bar$new(format = "  [:bar] :percent :eta", total = 100)

GA_result <- ga(type = "permutation", fitness = fitness, 
                lower = 1, upper = nrow(plantlist_names), maxiter = 100)

# Extract optimal subset of species

optimal_species_names <- plantlist_names$plant_name_id[GA_result@solution]

xx <- GA_result@solution 
optimal_species_indices <- GA_result@solution

# Get the corresponding species data
optimal_species <- plantlist_names[plantlist_names$plant_name_id %in% optimal_species_names, ]

summary(GA_result)

min(xx)

fitness[3]
fitness(3)

min(GA_result)

which(xx == min(xx), arr.ind = TRUE)

simulated_annealing <- function(initial_solution, max_iterations, temperature, cooling_rate) {
  current_solution <- initial_solution
  best_solution <- current_solution
  for (i in 1:max_iterations) {
    temperature <- temperature * cooling_rate
    new_solution <- perturb(current_solution)  # Function to generate a new solution
    if (accept(new_solution, current_solution, temperature)) {
      current_solution <- new_solution
      if (fitness(new_solution) < fitness(best_solution)) {
        best_solution <- new_solution
      }
    }
  }
  return(best_solution)
}
perturb <- function(current_solution) {
  # Generate a new solution based on your problem-specific perturbation logic
  # Example: Randomly select a subset of species and return the indices
}


accept <- function(new_solution, current_solution, temperature) {
  delta_fitness <- fitness(new_solution) - fitness(current_solution)
  if (delta_fitness < 0) {
    return(TRUE)
  } else {
    acceptance_prob <- exp(-delta_fitness / temperature)
    return(runif(1) < acceptance_prob)
  }
}

initial_solution <- # Define your initial solution (e.g., randomly selected species)
  max_iterations <- 10000  # Adjust as needed
temperature <- 1.0  # Initial temperature
cooling_rate <- 0.999  # Adjust as needed

best_solution <- simulated_annealing(initial_solution, max_iterations, temperature, cooling_rate)







# Load DEoptim package (if not already loaded)
library(DEoptim)


fitness <- function(species_indices) {
  species <- plantlist_names[species_indices, ]
  dist <- dist_native %>%
    filter(plant_name_id %in% species$plant_name_id)
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
  rich_rel <- rich_overall_bru %>%
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    replace(is.na(.), 0) %>%
    mutate(sp = length(species_indices))
  cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "spearman", exact = FALSE)[[4]]
  cor.pe <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "pearson", exact = FALSE)[[4]]
  abs(0.95 - max(cor.pe))  # Minimize the absolute difference from 0.95
}


# Define optimization problem
DE_result <- DEoptim(fitness, lower = rep(1, nrow(plantlist_names)), upper = rep(nrow(plantlist_names), nrow(plantlist_names)))

# Extract optimal solution
optimal_species_indices <- DE_result$optim$bestmem
optimal_species <- plantlist_names[optimal_species_indices, ]






subsample_greedy <- function(spec_n, target_correlation) {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
  # ... (rest of your code)
  
  # Define a function to calculate the correlation given a subset of species
  calculate_correlation <- function(subset_indices) {
    species <- plantlist_names[subset_indices, ]
    species <- plantlist_names[species_indices, ]
    dist <- dist_native %>%
      filter(plant_name_id %in% species$plant_name_id)
    sample_rich_bru <- dist %>%
      group_by(area_code_l3) %>%
      summarise(richness_sample = n_distinct(plant_name_id)) %>%
      rename(LEVEL_COD = area_code_l3)
    rich_rel <- rich_overall_bru %>%
      left_join(sample_rich_bru, by = "LEVEL_COD") %>%
      replace(is.na(.), 0) %>%
      mutate(sp = length(species_indices))
    cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                       method = "spearman", exact = FALSE)[[4]]
    cor.pe <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                       method = "pearson", exact = FALSE)[[4]]
    cor.sp
    
    # ... (Calculate correlation based on your criteria)
  }
  
  # Define a function for the greedy algorithm
  greedy_algorithm <- function(target_correlation) {
    selected_indices <- numeric(0)
    current_correlation <- 0
    
    while (current_correlation < target_correlation) {
      best_species <- NULL
      best_correlation <- -Inf
      
      for (i in setdiff(1:nrow(plantlist_names), selected_indices)) {
        temp_subset <- c(selected_indices, i)
        temp_correlation <- calculate_correlation(temp_subset)
        
        if (temp_correlation > best_correlation) {
          best_species <- i
          best_correlation <- temp_correlation
        }
      }
      
      if (!is.null(best_species)) {
        selected_indices <- c(selected_indices, best_species)
        current_correlation <- best_correlation
      } else {
        break  # No species improved correlation, exit loop
      }
    }
    
    return(list(selected_indices = selected_indices, correlation = current_correlation))
  }
  
  # Apply the greedy algorithm
  result <- greedy_algorithm(target_correlation)
  
  # Extract the selected species
  selected_species <- plantlist_names[result$selected_indices, ]
  
  # ... (rest of your code)
}




subsample_greedy <- function() {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
  # ... (rest of your code)
  
  # Define a function to calculate the correlation given a subset of species
  calculate_correlation <- function(subset_indices) {
    species <- plantlist_names[species_indices, ]
    dist <- dist_native %>%
      filter(plant_name_id %in% species$plant_name_id)
    sample_rich_bru <- dist %>%
      group_by(area_code_l3) %>%
      summarise(richness_sample = n_distinct(plant_name_id)) %>%
      rename(LEVEL_COD = area_code_l3)
    rich_rel <- rich_overall_bru %>%
      left_join(sample_rich_bru, by = "LEVEL_COD") %>%
      replace(is.na(.), 0) %>%
      mutate(sp = length(species_indices))
    cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                       method = "spearman", exact = FALSE)[[4]]
    cor.pe <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                       method = "pearson", exact = FALSE)[[4]]
    cor.sp
    # ... (Calculate correlation based on your criteria)
  }
  
  # Define a function for the greedy algorithm
  greedy_algorithm <- function(target_correlation) {
    selected_indices <- numeric(0)
    current_correlation <- 0
    
    while (current_correlation < target_correlation) {
      best_species <- NULL
      best_correlation <- -Inf
      
      for (i in setdiff(1:nrow(plantlist_names), selected_indices)) {
        temp_subset <- c(selected_indices, i)
        temp_correlation <- calculate_correlation(temp_subset)
        
        if (temp_correlation > best_correlation) {
          best_species <- i
          best_correlation <- temp_correlation
        }
      }
      
      if (!is.null(best_species)) {
        selected_indices <- c(selected_indices, best_species)
        current_correlation <- best_correlation
      } else {
        break  # No species improved correlation, exit loop
      }
    }
    
    return(list(selected_indices = selected_indices, correlation = current_correlation))
  }
  
  # Calculate correlation for increasing numbers of species
  correlations <- numeric(nrow(plantlist_names))
  for (i in 1:nrow(plantlist_names)) {
    result <- greedy_algorithm(0.95)  # Assuming target correlation is 0.95
    correlations[i] <- result$correlation
  }
  
  # ... (rest of your code)
}


install.packages("greedy")
library(greedy)

# Define a function to calculate the correlation given a subset of species

calculate_correlation <- function(subset_indices) {
  species <- plantlist_names[subset_indices, ]
  dist <- dist_native %>%
    filter(plant_name_id %in% species$plant_name_id)
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
  rich_rel <- rich_overall_bru %>%
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    replace(is.na(.), 0) %>%
    mutate(sp = length(subset_indices))
  cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "spearman", exact = FALSE)[[4]]
  cor.sp
  
  # ... (Calculate correlation based on your criteria)
}

# Define a function to perform greedy optimization
# Define a function to calculate the correlation given a subset of species
calculate_correlation <- function(subset_indices) {
  # (Existing code for calculating correlation)
}

# Define a function to perform greedy optimization
subsample_greedy <- function() {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
  # Initialize variables
  selected_indices <- numeric(0)
  current_correlation <- 0
  target_correlation <- 0.95  # Set your desired target correlation here
  
  # Main loop for greedy optimization
  while (length(selected_indices) < nrow(plantlist_names)) {
    best_species <- NULL
    best_correlation <- -Inf
    
    # Iterate through species not yet selected
    for (i in setdiff(1:nrow(plantlist_names), selected_indices)) {
      temp_subset <- c(selected_indices, i)
      temp_correlation <- calculate_correlation(temp_subset)
      
      # Check if this species improves correlation and meets target
      if (temp_correlation > best_correlation && temp_correlation >= target_correlation) {
        best_species <- i
        best_correlation <- temp_correlation
      }
    }
    
    # If a species is found that improves correlation and meets target, add it to selection
    if (!is.null(best_species)) {
      selected_indices <- c(selected_indices, best_species)
      current_correlation <- best_correlation
      
      # Print the number of species included in the sample
      cat("Added species:", i, "\n")
      cat("Number of species included:", length(selected_indices), "\n")
    } else {
      break  # No species improved correlation or reached the target, exit loop
    }
  }
  
  # The selected_indices now contain the optimal subset of species
  
  # ... (rest of your code)
}
subsample_greedy()



# Load necessary libraries
install.packages("dplyr")
install.packages("corrplot")
library(dplyr)
library(corrplot)

# Generate sample data (replace this with your actual dataframe)
set.seed(123)
species <- sample(c("Species1", "Species2", "Species3", "Species4"), 100, replace = TRUE)
countries <- sample(c("CountryA", "CountryB", "CountryC"), 100, replace = TRUE)
df <- data.frame(Species = species, Country = countries)

# Function to calculate Spearman correlation coefficient
calculate_spearman <- function(df) {
  return(cor.test(df$SpeciesRichness, df$SampledRichness, method = "spearman")$estimate)
}

# Initialize variables
original_species_richness <- df %>%
  group_by(Country) %>%
  summarise(SpeciesRichness = n_distinct(Species)) %>%
  ungroup()

sampled_species_richness <- data.frame(Country = character(0), SampledRichness = numeric(0))

# Initialize correlation
current_correlation <- 0

# Iterate until desired correlation is reached
while (current_correlation < 0.95) {
  remaining_species <- setdiff(df$Species, sampled_species_richness$Species)
  
  # Calculate correlation for each remaining species
  correlations <- sapply(remaining_species, function(species) {
    sampled_df <- bind_rows(sampled_species_richness, data.frame(Species = species))
    full_df <- inner_join(original_species_richness, sampled_df, by = "Country")
    new_correlation <- calculate_spearman(full_df)
    return(new_correlation)
  })
  
  # Find the species that increases the correlation the most
  best_species <- remaining_species[which.max(correlations)]
  
  # Add the best species to the sampled list
  sampled_species_richness <- bind_rows(sampled_species_richness, data.frame(Species = best_species))
  
  # Calculate new correlation
  full_df <- inner_join(original_species_richness, sampled_species_richness, by = "Country")
  current_correlation <- calculate_spearman(full_df)
}

# Print results
cat("Desired correlation (0.95) reached with", nrow(sampled_species_richness), "species.\n")
print(sampled_species_richness)



















# Load necessary libraries
install.packages("dplyr")
install.packages("corrplot")
library(dplyr)
library(corrplot)

# Generate sample data (replace this with your actual dataframe)
set.seed(123)
species <- sample(c("Species1", "Species2", "Species3", "Species4"), 100, replace = TRUE)
countries <- sample(c("CountryA", "CountryB", "CountryC"), 100, replace = TRUE)
df <- data.frame(Species = species, Country = countries)

# Function to calculate Spearman correlation coefficient
calculate_spearman <- function(df) {
  return(cor.test(df$SpeciesRichness, df$SampledRichness, method = "spearman")$estimate)
}

# Initialize variables
original_species_richness <- df %>%
  group_by(Country) %>%
  summarise(SpeciesRichness = n_distinct(Species)) %>%
  ungroup()

sampled_species_richness <- data.frame(Country = character(0), SampledRichness = numeric(0))

# Initialize correlation
current_correlation <- 0

# Iterate until desired correlation is reached
while (current_correlation < 0.95) {
  remaining_species <- setdiff(df$Species, sampled_species_richness$Species)
  
  # Calculate correlation for each remaining species
  correlations <- sapply(remaining_species, function(species) {
    sampled_df <- bind_rows(sampled_species_richness, data.frame(Species = species))
    full_df <- inner_join(original_species_richness, sampled_df, by = "Country")
    new_correlation <- calculate_spearman(full_df)
    return(new_correlation)
  })
  
  # Find the species that increases the correlation the most
  best_species <- remaining_species[which.max(correlations)]
  
  # Add the best species to the sampled list
  sampled_species_richness <- bind_rows(sampled_species_richness, data.frame(Species = best_species))
  
  # Calculate new correlation
  full_df <- inner_join(original_species_richness, sampled_species_richness, by = "Country")
  current_correlation <- calculate_spearman(full_df)
}

# Print results
cat("Desired correlation (0.95) reached with", nrow(sampled_species_richness), "species.\n")
print(sampled_species_richness)





# Define the objective function (correlation to be maximized)
objective <- function(species_indices) {
  # Your code to calculate correlation based on selected species
}

# Define lower and upper bounds for each species (dimension)
lower_bound <- rep(1, length(remaining_species))  # Start from the first species
upper_bound <- rep(length(remaining_species), length(remaining_species))  # Go up to the last species

# Initialize DE parameters
np <- 10  # Start with a population of 10 candidate solutions
iter <- 100  # Allow the algorithm to run for up to 100 iterations

# Run DE algorithm
result <- DEoptim(objective, lower = lower_bound, upper = upper_bound,
                  DEoptim.control(NP = np, itermax = iter))



objective <- function(subset_indices) {
  species <- plantlist_names[subset_indices, ]
  dist <- dist_native %>%
    filter(plant_name_id %in% species$plant_name_id)
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
  rich_rel <- rich_overall_bru %>%
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    replace(is.na(.), 0) %>%
    mutate(sp = length(subset_indices))
  cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "spearman", exact = FALSE)[[4]]
  cor.sp
  
  # ... (Calculate correlation based on your criteria)
}



library(DEoptim)

objective <- function(subset_indices) {
  species <- slice_sample(n = round(subset_indices),plantlist_names)
  dist <- plantlist_dist %>%
    filter(plant_name_id %in% species$plant_name_id)
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
  rich_rel <- rich_overall_bru %>%
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    replace(is.na(.), 0) %>%
    mutate(sp = length(species))
  cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "spearman", exact = FALSE)[[4]]
    return(-cor.sp)  # Return negative correlation (since DEoptim minimizes)
    #return(ifelse(correlation >= 0.95, correlation, 0))
  }






# Define lower and upper bounds for each species (dimension)
lower_bound <- 1  # Start from the first species
upper_bound <- nrow(plantlist_names) # Go up to the last species

# Initialize DE parameters
np <- 10  # Start with a population of 10 candidate solutions
iter <- 100  # Allow the algorithm to run for up to 100 iterations

# Run DE algorithm

frame <- data.frame(matrix(ncol = 3, nrow = 1391))
names(frame) <- c("Species", "bestcorr", "members")

for (i in 1:1390) {
  frame$Species <- i
  result <- DEoptim(objective, lower = i, upper = i,
                  DEoptim.control(NP = np, itermax = 10))
  frame$bestcorr <- result$optim$bestval
  frame$members <- result$optim$bestmemit
}



library(foreach)
library(doParallel)
# Define lower and upper bounds for each species (dimension)
lower_bound <- 1  # Start from the first species
upper_bound <- nrow(plantlist_names) # Go up to the last species

# Initialize DE parameters
np <- 10  # Start with a population of 10 candidate solutions
iter <- 100  # Allow the algorithm to run for up to 100 iterations

# Run DE algorithm

results_list <- vector("list", length = 1390)

# Define the number of cores you want to use
num_cores <- 3  # Adjust this based on your system

# Create a cluster
cl <- makeCluster(num_cores)

# Register the cluster with foreach
registerDoParallel(cl)

# Rest of your code using foreach...


# Create a foreach loop to parallelize the loop
foreach(i = 1:1390) %do% {
  frame <- data.frame(matrix(ncol = 3, nrow = 1))
  names(frame) <- c("Species", "bestcorr", "members")
  de_result <- DEoptim(objective, lower = i, upper = i,
                    DEoptim.control(NP = np, itermax = 5))
  
  frame$Species <- i
  frame$bestcorr <- de_result$optim$bestval * -1
  frame$members <- de_result$optim$bestmemit
  
  # Store the result in the list
  results_list[[i]] <- frame
}


# Combine the results into a single data frame
final_frame <- do.call(rbind, results_list)

# Clean up the cluster after finishing the loop
stopCluster(cl)


library(GA)
list()
for (i in 1:1390) {
  ga_result <- ga(type = "real-valued", 
                fitness = objective, 
                lower = i,  # Minimum sample size
                upper = i,  # Maximum sample size
                maxiter = 10)  # Number of iteration
  list[[i]] <- ga_result
}

# Initialize variables to keep track of best correlation and corresponding subset
best_correlation <- 0
best_subset <- NULL
min_species <- Inf  # Initialize with a high value

# Run DE algorithm
for (i in 1:1) {
  result <- DEoptim(objective, lower = lower_bound, upper = upper_bound,
                    DEoptim.control(NP = np, itermax = 1))  # Run DE for one iteration
  
  # Extract the best solution and its correlation
  best_species_indices <- result$optim$bestmem
  selected_species <- plantlist_names[best_species_indices, ]
  best_correlation <- -result$optim$value  # Convert to positive since DEoptim minimizes
  
  if (best_correlation >= 0.95 && length(selected_species) < min_species) {
    best_subset <- selected_species
    min_species <- length(selected_species)
  }
  
  # Print the best correlation at this iteration
  cat("Iteration", i, "Best Correlation:", best_correlation, "\n")
}

# Print the minimum number of species needed to achieve desired correlation
cat("Minimum Number of Species Needed:", length(best_subset), "\n")
cat("Corresponding Subset of Species:\n")
print(best_subset)

objective(1000)



library(GA)




objective <- function(subset_indices) {
  #print(paste0(subset_indices))
  subset_indices <- round(subset_indices)
  species <- slice_sample(n = subset_indices,plantlist_names)
  dist <- plantlist_dist %>%
    filter(plant_name_id %in% species$plant_name_id)
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
  rich_rel <- rich_overall_bru %>%
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    replace(is.na(.), 0) %>%
    mutate(sp = length(species))
  cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "spearman", exact = FALSE)[[4]]
  if (cor.sp >= 0.95) {
    return(-cor.sp)  # Return negative correlation (since DEoptim minimizes)
  } else {
    return(1000)  # Return a high cost for correlations below 0.95
  }
}





fitness <- function(subset_indices) {
  cor_values <- apply(as.matrix(subset_indices), 1, function(indices) {
    correlation <- objective(indices)
    return(ifelse(correlation >= 0.95, correlation, 0))
  })
  return(cor_values)
}

# Set up GA parameters
num_species <- 10000  # Assuming you have 100 species, adjust this value accordingly
GA_params <- list(
  type = "permutation",  # Discrete optimization
  fitness = fitness,
  nBits = num_species,  # Number of species
  popSize = 50,  # Population size
  maxiter = 10,  # Maximum number of iterations
  run = 1,  # Number of runs
  seed = 123  # Seed for reproducibility
)


# Run the GA
ga_result <- ga(type = GA_params$type, fitness = GA_params$fitness, 
                nBits = GA_params$nBits, popSize = GA_params$popSize, 
                maxiter = GA_params$maxiter, run = GA_params$run, 
                seed = GA_params$seed, lower = 1, upper = 10000, keepBest = T)


# Extract the best solution and its correlation
best_species_indices <- ga_result@solution
best_species <- plantlist_names[best_species_indices, ]
best_correlation <- objective(best_species_indices)

# Print the results
cat("Best Correlation:", best_correlation, "\n")
cat("Corresponding Subset of Species:\n")
print(best_species)

