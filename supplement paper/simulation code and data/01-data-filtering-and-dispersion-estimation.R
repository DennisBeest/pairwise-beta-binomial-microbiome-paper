# rm(list = ls())
library(tidyverse)
library(corncob)
source("custom-functions/functions-data-selection.R")

## Load Data -------------------------------------------------------------------
(infant_data <- readRDS(file = "data/infant_data.RDS"))
# Infant data have been arranged in a long-format data frame. This data-frame
# contains all information used throughout the simulations. 

## Visualize Sparsity Across Both studies --------------------------------------
# sparsity_summary_per_study <- 
#   infant_data %>% 
#   group_by(taxon, study_id) %>% 
#   summarize(zero_rate = sum(count == 0) / n()) %>%
#   arrange(., zero_rate)
# ggplot(sparsity_summary_per_study) +
#   geom_point(aes(y = as.factor(taxon), x = zero_rate)) +
#   facet_wrap(~study_id) + labs("Zero Rate per taxon across the studies") +
#   theme(axis.text.y = element_text(size = 5)) 

## Create Sparsity summary pooled over both studies ----------------------------
(sparsity_summary_pooled <- 
  infant_data %>% 
  group_by(taxon) %>% 
  summarize(zero_rate = sum(count == 0) / n()) %>%
  mutate(., zero_rate = round(zero_rate, 2)) %>%
  arrange(., zero_rate))

## Filter Taxa -----------------------------------------------------------------
## Any taxa at more than 80% zeros are removed from the dataset. 
selected_taxa <- filter(sparsity_summary_pooled, zero_rate <= 0.8) %>% 
  pull(taxon)
infant_data_filtered <- filter(infant_data, taxon %in% selected_taxa)
saveRDS(object = infant_data_filtered, file = "data/infant_data_filtered.RDS")

## Create Pairwise Taxon Tables ------------------------------------------------
## For Bifidobacterium
ratio_table_bifido <- 
  .generate_pairwise_table_reference_category(infant_data_filtered, 
                                              "Bifidobacterium")
saveRDS(object = ratio_table_bifido, file = "data/ratio_table_bifido.RDS")

## For Anaerostipes
ratio_table_anaero <- 
  .generate_pairwise_table_reference_category(infant_data_filtered, 
                                              "Anaerostipes")
saveRDS(object = ratio_table_anaero, file = "data/ratio_table_anaero.RDS")


## Dispersion Estimation - Bifidobacterium -------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(7347843)

# Turn study_id into a factor for model fitting
dataset <- 
  mutate(ratio_table_bifido, study = as_factor(study_id)) %>%
  select(-c(study_id))
n_ratios <- length(unique(dataset$ratio_id))
out_list <- vector(mode = "list", length = n_ratios)
for (i in 1:n_ratios) {
  # Select Data for a specific taxon pair
  tmp_data <- filter(dataset, ratio_id == i)
  
  # Save random seed state
  tmp_seed <- .Random.seed
  
  # Fit Model
  # More elaborate try-catch: only returns fitted parameters if model 
  # trustworthy in the sense of no separation detected. Moreover, a non-critical
  # version warning is suppressed. Initial values are provided to avoid
  # common inital value finding problems.
  fit  <- tryCatch({
    pkgcond::suppress_warnings({bbdml(formula = cbind(num_count, denom_count) ~ study,
                                      phi.formula = ~ 1, data = tmp_data,
                                      inits = matrix(data = c(0,0,0),
                                                     nrow = 1, ncol = 3),
                                      link = "logit", phi.link = "log")}, 
                               pattern = "(will be removed from 'brglm2' at version 0.8.)")},
    error = function(e){return(NA)})
  
  # Save Outputs if available
  if (anyNA(summary(fit)$coefficients)) {
    out_df <- tibble(mu_1 = NA, mu_2 = NA, mu_average = NA, 
                     dispersion = NA, 
                     precision = NA, alpha = NA, beta = NA, ratio_id = i)
    out_df$model_fit <- list(NA)
    out_df$used_data <- list(tmp_data)
    out_df$random_seed <- list(tmp_seed)
    out_list[[i]] <- out_df
  } else {
    # Add output vector to list
    tmp <- summary(fit)
    tmp$coefficients
    # Compute model parameters & dispersion
    mu_1 <- invlogit(tmp$coefficients[1,1])
    mu_2 <- invlogit(tmp$coefficients[1,1] + tmp$coefficients[2,1])
    
    mu_average <- (mu_1 + mu_2) / 2
    dispersion <- exp(tmp$coefficients[3,1])
    precision <- dispersion^-1
    alpha <- mu_average * precision
    beta <- (1-mu_average) * precision
    # numeric outcomes + seed + fitted model
    out_df <- tibble(mu_1, mu_2, mu_average, dispersion, 
                     precision, alpha, beta, ratio_id = i)
    out_df$model_fit <- list(fit)
    out_df$used_data <- list(tmp_data)
    out_df$random_seed <- list(tmp_seed)
    out_list[[i]] <- out_df
    # Remove temporary variables
    rm(list = c("out_df", "tmp", "mu_1", "mu_2", "mu_average",
                "dispersion", "precision", "alpha", "beta", "tmp_seed"))
  }
  rm(list = c("fit"))
}
fitted_parameters_bifido <- bind_rows(out_list)
saveRDS(object = fitted_parameters_bifido, file = "data/fitted-parameters-bifido.RDS")

# Compute and save dispersion quartiles
(quartiles_bifido <- tibble(dispersion_value = quantile(fitted_parameters_bifido$dispersion, probs = c(0.25, 0.5, 0.75)),
                           quartile = c("lower", "median", "upper")))
saveRDS(object = quartiles_bifido, file = "data/quartiles-bifido.RDS")

## Dispersion Estimation - Anaerostipes ----------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(7347843)

# Turn study_id into a factor for model fitting
dataset <- 
  mutate(ratio_table_anaero, study = as_factor(study_id)) %>%
  select(-c(study_id))
n_ratios <- length(unique(dataset$ratio_id))
out_list <- vector(mode = "list", length = n_ratios)
unique(dataset$study)
for (i in 1:n_ratios) {
  # Select Data for a specific taxon pair
  tmp_data <- filter(dataset, ratio_id == i)
  
  # Save random seed state
  tmp_seed <- .Random.seed
  
  # Fit Model
  # More elaborate try-catch: only returns fitted parameters if model 
  # trustworthy in the sense of no separation detected. Moreover, a non-critical
  # version warning is suppressed. Initial values are provided to avoid
  # common inital value finding problems.
  fit  <- tryCatch({
    pkgcond::suppress_warnings({bbdml(formula = cbind(num_count, denom_count) ~ study,
                                      phi.formula = ~ 1, data = tmp_data,
                                      inits = matrix(data = c(0,0,0),
                                                     nrow = 1, ncol = 3),
                                      link = "logit", phi.link = "log")}, 
                               pattern = "(will be removed from 'brglm2' at version 0.8.)")},
    error = function(e){return(NA)})
  
  # Save output if available
  if (anyNA(summary(fit)$coefficients)) {
    out_df <- tibble(mu_1 = NA, mu_2 = NA, mu_average = NA, 
                     dispersion = NA, 
                     precision = NA, alpha = NA, beta = NA, ratio_id = i)
    out_df$model_fit <- list(NA)
    out_df$used_data <- list(tmp_data)
    out_df$random_seed <- list(tmp_seed)
    out_list[[i]] <- out_df
  } else {
    # Add output vector to list
    tmp <- summary(fit)
    tmp$coefficients
    # Compute model parameters & dispersion
    mu_1 <- invlogit(tmp$coefficients[1,1])
    mu_2 <- invlogit(tmp$coefficients[1,1] + tmp$coefficients[2,1])
    
    mu_average <- (mu_1 + mu_2) / 2
    dispersion <- exp(tmp$coefficients[3,1])
    precision <- dispersion^-1
    alpha <- mu_average * precision
    beta <- (1-mu_average) * precision
    # numeric outcomes + seed + fitted model
    out_df <- tibble(mu_1, mu_2, mu_average, dispersion, 
                     precision, alpha, beta, ratio_id = i)
    out_df$model_fit <- list(fit)
    out_df$used_data <- list(tmp_data)
    out_df$random_seed <- list(tmp_seed)
    out_list[[i]] <- out_df
    # Remove temporary variables
    rm(list = c("out_df", "tmp", "mu_1", "mu_2", "mu_average",
                "dispersion", "precision", "alpha", "beta", "tmp_seed"))
  }
  rm(list = c("fit"))
}
fitted_parameters_anaero <- bind_rows(out_list)
saveRDS(object = fitted_parameters_anaero, file = "data/fitted-parameters-anaero.RDS")

# Compute and save dispersion quartiles
(quartiles_anaero <- tibble(dispersion_value = 
                             quantile(fitted_parameters_anaero$dispersion, 
                                      probs = c(0.25, 0.5, 0.75)),
                           quartile = c("lower", "median", "upper")))
saveRDS(object = quartiles_anaero, file = "data/quartiles-anaero.RDS")






