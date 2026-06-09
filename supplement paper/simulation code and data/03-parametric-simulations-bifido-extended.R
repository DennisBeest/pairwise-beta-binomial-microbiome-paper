## Load packages and custom functions -------------------------------------------
# rm(list = ls())
library(tidyverse)
library(corncob)
library(doParallel)
source("custom-functions/functions-set-up-sim.R")
source("custom-functions/functions-run-sim-parametric.R")

## Initialize Simulaton study --------------------------------------------------
## Load needed data & specify simulation parameters
ratio_table <- readRDS("data/ratio_table_bifido.RDS")
totals_list <- .generate_pairwise_totals_list(ratio_table)
dispersion_quartile_table <- readRDS("data/quartiles-bifido.RDS")
disp <- dispersion_quartile_table$dispersion_value
ratio_ids <- unique(ratio_table$ratio_id)
n_samples = c(20,40,60,80,100,200)
methods = c("bbglm", "perm", "wilcox", "lrlm", "lrlm_2")
n_sim = 1000
odds <- c(1/1, 1/9, 1/99, 1/999)
ratio_ids <- unique(ratio_table$ratio_id)
## Create scenario table with details on each simulation
scenario_table <- .make_scenario_table_parametric(disp, odds, n_sim, n_samples, 
                                                  methods, ratio_ids, 
                                                  stepsize = 0.005, 
                                                  max_steps = 10)
# Modify scenario_ids and sim_ids to avoid any overlap with primary simulations
scenario_table <- mutate(scenario_table, 
                         scenario_id = scenario_id + as.integer(10000))
scenario_table <- mutate(scenario_table, 
                         sim_id = sim_id + as.integer(1000000000))

# Initialize output filenames (output saved by scenario)
filename_1 <- "output/parametric-scenario-table-bifido-extended.RDS"
filename_2_prefix <- "output/parametric-bifido-extended-scenario-"
saveRDS(object = scenario_table, file = filename_1) 

# Run Simulations --------------------------------------------------------------
# Set-up random streams over 6 cores
RNGkind("L'Ecuyer-CMRG")
registerDoParallel(cores = 6)
set.seed(15548)

# For memory optimization and speed, each scenario is run separately and saved
scenarios <- unique(scenario_table$scenario_id)
for (scenario_index in 1:length(scenarios)){
  # Access scenario_index
  selected_scenario <- scenarios[scenario_index]
  
  # Make settings list for current scenario
  tmp <- filter(scenario_table, scenario_id == selected_scenario) %>%
    group_split(., sim_id)
  
  scenario_output <- foreach(j=1:length(tmp)) %dopar% {
    out <- .run_sim_para(alpha_control = tmp[[j]]$alpha_control, 
                         beta_control  = tmp[[j]]$beta_control, 
                         alpha_new = tmp[[j]]$alpha_new, 
                         beta_new = tmp[[j]]$beta_new,
                         n_samples  = tmp[[j]]$n_samples,
                         totals = totals_list[[tmp[[j]]$ratio_id]], # --> 
                         method = tmp[[j]]$method)
    out$scenario_id <- tmp[[j]]$scenario_id
    out$sim_id <- tmp[[j]]$sim_id
    out$n_sim <- tmp[[j]]$n_sim
    return(out)
  }
  saveRDS(bind_rows(scenario_output), 
          file = paste(filename_2_prefix, scenario_index, ".RDS", sep = ""))
}
rm(list = ls())
