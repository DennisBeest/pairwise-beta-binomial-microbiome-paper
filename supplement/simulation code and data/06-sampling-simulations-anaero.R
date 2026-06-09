# Load packages and custom functions -------------------------------------------
library(tidyverse)
library(corncob)
library(doParallel)
source("custom-functions/functions-set-up-sim.R")
source("custom-functions/functions-run-sim-sampling.R")

# Load reference data ----------------------------------------------------------
ratio_table <- readRDS("data/ratio_table_anaero.RDS")
totals_list <- .generate_pairwise_totals_list(ratio_table)
raw_data_list <- .generate_pairwise_raw_data_list(ratio_table)

# Generate simulation input objects --------------------------------------------
scenario_table <- 
  .make_scenario_table_resampling(ratio_ids = seq(1:50), 
                                  n_samples = c(20, 40, 60, 100, 200), 
                                  methods = c("bbglm", "lrlm", 
                                              "wilcox", "perm"),
                                  n_sim = 10000)

filename_1 <- "output/sampling-scenario-table-anaero.RDS"
filename_2_prefix <- "output/sampling-anaero-scenario-"

saveRDS(file = filename_1, object = scenario_table)

# Run simulations --------------------------------------------------------------
# Set-up random streams over 6 cores
RNGkind("L'Ecuyer-CMRG")
registerDoParallel(cores = 6)
set.seed(551514845)

scenarios <- unique(scenario_table$scenario_id)
for (scenario_index in 1:length(scenarios)){
  selected_scenario <- scenarios[scenario_index]
  
  # For memory optimization and speed, each scenario is run separately and saved
  tmp <-
    filter(scenario_table, scenario_id == selected_scenario) %>%
    group_split(., sim_id)
  
  scenario_output <- foreach(j=1:length(tmp)) %dopar% {
    out <- .run_sim_resample(n_samples = tmp[[j]]$n_samples,
                             raw_ratio_data = raw_data_list[[tmp[[j]]$ratio_id]],
                             method = tmp[[j]]$method)
    out$scenario_id <- tmp[[j]]$scenario_id
    out$sim_id <- tmp[[j]]$sim_id
    out$n_sim <- tmp[[j]]$n_sim
    return(out)
  }
  saveRDS(bind_rows(scenario_output), 
          file = paste(filename_2_prefix,
                       scenario_index, ".RDS", sep = ""))
}
rm(list = ls())

