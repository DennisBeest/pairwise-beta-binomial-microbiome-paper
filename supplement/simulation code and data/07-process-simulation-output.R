rm(list = ls())
library("tidyverse")
source("functions/functions-simulation-output-processing.R")

# Compute Summary table for Anaerostipes reference run -------------------------
scenario_table <- 
  readRDS("output/parametric-scenario-table-anaero.RDS") %>%
  select(., -c("sim_id", "n_sim", "ratio_id")) %>% distinct(.)

prefix <- "output/parametric-anaero-scenario-"
scenarios <- seq(1,length(unique(scenario_table$scenario_id)))
filenames <- paste(prefix, scenarios, ".RDS", sep = "")
output_summary_list <- list()

for (scenario in 1:length(filenames)){
  # Load scenario specific simulation output
  file <- filenames[scenario]
  tmp <- readRDS(file) %>% left_join(., scenario_table, by = "scenario_id")
  # Make coverage summaries only for methods with confidence intervals
  method <- unique(tmp$method.x)
  # Compute summaries (method specific since coverage only computable for parametric methods)
  if (method %in% c("bbglm", "lrlm", "lrlm_2")){
    output_summary_list[[scenario]] <- .compute_summaries_parametric(tmp)
  } else {
    output_summary_list[[scenario]] <- .compute_summaries_nonparametric(tmp)
  }
}
output_summary_table <- bind_rows(output_summary_list)
output_summary_table <- left_join(output_summary_table, scenario_table, 
                                  by = "scenario_id")
saveRDS(output_summary_table, "output/anaero-parametric-summary.RDS")


# Compute Summary table for Bifidobacterium reference run ----------------------
scenario_table <- 
  readRDS("output/parametric-scenario-table-bifido.RDS") %>%
  select(., -c("sim_id", "n_sim", "ratio_id")) %>% distinct(.)

prefix <- "output/parametric-bifido-scenario-"
scenarios <- seq(1,length(unique(scenario_table$scenario_id)))
filenames <- paste(prefix, scenarios, ".RDS", sep = "")
output_summary_list <- list()

for (scenario in 1:length(filenames)){
  # Load scenario specific simulation output
  file <- filenames[scenario]
  tmp <-  readRDS(file) %>% left_join(., scenario_table, by = "scenario_id")
  # Make coverage summaries only for methods with confidence intervals
  method <- unique(tmp$method.x)
  # Compute summaries (method specific since coverage only computable for parametric methods)
  if (method %in% c("bbglm", "lrlm", "lrlm_2")){
    output_summary_list[[scenario]] <- .compute_summaries_parametric(tmp)
  } else {
    output_summary_list[[scenario]] <- .compute_summaries_nonparametric(tmp)
  }
}
output_summary_table <- bind_rows(output_summary_list)
output_summary_table <- left_join(output_summary_table, scenario_table, 
                                  by = "scenario_id")
saveRDS(output_summary_table, "output/bifido-parametric-summary.RDS")

# Compute Summary table for Bifidobacterium reference extended run -------------
scenario_table <- 
  readRDS("output/parametric-scenario-table-bifido-extended.RDS") %>%
  select(., -c("sim_id", "n_sim", "ratio_id")) %>% distinct(.)

prefix <- "output/parametric-bifido-extended-scenario-"
scenarios <- seq(10001,13600) # hardcoded based on sim-id modifications in script 03-par...
filenames <- paste(prefix, scenarios, ".RDS", sep = "")
output_summary_list <- list()

for (scenario in 1:length(filenames)){
  # Load scenario specific simulation output
  file <- filenames[scenario]
  tmp <- readRDS(file) %>% left_join(., scenario_table, by = "scenario_id")
  # Make coverage summaries only for methods with confidence intervals
  method <- unique(tmp$method.x)
  # Compute summaries (method specific since coverage only computable for parametric methods)
  if (method %in% c("bbglm", "lrlm", "lrlm_2")){
    output_summary_list[[scenario]] <- .compute_summaries_parametric(tmp)
  } else {
    output_summary_list[[scenario]] <- .compute_summaries_nonparametric(tmp)
  }
}
output_summary_table <- bind_rows(output_summary_list)
output_summary_table <- left_join(output_summary_table, scenario_table, 
                                  by = "scenario_id")
saveRDS(output_summary_table, "output/bifido-parametric-extended-summary.RDS")

# Compute Summary for Bifido-resampling run ------------------------------------
scenario_table <- 
  readRDS("output/sampling-scenario-table-bifido.RDS") %>%
  select(., -c("sim_id", "n_sim")) %>% distinct(.)

prefix <- "output/sampling-bifido-scenario-"
scenarios <- seq(1,length(unique(scenario_table$scenario_id)))
filenames <- paste(prefix, scenarios, ".RDS", sep = "")
output_summary_list <- list()

for (scenario in 1:length(filenames)){
  # Load scenario specific simulation output
  file <- filenames[scenario]
  tmp <- readRDS(file)
  output_summary_list[[scenario]] <- .compute_summaries_sampling(tmp)
}
output_summary_table <- bind_rows(output_summary_list)
output_summary_table <- left_join(output_summary_table, scenario_table, 
                                  by = "scenario_id")
saveRDS(output_summary_table, "output/bifido-sampling-summary.RDS")

# Compute Summary for Anaero-sampling run ---------------------------------------
scenario_table <- 
  readRDS("output/sampling-scenario-table-anaero.RDS") %>%
  select(., -c("sim_id", "n_sim")) %>% distinct(.)

prefix <- "output/sampling-anaero-scenario-"
scenarios <- seq(1,1000)
filenames <- paste(prefix, scenarios, ".RDS", sep = "")
output_summary_list <- list()

for (scenario in 1:length(filenames)){
  # Load scenario specific simulation output
  file <- filenames[scenario]
  tmp <- readRDS(file)
  output_summary_list[[scenario]] <- .compute_summaries_sampling(tmp)
}
output_summary_table <- bind_rows(output_summary_list)
output_summary_table <- left_join(output_summary_table, scenario_table, 
                                  by = "scenario_id")
saveRDS(output_summary_table, "output/anaero-sampling-summary.RDS")

