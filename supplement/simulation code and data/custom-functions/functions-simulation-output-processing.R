.compute_summaries_parametric <- function(output_df){
  # Number of simulations
  n_sim <- dim(output_df)[1]
  
  # Number/rate of failed tests
  test_na_number <- sum(if_else(is.na(output_df$p_value) | 
                                  is.infinite(output_df$p_value), TRUE, FALSE))
  test_na_rate <- test_na_number / n_sim
  
  # Power ; conditional on model success
  power_conditional <- 
    sum(output_df$p_value <= 0.05, na.rm = TRUE) / (n_sim - test_na_rate)
  # Power ; based on total number of runs, model failures are failure to reject
  power_total <- sum(output_df$p_value <= 0.05, na.rm = TRUE) / (n_sim)
  
  # Number/rate of failed confidence interval estimations
  ci_na_number <- sum(if_else(is.na(output_df$effect_se) | 
                                is.infinite(output_df$effect_se), TRUE, FALSE))
  ci_na_rate <- ci_na_number / n_sim
  
  # Coverage
  mu_new <- unique(output_df$mu_new)
  mu_control <- unique(output_df$mu_control)
  true_effect <- log((mu_new / mu_control) / ((1-mu_new) / (1-mu_control)))
  
  tmp2 <- na.omit(select(output_df, effect_ci_lower, effect_ci_upper)) %>%
    mutate(., effect_covered = 
             if_else(true_effect >= effect_ci_lower & true_effect <= 
                       effect_ci_upper, TRUE, FALSE))
  
  # Coverage; conditional on model success
  coverage_conditional <- sum(tmp2$effect_covered) / (n_sim - ci_na_rate)
  
  # Coverage; based on total number of runs
  coverage_total <- sum(tmp2$effect_covered) / n_sim
  
  # Median ci_width
  median_ci_width <- tmp2 %>% 
    summarize(., median_ci_width = 
                median(effect_ci_upper - effect_ci_lower)) %>% 
    pull(median_ci_width)
  return(tibble(scenario_id = unique(output_df$scenario_id),
                test_na_rate, power_conditional, 
                power_total, ci_na_rate, 
                coverage_conditional, 
                coverage_total, median_ci_width))
}

.compute_summaries_nonparametric <- function(output_df){
  # Number of simulations
  n_sim <- dim(output_df)[1]
  
  # Number/rate of failed tests
  test_na_number <- sum(if_else(is.na(output_df$p_value) | 
                                  is.infinite(output_df$p_value), TRUE, FALSE))
  test_na_rate <- test_na_number / n_sim
  
  # Power ; conditional on model success
  power_conditional <- sum(output_df$p_value <= 0.05, na.rm = TRUE) / 
    (n_sim - test_na_rate)
  # Power ; based on total number of runs, model failures are failure to reject
  power_total <- sum(output_df$p_value <= 0.05, na.rm = TRUE) / (n_sim) 
  return(tibble(scenario_id = unique(output_df$scenario_id),
                test_na_rate, power_conditional, 
                power_total, ci_na_rate = NA, 
                coverage_conditional = NA, 
                coverage_total = NA,
                median_ci_width = NA))
}
.compute_summaries_sampling <- function(output_df){
  # Number of simulations
  n_sim <- dim(output_df)[1]
  
  # Number/rate of failed tests
  test_na_number <- sum(if_else(is.na(output_df$p_value) | 
                                  is.infinite(output_df$p_value), TRUE, FALSE))
  test_na_rate <- test_na_number / n_sim
  
  # Power ; conditional on model success
  type_I_error_rate_conditional <- 
    sum(output_df$p_value <= 0.05, na.rm = TRUE) / (n_sim - test_na_rate)
  # Power ; based on total number of runs, model failures are failure to reject
  type_I_error_rate_total <- 
    sum(output_df$p_value <= 0.05, na.rm = TRUE) / (n_sim)
  return(tibble(scenario_id = scenario,
                test_na_rate, type_I_error_rate_conditional, 
                type_I_error_rate_total))
}