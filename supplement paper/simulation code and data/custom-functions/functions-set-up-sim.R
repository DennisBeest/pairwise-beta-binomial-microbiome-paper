.assert <- function (expr, error){if (! expr) stop(error, call. = FALSE)}
.odds_to_prob <- function(odds){return(prob <- odds / (1 + odds))}
.prob_to_odds <- function(prob){return(odds <- prob / (1 - prob))}
.make_scenario_table_parametric <- function(dispersion, odds, n, n_samples, 
                                            methods, ratio_ids,
                                            stepsize, max_steps){
  # Sequence of commands that expand a given baseline odds and dispersion 
  # vector into a scenario table.
  # Note:
  # ratio_ids are added for total sampling. if n_sim not a total of ratio_id,
  # the ratio ids will not be balanced, i.e. some ratio_ids will be used less
  # often than others. 
  if ((n %% length((ratio_ids)) != 0)){
    message("n_sim not a multiple of ratio_id. ratio_id usage will not be balanced.")}

  
  # Make Table of dispersions and odds combinations
  tmp_control_scenarios <- expand_grid(disp = dispersion, 
                                       odds = odds) %>%
    rowid_to_column(., "baseline_id")
  # Make power scenario tables via stepwise probability increasing
  tmp_scenario_list <- mapply(.make_alpha_beta_range, tmp_control_scenarios$odds, 
                              tmp_control_scenarios$disp, stepsize, max_steps, 
                              SIMPLIFY = FALSE)
  # Expand power scenarios tables to all methods and sample sizes and assign
  # unqiue scenario (unique combination of simulation setting and method) and 
  # sim identifiers (each scenario is repeated n_samples times)
  scenario_table <- bind_rows(tmp_scenario_list, .id = "baseline_id") %>%
    mutate(., baseline_id = as.integer(baseline_id)) %>%
    expand_grid(., method = methods, n_samples = n_samples) %>%
    rowid_to_column(., var = "scenario_id") %>%
    expand_grid(., n_sim = c(1:n)) %>%
    rowid_to_column(., var = "sim_id") %>%
    # add ratio_ids evenly to each scenario for total sampling
    group_by(., scenario_id) %>%
    mutate(., ratio_id = rep(unique(ratio_ids), length.out = n)) %>%
    ungroup()
  return(scenario_table)
}
.make_alpha_beta_range <- function(odds, disp, step = 0.05, max_steps = 10){
  # Function Outline
  # 1. Convert odds to probabilities
  # 2. Grow probabilities vector from value to boundary by stepsize step
  # 3. Cut probabilities vector to max_steps, default 10 <=> 0 to 0.45
  # 4. Place new probabilities into tibble as alpha and beta parameters using
  #    the provided dispersion. Add implied odds.
  #
  # Function checks validity of mu, and direction in which to grow such that
  # a maximal delta-probabilites can be achieved.
  mu <- .odds_to_prob(odds)
  disp <- disp
  
  # the conditional ensures that a range can actually be established
  .assert( mu >= 0 & mu <= 1, "mu must be a valid proportion between 0 and 1.")
  
  if (mu == 0.5) {
    # If mu == 0.5, make step range increasing until 1-step, starting at mu
    mu_list <- seq(mu, 1-step, by = step)
  } else if (mu <= 0.5) {
    # .assert( mu < (1*step) , "mu too close to boundary for given step size.")
    mu_list <- seq(mu, 1-step, by = step)
  } else if (mu > 0.5) {
    # If my greater than 0.5, 
    .assert( mu < 1- (1*step) , "mu too close to boundary for given step size.")
    mu_list <- seq(mu, 0 + step, by = -step)
  } 
  .assert(length(mu_list) != 1, 
          "Failed to grow mu_list. Check input configuration. Step may be too large.")
  # Cut mu_list down to max_steps
  if (length(mu_list) > max_steps) {
    mu_list <- mu_list[1:max_steps]
  }
  outcome_table <- tibble(mu_new = mu_list) %>%
    mutate(., 
           mu_control = mu,
           odds_control = odds, 
           disp_control = disp,
           alpha_control = mu * disp^-1,
           beta_control = (1-mu)*disp^-1,
           
           disp_new = disp, 
           alpha_new = mu_new * disp^-1, 
           beta_new = (1-mu_new) * disp^-1,
           mu_new = mu_new, 
           odds_new = .prob_to_odds(mu_new),
           odds_ratio = odds_new/odds_control)
  return(outcome_table)
}
.make_scenario_table_resampling <- function(ratio_ids, n_samples, 
                                            methods, n_sim){
  tmp_scenario_table <- expand_grid(ratio_id = ratio_ids,
                                    n_samples = as.integer(n_samples), 
                                    method = methods)
  output_table <- 
    rowid_to_column(tmp_scenario_table, var = "scenario_id") %>%
    expand_grid(., n_sim = c(1:n_sim)) %>%
    rowid_to_column(., var = "sim_id")
  return(output_table)
}
.generate_pairwise_totals_list <- function(pairwise_table){
  pairwise_table <- mutate(pairwise_table, ratio_id = as.integer(ratio_id),
                           totals = num_count + denom_count)
  ratio_ids <- unique(pairwise_table$ratio_id)
  n_ratios <- length(ratio_ids)
  out <- vector(mode = "list", length = n_ratios)
  for (i in 1:n_ratios){
    out[[i]] <- filter(pairwise_table, ratio_id == i) %>% pull(totals)
  }
  return(out)
}
.generate_pairwise_raw_data_list <- function(pairwise_table){
  pairwise_table <- mutate(pairwise_table, ratio_id = as.integer(ratio_id)) %>%
    mutate(., study_id = as.integer(as.factor(study_id)))
  ratio_ids <- unique(pairwise_table$ratio_id)
  n_ratios <- length(ratio_ids)
  out <- vector(mode = "list", length = n_ratios)
  for (i in 1:n_ratios){
    out[[i]] <- filter(pairwise_table, ratio_id == i) %>% 
      select(., c(sample_id, ratio_id, num_count, denom_count, study_id))
  }
  return(out)
}


