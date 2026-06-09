#' Notes
#' 
#' --> .make_beta_bin_data assumes that total vector can be divided into two 
#' equal length parts; e.g. n_samples = 11 won't work.
#' 
#' In order for the list to table conversion to work, all simulation methods 
#' must produce the same output formats; if a lm produces an intercept 
#' estimate, then the permutation test must have a NaN or NA value for this
#' also.

.make_beta_bin_data <- function(alpha_control, beta_control, alpha_treat, beta_treat, n_samples, totals) {
  # Utilization:
  # Functions creates beta-binomial dataset using control and treatment specific
  # alpha and beta parameters.
  # Assumptions:
  # Assumes even number of samples in each group
  # Assumes supplied permuation of real totals or fixed total scalar
  tmp_control <- 
    .draw_beta_bin(alpha_control, beta_control, n_samples/2, 
                   totals[1:(n_samples/2)]) %>% 
    mutate(., treat = 0)
  tmp_treat <- 
    .draw_beta_bin(alpha_treat, beta_treat, n_samples/2, 
                   totals[((n_samples/2)+1):length(totals)]) %>% 
    mutate(., treat = 1)
  data <- dplyr::bind_rows(tmp_treat, tmp_control) %>% 
    mutate(., treat = as.factor(treat))
  return(data)
}
.draw_beta_bin <- function(alpha, beta, n_samples, totals){
  # Utilization:
  # Function draws a set of beta-binomial draws using supplied parameters and
  # totals. Output organized as a tibble with columns num_count & denom_count
  # Notes:
  # can deal with totals = 0: rbinom return zeros for totals = 0.
  # note that, if the n_samples < length(totals), rbinom reuses the specified
  # totals over; runs through the totals vector over and over.
  # consider adding error if n_samples != length(totals)
  probabilities <-  rbeta(n = n_samples, shape1 = alpha, shape2 = beta)
  successes   <- rbinom(n = n_samples, size = totals, prob = probabilities)
  bb_draws <- tibble(num_count = successes, denom_count = totals - successes)
  return(bb_draws)
}
.run_sim_para <- function(alpha_control, beta_control, alpha_new, beta_new, 
                          n_samples, totals, method){
  tmp_seed <- .Random.seed
  # Resample Totals
  resampled_totals <- sample(totals, size = n_samples, replace = TRUE)
  # Generate Dataset
  data <- .make_beta_bin_data(alpha_control, beta_control, alpha_new, beta_new,
                              n_samples, totals = resampled_totals)
  if (method == "lrlm") {
    out <- .fit_two_group_log_ratio_lm(data)
  } else if (method == "lrlm_2") {
    out <- .fit_two_group_log_ratio_lm_no_imputation(data)
  } else if (method == "bbglm") {
    out <- .fit_two_group_bb_lr_corncob(data)
  } else if (method == "perm") {
    out <- .perform_two_group_perm(data)
  }  else if (method == "wilcox") {
    out <- .perform_two_group_wilcox(data)
  } else {stop("No valid method specified.")}
  out$seed <- list(tmp_seed)
  return(out)
}
.perform_log_ratio_transformation <- function(input_data){
  # Imputation Performed
  out_data <- mutate(input_data, 
                     num_count = if_else(num_count == 0, 
                                         num_count + as.integer(1), 
                                         num_count),
                     denom_count = if_else(denom_count == 0, 
                                           denom_count + as.integer(1), 
                                           denom_count),
                     log_ratio = log(num_count / denom_count))
  return(out_data)
}
.perform_ratio_transformation <- function(input_data){
  # Imputation Performed
  out_data <- mutate(input_data, 
                     num_count = if_else(num_count == 0, 
                                         num_count + as.integer(1), 
                                         num_count),
                     denom_count = if_else(denom_count == 0, 
                                           denom_count + as.integer(1), 
                                           denom_count),
                     ratio = (num_count / denom_count))
  return(out_data)
}
.perform_log_ratio_transformation_no_imputation <- function(input_data){
  out_data <- mutate(input_data, log_ratio = log(num_count / denom_count))
  out_data <- filter(out_data, !is.na(log_ratio) & !is.infinite(log_ratio)) 
  return(out_data)
}

# Model Fitting ----------------------------------------------------------------
.fit_two_group_log_ratio_lm_no_imputation <-  function(data, nominal_alpha = 0.05){
  # Pseudo-count impute data
  data <- .perform_log_ratio_transformation_no_imputation(data)
  fit <- tryCatch({fit <- lm(formula = log_ratio ~ treat, data = data)},
                  error = function(e){return(NA)})
  # if model successful, make outcome table; otherwise make NA table
  if (!identical(fit, NA)) {
    # Create summary objects for accessing fitted model components
    tmp_summary <- summary(fit)$coefficients
    tmp_conf <- stats::confint.lm(fit, level = 1 - nominal_alpha)
    p_value <- tmp_summary[2 , 4]
    
    intercept_est <- tmp_summary[1, 1]
    intercept_se  <- tmp_summary[1, 2]
    intercept_ci_lower <- tmp_conf[1,1]
    intercept_ci_upper <- tmp_conf[1,2]
    
    effect_est <- tmp_summary[2, 1]
    effect_se  <- tmp_summary[2, 2]
    effect_ci_lower <- tmp_conf[2,1]
    effect_ci_upper <- tmp_conf[2,2]
    
    residual_se <- summary(fit)$sigma
    output <- tibble("p_value" = p_value,
                     "intercept_est" = intercept_est,
                     "intercept_se" = intercept_se, 
                     "intercept_ci_lower" = intercept_ci_lower,
                     "intercept_ci_upper" = intercept_ci_upper,
                     "effect_est" = effect_est, 
                     "effect_se" = effect_se,
                     "effect_ci_lower" = effect_ci_lower, 
                     "effect_ci_upper" = effect_ci_upper,
                     "residual_se" = residual_se,
                     "phi_est" = NA,
                     "phi_se" = NA,
                     "phi_ci_lower" = NA,
                     "phi_ci_upper" = NA,
                     "method" = "lrlm_2")
  } else {
    output <- tibble("p_value" = NA,
                     "intercept_est" = NA,
                     "intercept_se" = NA, 
                     "intercept_ci_lower" = NA,
                     "intercept_ci_upper" = NA,
                     "effect_est" = NA, 
                     "effect_se" = NA,
                     "effect_ci_lower" = NA, 
                     "effect_ci_upper" = NA,
                     "residual_se" = NA, 
                     "phi_est" = NA,
                     "phi_se" = NA,
                     "phi_ci_lower" = NA,
                     "phi_ci_upper" = NA,
                     "method" = "lrlm_2")
  }
  return(output)
}  
.fit_two_group_log_ratio_lm <- function(data, nominal_alpha = 0.05){
  # Pseudo-count impute data
  data <- .perform_log_ratio_transformation(data)
  fit <- tryCatch({fit <- lm(formula = log_ratio ~ treat, data = data)},
                  error = function(e){return(NA)})
  # if model successful, make outcome table; otherwise make NA table
  if (!identical(fit, NA)) {
    # Create summary objects for accessing fitted model components
    tmp_summary <- summary(fit)$coefficients
    tmp_conf <- stats::confint.lm(fit, level = 1 - nominal_alpha)
    p_value <- tmp_summary[2 , 4]
    
    intercept_est <- tmp_summary[1, 1]
    intercept_se  <- tmp_summary[1, 2]
    intercept_ci_lower <- tmp_conf[1,1]
    intercept_ci_upper <- tmp_conf[1,2]
    
    effect_est <- tmp_summary[2, 1]
    effect_se  <- tmp_summary[2, 2]
    effect_ci_lower <- tmp_conf[2,1]
    effect_ci_upper <- tmp_conf[2,2]
    
    residual_se <- summary(fit)$sigma
    output <- tibble("p_value" = p_value,
                     "intercept_est" = intercept_est,
                     "intercept_se" = intercept_se, 
                     "intercept_ci_lower" = intercept_ci_lower,
                     "intercept_ci_upper" = intercept_ci_upper,
                     "effect_est" = effect_est, 
                     "effect_se" = effect_se,
                     "effect_ci_lower" = effect_ci_lower, 
                     "effect_ci_upper" = effect_ci_upper,
                     "residual_se" = residual_se,
                     "phi_est" = NA,
                     "phi_se" = NA,
                     "phi_ci_lower" = NA,
                     "phi_ci_upper" = NA,
                     "method" = "lrlm")
  } else {
    output <- tibble("p_value" = NA,
                     "intercept_est" = NA,
                     "intercept_se" = NA, 
                     "intercept_ci_lower" = NA,
                     "intercept_ci_upper" = NA,
                     "effect_est" = NA, 
                     "effect_se" = NA,
                     "effect_ci_lower" = NA, 
                     "effect_ci_upper" = NA,
                     "residual_se" = NA, 
                     "phi_est" = NA,
                     "phi_se" = NA,
                     "phi_ci_lower" = NA,
                     "phi_ci_upper" = NA,
                     "method" = "lrlm")
  }
  return(output)
}  
.fit_two_group_bb_lr_corncob <- function(data, nominal_alpha = 0.05){
  # suppressed a warning message regarding a dependency. no impact on run. 
  data <- data %>% filter(., (num_count + denom_count) != 0)
  fit_null <- tryCatch({
    pkgcond::suppress_warnings({
      fit_null <- bbdml(formula = cbind(num_count, denom_count) ~ 1, 
                        phi.formula = ~ 1, data = data,  
                        inits = matrix(data = c(0,0), 
                                       nrow = 1, ncol = 2),
                        link = "logit", phi.link = "log")},
      pattern = "(will be removed from 'brglm2' at version 0.8.)")}, 
                       error = function(e){return(NA)})
  fit_alt <- tryCatch({
    pkgcond::suppress_warnings({
      fit_alt  <- bbdml(formula = cbind(num_count, denom_count) ~ treat, 
                        phi.formula = ~ 1, data = data,
                        inits = matrix(data = c(0,0,0), 
                                       nrow = 1, ncol = 3),
                        link = "logit", phi.link = "logit")},
      pattern = "(will be removed from 'brglm2' at version 0.8.)")}, 
    error = function(e){return(NA)})
  
  p_value <- tryCatch({p_value <- lrtest(mod_null = fit_null, mod = fit_alt)}, 
                      error = function(e){return(NA)})
  # if model fitting & LR successful, make outcome table; otherwise make NA 
  # table. 
  # P-value always exists as either NA or scalar value. 
  # Boolean thus always functional
  if (!(identical(fit_alt, NA)) && !is.na(p_value)) {
    # Estimate assumed to be on the logit scale, as is the standard error.
    # Asymptotic normality wald-approximation to CI
    # Alternative: bootstrap confidence intervals
    intercept_est <- summary(fit_alt)$coefficients[1,1]
    intercept_se  <- summary(fit_alt)$coefficients[1,2]
    intercept_ci_lower <- 
      intercept_est - (qnorm(nominal_alpha/2, lower.tail = FALSE) * intercept_se)
    intercept_ci_upper <- 
      intercept_est + (qnorm(nominal_alpha/2, lower.tail = FALSE) * intercept_se)
    
    effect_est <- summary(fit_alt)$coefficients[2,1]
    effect_se  <- summary(fit_alt)$coefficients[2,2]
    effect_ci_lower <- 
      effect_est - (qnorm(nominal_alpha/2, lower.tail = FALSE) * effect_se)
    effect_ci_upper <- 
      effect_est + (qnorm(nominal_alpha/2, lower.tail = FALSE) * effect_se)
    
    phi_est <- summary(fit_alt)$coefficients[3,1]
    phi_se  <- summary(fit_alt)$coefficients[3,2]
    phi_ci_lower <- 
      phi_est - (qnorm(nominal_alpha/2, lower.tail = FALSE) * phi_se)
    phi_ci_upper <- 
      phi_est + (qnorm(nominal_alpha/2, lower.tail = FALSE) * phi_se)
    
    output <- tibble("p_value" = p_value,
                     "intercept_est" = intercept_est,
                     "intercept_se" = intercept_se, 
                     "intercept_ci_lower" = intercept_ci_lower,
                     "intercept_ci_upper" = intercept_ci_upper,
                     "effect_est" = effect_est, 
                     "effect_se" = effect_se,
                     "effect_ci_lower" = effect_ci_lower, 
                     "effect_ci_upper" = effect_ci_upper,
                     "residual_se" = NA, 
                     "phi_est" = phi_est,
                     "phi_se" = phi_se,
                     "phi_ci_lower" = phi_ci_lower,
                     "phi_ci_upper" = phi_ci_upper,
                     "method" = "bbglm")
  } else {
    output <- tibble("p_value" = NA,
                     "intercept_est" = NA,
                     "intercept_se" = NA, 
                     "intercept_ci_lower" = NA,
                     "intercept_ci_upper" = NA,
                     "effect_est" = NA, 
                     "effect_se" = NA,
                     "effect_ci_lower" = NA, 
                     "effect_ci_upper" = NA,
                     "residual_se" = NA, 
                     "phi_est" = NA,
                     "phi_se" = NA,
                     "phi_ci_lower" = NA,
                     "phi_ci_upper" = NA,
                     "method" = "bbglm")
  }
  return(output)
}
.perform_two_group_wilcox <- function(data){
  data <- .perform_log_ratio_transformation(data)
  control <- filter(data, treat == 0) %>% pull(log_ratio)
  treat <- filter(data, treat == 1) %>% pull(log_ratio)
  tryCatch({fit <- wilcox.test(x = control, y = treat, 
                               alternative = "two.sided", correct = TRUE)}, 
           error = function(e){return(NA)})
  # warnings are acceptable; program goes for approximate p-val in case of ties
  p_value <- tryCatch({p_value <- fit$p.value}, 
                      error = function(e){return(NA)})
  # NA columns are created for bind_rows compatibility across methods
  output <- tibble("p_value" = p_value,
                   "intercept_est" = NA,
                   "intercept_se" = NA, 
                   "intercept_ci_lower" = NA,
                   "intercept_ci_upper" = NA,
                   "effect_est" = NA, 
                   "effect_se" = NA,
                   "effect_ci_lower" = NA, 
                   "effect_ci_upper" = NA,
                   "residual_se" = NA, 
                   "phi_est" = NA,
                   "phi_se" = NA,
                   "phi_ci_lower" = NA,
                   "phi_ci_upper" = NA,
                   "method" = "wilcox")
  return(output)
}
.perform_two_group_perm <- function(data){
  # Permutation test applied on log-ratio of components. 
  # Function from coin package.
  data <- .perform_log_ratio_transformation(data)
  tryCatch({fit <- coin::independence_test(formula = log_ratio ~ treat, data = data,
                                     distribution = 
                                       coin::approximate(nresample = 10000))}, 
           error = function(e){return(NA)})
  p_value <- tryCatch({p_value <- as.numeric(coin::pvalue(fit))}, 
                      error = function(e){return(NA)})
  # NA columns are created for bind_rows compatibility across methods
  output <- tibble("p_value" = p_value,
                   "intercept_est" = NA,
                   "intercept_se" = NA, 
                   "intercept_ci_lower" = NA,
                   "intercept_ci_upper" = NA,
                   "effect_est" = NA, 
                   "effect_se" = NA,
                   "effect_ci_lower" = NA, 
                   "effect_ci_upper" = NA,
                   "residual_se" = NA, 
                   "phi_est" = NA,
                   "phi_se" = NA,
                   "phi_ci_lower" = NA,
                   "phi_ci_upper" = NA,
                   "method" = "perm")
  return(output)
}

if (FALSE){
  library(corncob)
  # Example Run
  scenario <- list()
  scenario$alpha_control <- 0.1575591
  scenario$beta_control <- 1.044567
  scenario$n_samples <- 40
  scenario$mu_new <- 0.1310671
  scenario$disp_new <- 0.8318595
  scenario$alpha_new <- 0.1575591
  scenario$beta_new <- 1.044567
  scenario$n_sim <- 1
  scenario$method <- "wilcox"
  scenario$method <- "perm"
  scenario$method <- "bbglm"
  scenario$method <- "lrlm"
  scenario$method <- "lrlm_2"
  totals <- rep(2000, times = scenario$n_samples)
  set.seed(10)
  out <- .run_sim_para(alpha_control = scenario$alpha_control, 
                       beta_control  = scenario$beta_control, 
                       alpha_new = scenario$alpha_new, 
                       beta_new = scenario$beta_new,
                       n_samples  = scenario$n_samples,
                       totals = totals,
                       method = scenario$method)
  print(out, width = Inf)
  saved_out <- out
  rm(list = c("out"))
  
  # Re-running inner function part with saved seed to check reproducibility
  .Random.seed <- saved_out$seed[[1]]
  alpha_control = scenario$alpha_control
  beta_control  = scenario$beta_control
  alpha_new = scenario$alpha_new
  beta_new = scenario$beta_new
  n_samples  = scenario$n_samples
  totals = totals
  method = scenario$method
  # Resample Totals
  resampled_totals <- sample(totals, size = n_samples, replace = TRUE)
  # Generate Dataset
  data <- .make_beta_bin_data(alpha_control, beta_control, alpha_new, beta_new,
                              n_samples, totals = resampled_totals)
  if (method == "lrlm") {
    out <- .fit_two_group_log_ratio_lm(data)
  } else if (method == "bbglm") {
    out <- .fit_two_group_bb_lr_corncob(data)
  } else {stop("No valid method specified.")}
  
  print(out, width = Inf)
  print(saved_out, width = Inf)
  
  #' Example is fully reproducible given the input data and the seed. 
}
  