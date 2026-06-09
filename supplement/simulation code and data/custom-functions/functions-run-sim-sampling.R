.run_sim_resample <- function(n_samples, raw_ratio_data, method){
  tmp_seed <- .Random.seed
  # Generate Dataset through resampling
  data <- .sample_data(raw_ratio_data, n_samples)
  if (method == "lrlm") {
    out <- .fit_two_group_log_ratio_lm_resample(data)
  } else if (method == "bbglm") {
    out <- .fit_two_group_bb_lr_corncob_resample(data)
  } else if (method == "wilcox") {
      out <- .perform_two_group_wilcox(data)
  } else if (method == "perm") {
    out <- .perform_two_group_perm(data)
  } else {stop("No valid method specified.")}
  out$seed <- list(tmp_seed)
  return(out)
}
.sample_data <- function(ratio_data, n_samples){
  studies <- unique(ratio_data$study_id)
  n_studies <- length(studies)
  samples_per_study <- n_samples / n_studies
  tmp <- list()
  for (i in studies){
    # filter to selected study, randomly sample rows, assign treatment
    tmp[[i]] <- filter(ratio_data, study_id == i) %>%
      sample_n(., size = samples_per_study, replace = TRUE) %>%
      mutate(., treat = as.integer(rep(c(0,1), length.out = samples_per_study)))
  }
  dataset <- bind_rows(tmp)
  return(dataset)
}
.perform_log_ratio_transformation <- function(input_data){
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
                   "method" = "wilcox")
  return(output)
}
.perform_two_group_perm <- function(data){
  # Permutation test applied on log-ratio of components. 
  # Function from coin package.
  data <- .perform_log_ratio_transformation(data)
  tryCatch({fit <- coin::independence_test(formula = log_ratio ~ treat, 
                                           data = data,
                                           distribution = 
                                             coin::approximate(nresample = 
                                                                 10000))}, 
           error = function(e){return(NA)})
  p_value <- tryCatch({p_value <- as.numeric(coin::pvalue(fit))}, 
                      error = function(e){return(NA)})
  # NA columns are created for bind_rows compatibility across methods
  output <- tibble("p_value" = p_value,
                   "method" = "perm")
  return(output)
}
.fit_two_group_log_ratio_lm_resample <- function(data, alpha = 0.05){
  # Pseudo-count impute data
  data <- .perform_log_ratio_transformation(data)
  data <- mutate(data, study_id = as.factor(study_id))
  
  fit <- tryCatch({fit <- lm(formula = log_ratio ~ treat + study_id, data = data)},
                  error = function(e){return(NA)},  warning = function(w){NA})
  # if model successful, make outcome table; otherwise make NA table
  if (!identical(fit, NA)) {
    # Create summary objects for accessing fitted model components
    tmp_summary <- summary(fit)$coefficients
    p_value <- tmp_summary[2 , 4]
    output <- tibble("p_value" = p_value, "method" = "lrlm")
  } else {
    output <- tibble("p_value" = NA, "method" = "lrlm")
  }
  return(output)
}  
.fit_two_group_bb_lr_corncob_resample <- function(data, 
                                                  nominal_alpha = 0.05){
  data <- data %>% filter(., (num_count + denom_count) != 0)
  data <- mutate(data, study_id = as.factor(study_id))
  
  fit_null <- tryCatch({
    pkgcond::suppress_warnings({
      fit_null <- 
        bbdml(formula = cbind(num_count, denom_count) ~ 1 + study_id, 
              phi.formula = ~ 1, 
              data = data, inits = matrix(data = c(0,0,0), nrow = 1, ncol = 3),
              link = "logit", phi.link = "log")},
      pattern = "(will be removed from 'brglm2' at version 0.8.)")}, 
    error = function(e){return(NA)})
  
  fit_alt <- tryCatch({
    pkgcond::suppress_warnings({
      fit_alt  <- 
        bbdml(formula = cbind(num_count, denom_count) ~ treat + study_id, 
              phi.formula = ~ 1, data = data,
              inits = matrix(data = c(0,0,0,0), nrow = 1, ncol = 4),
              link = "logit", phi.link = "logit")},
      pattern = "(will be removed from 'brglm2' at version 0.8.)")}, 
    error = function(e){return(NA)})
  
  p_value <- tryCatch({p_value <- lrtest(mod_null = fit_null, mod = fit_alt)}, 
                      error = function(e){return(NA)})
  # P-value always exists as either NA or scalar value. 
  output <- tibble("p_value" = p_value, "method" = "bbglm")
  return(output)
}