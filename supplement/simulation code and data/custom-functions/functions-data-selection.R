.generate_pairwise_table_reference_category <- function(dlong, reference){
  # Creates a table with pairwise ratio data 
  output <- list()
  n_samples <- length(unique(dlong$sample_id))
  
  # Make list of non-reference taxa to loop over
  taxa <- as.character(unique(dlong$taxon)) 
  taxa <- taxa[ -which(taxa == reference)]
  
  for (i in 1:length(taxa)) {
    pair <-  c(taxa[i], reference)
    tmp1 <- dlong %>% filter(., taxon == pair[1])
    tmp2 <- dlong %>% filter(., taxon == pair[2])
    
    # tmp1 taxon is more sparse -> numerator
    tmp1 <- tmp1 %>% 
      mutate(numerator = taxon, num_count = count) %>%
      select(sample_id, numerator, num_count)
    tmp2 <- tmp2 %>% 
      mutate(denominator = taxon, denom_count = count) %>%
      select(-taxon, -count)
    pair_table <- inner_join(tmp1, tmp2, by = "sample_id") %>%
      mutate(ratio_id = as.integer(i))
    
    # Add table to list
    output[[ i]] <- pair_table %>% select(sample_id, numerator, denominator, 
                                          num_count, denom_count, everything())
  }
  output
  output <- bind_rows(output, .id = NULL) %>%
    mutate(., ratio_id = as_factor(ratio_id)) %>%
    mutate(., taxon_pair = as_factor(paste(numerator, "-to-", 
                                           denominator, sep = "")))
  return(output)
}
.filter_by_prevalence <- function(dlong, prevalence_threshold){
  # Filtered such that only taxa with a zero percentage below or equal to the
  # threshold are kept in the table. Prevalence computed over entire dataset.
  # Assumes a long format data frame with a taxon identifier column and a count
  # column.
  dlong <- mutate(dlong, zero = if_else(count == 0, TRUE, FALSE)) %>%
    group_by(., taxon) %>%
    mutate(., zero_percentage = sum(zero) / n()) %>%
    ungroup() %>%
    filter(., zero_percentage <= prevalence_threshold)
  return(dlong)
}