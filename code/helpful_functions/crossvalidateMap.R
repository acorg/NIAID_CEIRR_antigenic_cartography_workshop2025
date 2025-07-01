crossvalidateMap <- function(
    map,
    test_proportion = 0.1,
    number_of_optimizations = 100,
    number_of_replicates = 100,
    optimization_number = 1,
    options = list(dim_annealing = TRUE, ignore_disconnected = TRUE)
) {
  
  # Perform the CV testing
  cv_results <- Racmacs:::runDimensionTestMap(
    map = map,
    dimensions_to_test = mapDimensions(map, optimization_number),
    test_proportion = test_proportion,
    minimum_column_basis = minColBasis(map, optimization_number),
    fixed_column_bases = fixedColBases(map, optimization_number),
    number_of_optimizations = number_of_optimizations,
    replicates_per_dimension = number_of_replicates,
    options = options
  )
  
  # Summarise the results
  do.call(
    bind_rows,
    lapply(seq_along(cv_results$results), \(n) {
      tibble(
        measured_titer = cv_results$titers[cv_results$results[[n]]$test_indices],
        predicted_logtiter = cv_results$results[[n]]$predictions[[1]],
        titer_index = as.vector(cv_results$results[[n]]$test_indices),
        run = n
      )
    })
  ) %>% 
    mutate(
      ag_num = Racmacs:::agNumMatrix(map)[titer_index],
      sr_num = Racmacs:::srNumMatrix(map)[titer_index]
    ) %>%
    transmute(
      ag_name = agNames(map_op_10000)[ag_num],
      sr_name = srNames(map_op_10000)[sr_num],
      measured_logtiter = Racmacs:::log_titers(measured_titer, 0),
      predicted_logtiter = predicted_logtiter,
      measured_titer_type = c("numeric", "less than", "more than")[Racmacs:::titer_types_int(measured_titer)]
    )
}
