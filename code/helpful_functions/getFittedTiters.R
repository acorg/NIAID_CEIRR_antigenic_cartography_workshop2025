getFittedTiters <- function(map) {
  df <- tibble(
    "map_d" = as.numeric(as.vector(mapDistances(map))),
    "table_d" = as.vector(tableDistances(map)),
    "fitted_titer" = as.vector(
      sapply(
        1:ncol(mapDistances(map)), 
        function(x) {
          colBases(map)[x] - as.numeric(mapDistances(map)[, x])
        }
      )
    ),
    "fitted_titer_w_below" = as.vector(
      sapply(
        1:ncol(mapDistances(map)),
        function(x) {
          mapd <- as.numeric(gsub(">", "", mapDistances(map)[, x]))
          colBases(map)[x] - mapd
        }
      )
    ),
    "measured_titer" = as.vector(adjustedLogTiterTable(map)),
    sr = as.vector(matrix(
      srNames(map),
      nrow = Racmacs::numAntigens(map),
      ncol = Racmacs::numSera(map),
      byrow = T
    ))
  ) %>%
    mutate(
      below_detectable = grepl(">", table_d),
      table_d = as.numeric(gsub(">", "", table_d))
    ) %>%
    transmute(
      measured_titer = measured_titer,
      below_detectable = below_detectable,
      table_distance = table_d,
      map_distance = map_d,
      fitted_titer = fitted_titer_w_below,
      sr = sr
    )
  
  return(df)
}
