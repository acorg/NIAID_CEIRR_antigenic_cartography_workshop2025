.FULLMAP = Racmacs::read.acmap(here::here("data", "maps", "01_sars_cov_2_map_square.ace"))

makeCheckMap = function(serum_limit){
  f = function(
    serum_choices,
    triangulate_function,
    noisy_boostrap_function
  ){
    serum_choices = serum_choices[!duplicated(serum_choices)]
    
    if (!all(serum_choices %in% srNames(.FULLMAP))){
      stop("Some of the sera you specified are are not in the map: ", 
           paste0(serum_choices[!serum_choices %in% srNames(.FULLMAP)], collapse = ", "))
    }
    
    if (length(serum_choices) > serum_limit){
      stop("Too many sera! You have ", length(serum_choices), ", you should have <=", serum_limit)
    } else if (length(serum_choices) < serum_limit){
      warning("You could have more sera if you'd like... you have ", length(serum_choices), ", you can have up to ", serum_limit)
    }
    
    acmap_sr = removeSera(.FULLMAP, which(!srNames(.FULLMAP) %in% serum_choices))
    acmap_sr = removeOptimizations(acmap_sr)
    
    acmap_sr = optimizeMap(acmap_sr, 2, 1000) %>%
      realignMap(.FULLMAP)
    
    par(mfrow = c(1, 3), mai = rep(0.2, 4))
    plot(acmap_sr, margins = NULL)
    plot(triangulate_function(acmap_sr), margins = NULL)
    plot(noisy_boostrap_function(acmap_sr), margins = NULL)
  }
}

checkMap15 = makeCheckMap(15)
checkMap8 = makeCheckMap(8)

