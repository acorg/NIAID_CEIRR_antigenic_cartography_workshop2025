# plot the base map as r3js object
base_plot_data3js <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names,
                              add_border = TRUE, add_axis = TRUE, add_ag_labels = TRUE){
  
  # get coordinates for variants that should be plotted (highlighted ags)
  x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
  y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
  z_coords <- rep(0.02, length(highlighted_ags))
  ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
  ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
  ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
  labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  border_col <- "grey50"
  
  z_lims <- c(0,10)
  axis_at <- seq(z_lims[1], z_lims[2],2)
  
  # Setup plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid =  0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "log"
    ),
    show.axis = FALSE
  )
  
  # add z axis
  if(add_axis){
    
    axis_labels <- 2^axis_at*10
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
      cornerside = "f",
      size = 20,
      alignment = "right"
    )
  }
  
  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = c(0, 10),
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  # add variants to highlight
  data3js <- r3js::points3js(
    data3js,
    x          = x_coords,
    y          = y_coords,
    z          = z_coords,
    size       = ag_point_size,
    col        = ag_col,
    fill       = ag_fill,
    lwd        = 0.5,
    opacity    = 1,
    highlight  = list(col = "red"),
    label      = labels,
    toggle     = "Basepoints",
    depthWrite = FALSE,
    shape      = "circle filled"
  )
  
  if(add_ag_labels){
    text_x <- x_coords
    text_y <- c(y_coords - ag_point_size*0.2)
    
    data3js <- r3js::text3js(
      data3js,
      x          = text_x,
      y          = text_y,
      z          = z_coords,
      text       = labels,
      size       = c(rep(10*0.02, length(text_x))), 
      alignment  = "center"
    )
    
  }
  
  # thicker border of the base map
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[1]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    # y border
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[1]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[2], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    data3js <- r3js::box3js(
      data3js,
      col   = border_col
    )
    
  }
  
  
  return(data3js)
}


# plot landscapes from a list that contains the fits
plot_landscapes_from_list <- function(data3js, # base map as r3js object
                                      lndscp_fits, # landscape fits as list, named to match serum group in gmt_data
                                      map, # map for x,y coordinates of GMT points
                                      highlighted_ags, # which map antigens should be highlighted (GMT coordinates)
                                      lndscp_colors, # colors for the landscapes
                                      gmt_data = NULL, 
                                      show_gmts = TRUE, 
                                      show.individual.surfaces = FALSE, 
                                      options.individual.surfaces = list(opacity.surface.grid = 0.4,
                                                                         opacity.surface = 0.2, 
                                                                         col.surface = "grey70", 
                                                                         col.surface.grid = "grey70"), 
                                      show_landscapes = TRUE){
  
  
  # get x and y coords for gmt points
  x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
  y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
  coords <- cbind(x_coords, y_coords)
  coords <- coords[!is.na(x_coords),]
  
  
  # for each landscape fit object (serum group), do the plotting
  for (i in seq_along(lndscp_fits)) {
    
    # select sr group
    srg <- names(lndscp_fits)[i]
    lndscp_fit <- lndscp_fits[[i]]
    
    for (j in seq_len(nrow(coords))) {
      
      if(show_gmts){
        
        if(is.null(gmt_data)){
          warning("Error: You want to plot GMT but did not provide GMT data")
          return()
        }
        
        # select the gmt data for the serum group
        gmts <- gmt_data %>%
          filter(sr_group == srg)
        
        gmts <- gmts[match(rownames(coords), gmts$ag_name),]
        
        # plot GMTs
        data3js <- r3js::lines3js(
          data3js,
          x = rep(coords[j, 1], 2),
          y = rep(coords[j, 2], 2),
          z = c(0, gmts$gmt[j]),
          col = "grey50",
          toggle = sprintf("GMT, %s", srg),
          geometry = TRUE,
          opacity = 0.7,
          lwd = 0.2 
        )
        
        data3js <- r3js::points3js(
          data3js,
          x         = coords[j, 1],
          y         = coords[j, 2],
          z         = gmts$gmt[j],
          size      = 0.9,
          col  = lndscp_colors[srg, "Color"],
          toggle = sprintf("GMT, %s", srg),
          opacity   = 1 
          
        )
      }
      
    }
    
    
    if(show_landscapes){
      
      # Add GMT landscapes
      data3js <- lndscp3d_surface(
        data3js = data3js,
        object = lndscp_fit,
        crop2chull = FALSE,
        toggle = sprintf("GMT landscape, %s", srg),
        grid_spacing = 0.5,
        padding = 0.2,
        options = list(
          col.surface = lndscp_colors[srg, "Color"],
          opacity.surface = 0.8
        )
      )
      
      fit <- lndscp_fit
      
      # show individual landscapes
      if (show.individual.surfaces) {
        for (i in seq_len(nrow(fit$titers))) {
          individual_fit <- fit
          individual_fit$titers <- fit$titers[i, ]
          individual_fit$logtiters <- fit$logtiters[i, ]
          individual_fit$logtiters.upper <- individual_fit$logtiters.upper[i, 
          ]
          individual_fit$logtiters.lower <- individual_fit$logtiters.lower[i, 
          ]
          individual_fit$lessthans <- individual_fit$lessthans[i, 
          ]
          individual_fit$morethans <- individual_fit$morethans[i, 
          ]
          individual_fit$fitted.values <- NULL
          individual_fit$residuals <- NULL
          individual_fit$residuals.lessthan <- NULL
          individual_fit$residuals.morethan <- NULL
          if (!is.null(individual_fit$cone)) {
            individual_fit$cone$cone_coords <- individual_fit$cone$cone_coords[i, 
                                                                               , drop = F]
            individual_fit$cone$cone_heights <- individual_fit$cone$cone_heights[i]
          }
          data3js <- lndscp3d_surface(data3js = data3js, object = individual_fit, 
                                      crop2chull = FALSE, grid_spacing = 0.5, 
                                      padding = 0.1,
                                      options = options.individual.surfaces,
                                      toggle = "Individual landscapes")
        }
      }
      
    }
    
    
    
  }
  
  
  return(data3js)
}

# set orientation
set_r3js_orentation <- function(data3js, angle){
  r3js::r3js(
    data3js,
    rotation = angle$rotation,
    zoom = angle$zoom
  )
}

# bind landscape fit and calculated gmt
combine_landscape_and_calculated_gmt <- function(lndscp_fits, gmt_data){
  
  lndscp_gmts <- lapply(names(lndscp_fits), function(x){
    
    gmts <- data.frame(logtiter = lndscp_fits[[x]]$fitted.value,
                       ag_name = names(lndscp_fits[[x]]$fitted.value),
                       sr_group = x)
    return(gmts)
  })
  
  lndscp_gmts <- do.call(rbind, lndscp_gmts)
  
  
  ## compare lndscp gmts and claculated gmts
  comb_gmt <- rbind(lndscp_gmts %>%
                      mutate(Data = "Fitted Landscape GMT"),
                    gmt_data %>%
                      select(sr_group, ag_name, logtiter) %>%
                      unique() %>%
                      mutate(Data = "Calculated GMT"))
  
  
  return(comb_gmt)
}


# get residuals into long format
residuals_to_long <- function(residuals, values_name = "residuals"){
  
  ag_names <- colnames(residuals)
  as.data.frame(residuals) %>%
    rownames_to_column(var = "sr_name") %>%
    pivot_longer(cols = all_of(ag_names), names_to = "ag_name", values_to = values_name) -> residuals_long
  
  
  return(residuals_long)
}

# extract residuals from landscape fit function
combine_residuals <- function(fit, sr_group){
  
  
  residuals <- fit$residuals
  less_thans <- fit$residuals.lessthan
  more_thans <- fit$residuals.morethan
  
  long_res <- residuals_to_long(residuals, "residuals")
  long_less <- residuals_to_long(less_thans, "less_than")
  long_more <- residuals_to_long(more_thans, "more_than")
  
  serum_group <- sr_group
  
  comb <- long_res %>%
    left_join(long_less, by = c("ag_name", "sr_name")) %>%
    left_join(long_more, by = c("ag_name", "sr_name")) %>%
    mutate(measured = ifelse(is.na(residuals), ifelse(is.na(less_than), "more_than", "less_than"), "detectable"),
           residuals = ifelse(is.na(residuals), ifelse(is.na(less_than), more_than, less_than), residuals)) %>%
    select(!less_than:more_than) %>%
    mutate(sr_group = serum_group)
  
  return(comb)
}

# get root mean square error per variant
rmse_per_variant <- function(lndscp_fits){
  all_residuals <- lapply(names(lndscp_fits), function(x) combine_residuals(lndscp_fits[[x]], x))
  
  all_residuals <- do.call(rbind, all_residuals)
  
  all_residuals %>%
    filter(!is.na(residuals)) %>%
    group_by(sr_group) %>%
    summarize(ag_name = "Total", 
              rmse = sqrt(sum(residuals^2, na.rm = TRUE)/(length(residuals))), # Residual standard error
              residual_type = "All variants") -> total_ag
  
  all_residuals %>%
    filter(!is.na(residuals)) %>%
    group_by(ag_name, sr_group) %>%
    mutate(rmse = sqrt(sum(residuals^2, na.rm = TRUE)/(length(residuals))),
           residual_type = "By variant") %>%
    plyr::rbind.fill(., total_ag) -> ssr
  
  return(ssr)
}


return_idvl_fitted_value <- function(temp_fit, target_coords = NULL){
  
  cone_heights <- temp_fit$cone$cone_heights
  slope <- temp_fit$cone$cone_slope
  coords <- temp_fit$cone$cone_coords
  
  idvl_fitted <- list()
  for(x in 1:nrow(coords)){
    if(is.null(target_coords)){
      temp_dists <- dist(rbind(coords[x,], agCoords(temp_fit$acmap)), method = "euclidean")  
    } else {
      temp_dists <- dist(rbind(coords[x,], target_coords), method = "euclidean")  
    }
    
    target_dists <- as.matrix(temp_dists)[,1]
    
    idvl_fitted[[rownames(temp_fit$logtiters)[x]]] <- cone_heights[x] - slope*target_dists[2:length(target_dists)]
  }
  
  idvl_fitted <- do.call(rbind, idvl_fitted)
  
  return(idvl_fitted)
  
}


# calculate residuals for individual sera
calculate_idvl_residuals <- function(temp_fit){
  
  cone_heights <- temp_fit$cone$cone_heights
  slope <- temp_fit$cone$cone_slope
  coords <- temp_fit$cone$cone_coords
  
  idvl_residuals <- list()
  for(x in 1:nrow(coords)){
    temp_dists <- dist(rbind(coords[x,], agCoords(temp_fit$acmap)), method = "euclidean")
    target_dists <- as.matrix(temp_dists)[,1]
    
    fitted_vals <- cone_heights[x] - slope*target_dists[2:length(target_dists)]
    
    idvl_residuals[[rownames(temp_fit$logtiters)[x]]] <- fitted_vals - temp_fit$logtiters[x,]
  }
  
  idvl_residuals <- do.call(rbind, idvl_residuals)
  
  return(idvl_residuals)
}


