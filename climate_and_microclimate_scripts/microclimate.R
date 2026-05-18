library(ncdf4)
library(dplyr)
library(lubridate)
library(NicheMapR)


process_era5_site_with_combinations <- function(site_name, site_lon, site_lat, site_elev, year, 
                                                target_combinations = NULL,
                                                nc_dir = "era_data/nc_files",
                                                output_dir = "processed_era",
                                                local_temp_dir = NULL) {
  
  # Default combinations if not provided
  if (is.null(target_combinations)) {
    target_combinations <- data.frame(
      shade = c(0.0, 0.0, 0.8, 0.8, 0.3),
      height = c(0.01, 0.15, 0.01, 0.15, 0.03)
    )
  }
  
  # CHANGE 1: Identify unique heights that need unshaded runs
  unique_heights <- unique(target_combinations$height)
  cat(sprintf("Unique heights that need unshaded air temperature: %s\n", 
              paste(unique_heights, collapse = ", ")))
  
  # Use a local temporary directory if shared drive has permission issues
  if (is.null(local_temp_dir)) {
    local_temp_dir <- tempdir()
  }
  cat("Using local temp directory:", local_temp_dir, "\n")
  
  cat(sprintf("Processing site %s for year %d with %d combinations\n", 
              site_name, year, nrow(target_combinations)))
  
  # Define file paths
  file_prefix <- "era5_summer"
  
  # Try to find renamed NetCDF file
  renamed_nc_file <- file.path(nc_dir, paste0("renamed_", file_prefix, "_", year, ".nc"))
  orig_nc_file <- file.path(nc_dir, paste0(file_prefix, "_", year, ".nc"))
  
  # Check which file exists
  if (file.exists(renamed_nc_file)) {
    input_nc_file <- renamed_nc_file
    is_renamed <- TRUE
    cat("Using renamed NetCDF file:", renamed_nc_file, "\n")
  } else if (file.exists(orig_nc_file)) {
    input_nc_file <- orig_nc_file
    is_renamed <- FALSE
    cat("Using original NetCDF file:", orig_nc_file, "\n")
  } else {
    stop("Neither original nor renamed NetCDF file found!")
  }
  
  # Create temp directory for processing
  site_dir <- file.path(local_temp_dir, site_name)
  if (!dir.exists(site_dir)) {
    dir.create(site_dir, recursive = TRUE)
  }
  
  # Create final output directory
  final_site_dir <- file.path(output_dir, site_name)
  if (!dir.exists(final_site_dir)) {
    dir.create(final_site_dir, recursive = TRUE)
  }
  
  # Apply inverse distance weighting to create a local version of the NetCDF file
  cat("Extracting data for", site_name, "at", site_lon, ",", site_lat, "...\n")
  
  nc <- nc_open(input_nc_file)
  
  # Get dimensions
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  
  # Find time dimension
  time_dim_name <- NULL
  for (dim_name in names(nc$dim)) {
    if (grepl("time|valid_time", dim_name, ignore.case = TRUE)) {
      time_dim_name <- dim_name
      cat("Found time dimension:", time_dim_name, "\n")
      break
    }
  }
  
  if (is.null(time_dim_name)) {
    nc_close(nc)
    stop("Could not find time dimension in NetCDF file")
  }
  
  time <- ncvar_get(nc, time_dim_name)
  
  
  
  
  # Add this right after: nc <- nc_open(input_nc_file)
  
  cat("=== NETCDF DIAGNOSTIC ===\n")
  
  # Check file structure
  cat("NetCDF file:", input_nc_file, "\n")
  cat("Variables:", names(nc$var), "\n")
  cat("Dimensions:", names(nc$dim), "\n")
  
  # Check dimension order for first variable
  first_var <- names(nc$var)[1]
  var_dims <- nc$var[[first_var]]$dim
  dim_names <- sapply(var_dims, function(d) d$name)
  cat("First variable dimension order:", paste(dim_names, collapse=", "), "\n")
  
  # Check coordinate values
  cat("Longitude range:", range(lon), "\n")
  cat("Latitude range:", range(lat), "\n")
  cat("Time range:", range(time), "\n")
  
  # Test reading a sample value from the target location
  test_lon_idx <- which.min(abs(lon - site_lon))
  test_lat_idx <- which.min(abs(lat - site_lat))
  cat("Target site coords:", site_lon, site_lat, "\n")
  cat("Closest grid coords:", lon[test_lon_idx], lat[test_lat_idx], "\n")
  
  # Try reading t2m at this location for first few time steps
  if ("t2m" %in% names(nc$var)) {
    cat("Testing t2m extraction...\n")
    
    # Method 1: Fixed order (longitude, latitude, time)
    t2m_fixed <- tryCatch({
      ncvar_get(nc, "t2m", start=c(test_lon_idx, test_lat_idx, 1), count=c(1,1,5))
    }, error = function(e) { cat("Fixed order failed:", e$message, "\n"); NA })
    
    # Method 2: Dynamic order (your current method)
    lon_pos <- which(dim_names == "longitude")
    lat_pos <- which(dim_names == "latitude")
    time_pos <- which(dim_names == time_dim_name)
    
    start_vec <- rep(1, length(dim_names))
    count_vec <- rep(1, length(dim_names))
    start_vec[lon_pos] <- test_lon_idx
    start_vec[lat_pos] <- test_lat_idx
    start_vec[time_pos] <- 1
    count_vec[lon_pos] <- 1
    count_vec[lat_pos] <- 1
    count_vec[time_pos] <- 5
    
    t2m_dynamic <- tryCatch({
      ncvar_get(nc, "t2m", start=start_vec, count=count_vec)
    }, error = function(e) { cat("Dynamic order failed:", e$message, "\n"); NA })
    
    cat("t2m values (fixed order):", t2m_fixed, "\n")
    cat("t2m values (dynamic order):", t2m_dynamic, "\n")
    cat("Expected temperature range: 270-290K\n")
    
    # Check if these are reasonable temperatures
    if (!is.na(t2m_fixed[1])) {
      cat("Fixed order temp in Celsius:", t2m_fixed[1] - 273.15, "\n")
    }
    if (!is.na(t2m_dynamic[1])) {
      cat("Dynamic order temp in Celsius:", t2m_dynamic[1] - 273.15, "\n")
    }
  }
  
  cat("=== END DIAGNOSTIC ===\n")
  
  
  
  # Find closest grid cell to target location
  lon_idx <- which.min(abs(lon - site_lon))
  lat_idx <- which.min(abs(lat - site_lat))
  
  closest_lon <- lon[lon_idx]
  closest_lat <- lat[lat_idx]
  
  cat("Target coordinates:", site_lon, ",", site_lat, "\n")
  cat("Closest grid cell:", closest_lon, ",", closest_lat, "\n")
  
  # Get variables needed for inverse distance weighting
  # Find the 4 closest grid cells
  dist_matrix <- expand.grid(lon_idx = 1:length(lon), lat_idx = 1:length(lat))
  dist_matrix$distance <- sqrt((lon[dist_matrix$lon_idx] - site_lon)^2 + 
                                 (lat[dist_matrix$lat_idx] - site_lat)^2)
  
  closest_cells <- dist_matrix %>% 
    arrange(distance) %>% 
    head(4)
  
  # Extract data for all relevant variables
  extracted_data <- list()
  var_names <- names(nc$var)
  
  # Determine the dimension order by checking the variable dimensions
  # This is the critical change to handle different dimensional orderings
  sample_var <- var_names[1]
  var_dims <- nc$var[[sample_var]]$dim
  dim_names <- sapply(var_dims, function(d) d$name)
  
  # Find positions of longitude, latitude, and time
  lon_pos <- which(dim_names == "longitude")
  lat_pos <- which(dim_names == "latitude")
  time_pos <- which(dim_names == time_dim_name)
  
  cat("NetCDF dimension order: ", paste(dim_names, collapse=", "), "\n")
  
  # Loop through all variables and extract data from the 4 closest grid cells
  for (var_name in var_names) {
    if (var_name %in% c("number", "expver")) next
    
    cat("Extracting variable:", var_name, "\n")
    
    # Initialize array to store weighted values
    weighted_values <- rep(0, length(time))
    weight_sum <- 0
    
    # Calculate inverse distance weighted values
    for (i in 1:nrow(closest_cells)) {
      lon_i <- closest_cells$lon_idx[i]
      lat_i <- closest_cells$lat_idx[i]
      dist_i <- closest_cells$distance[i]
      
      # Avoid division by zero
      if (dist_i == 0) {
        # If exact match, just use this cell
        weight_i <- 1
        weight_sum <- 1
        
        # Create dynamic start and count vectors based on dimension order
        start_vec <- rep(1, length(dim_names))
        count_vec <- rep(1, length(dim_names))
        
        # Set the indices for longitude, latitude, and time
        start_vec[lon_pos] <- lon_i
        start_vec[lat_pos] <- lat_i
        start_vec[time_pos] <- 1
        
        count_vec[lon_pos] <- 1
        count_vec[lat_pos] <- 1
        count_vec[time_pos] <- length(time)
        
        # Get data for this variable at this grid cell
        tryCatch({
          var_data <- ncvar_get(nc, var_name, 
                                start = start_vec, 
                                count = count_vec)
          
          weighted_values <- var_data
        }, error = function(e) {
          cat("Error getting variable", var_name, "at position", paste(start_vec, collapse=","), 
              "with count", paste(count_vec, collapse=","), "\n")
          cat("Error message:", e$message, "\n")
          # Return NA values if extraction fails
          return(rep(NA, length(time)))
        })
        break
      } else {
        # Calculate weight based on inverse distance
        weight_i <- 1 / dist_i
        weight_sum <- weight_sum + weight_i
        
        # Create dynamic start and count vectors based on dimension order
        start_vec <- rep(1, length(dim_names))
        count_vec <- rep(1, length(dim_names))
        
        # Set the indices for longitude, latitude, and time
        start_vec[lon_pos] <- lon_i
        start_vec[lat_pos] <- lat_i
        start_vec[time_pos] <- 1
        
        count_vec[lon_pos] <- 1
        count_vec[lat_pos] <- 1
        count_vec[time_pos] <- length(time)
        
        # Get data for this variable at this grid cell
        tryCatch({
          var_data <- ncvar_get(nc, var_name, 
                                start = start_vec, 
                                count = count_vec)
          
          # Add weighted contribution
          weighted_values <- weighted_values + weight_i * var_data
        }, error = function(e) {
          cat("Error getting variable", var_name, "at position", paste(start_vec, collapse=","), 
              "with count", paste(count_vec, collapse=","), "\n")
          cat("Error message:", e$message, "\n")
          # Continue with zero contribution
        })
      }
    }
    
    # Finalize weighted average if any weights were accumulated
    if (weight_sum > 0) {
      weighted_values <- weighted_values / weight_sum
    }
    
    # Store in our list
    extracted_data[[var_name]] <- weighted_values
  }
  
  # Create a copy of the original file with our IDW data
  idw_nc_file <- file.path(site_dir, paste0(file_prefix, "_", year, ".nc"))
  
  # First, copy the original file to this location
  file.copy(input_nc_file, idw_nc_file, overwrite = TRUE)
  
  # Open the file for modification
  nc_out <- nc_open(idw_nc_file, write = TRUE)
  
  # Replace data at the closest grid cell with our IDW values
  for (var_name in names(extracted_data)) {
    cat("Writing IDW data for:", var_name, "\n")
    
    # Create dynamic start and count vectors based on dimension order
    start_vec <- rep(1, length(dim_names))
    count_vec <- rep(1, length(dim_names))
    
    # Set the indices for longitude, latitude, and time
    start_vec[lon_pos] <- lon_idx
    start_vec[lat_pos] <- lat_idx
    start_vec[time_pos] <- 1
    
    count_vec[lon_pos] <- 1
    count_vec[lat_pos] <- 1
    count_vec[time_pos] <- length(time)
    
    # Write the data
    tryCatch({
      ncvar_put(nc_out, var_name, extracted_data[[var_name]], 
                start = start_vec, count = count_vec)
    }, error = function(e) {
      cat("Error writing variable", var_name, ":", e$message, "\n")
    })
  }
  
  # Close the NetCDF files
  nc_close(nc)
  nc_close(nc_out)
  
  # CHANGE 2: First, run micro_era5 for all unique heights with NO SHADE
  # Store unshaded air temperatures for each height
  unshaded_air_temps <- list()
  
  cat("\n=== Running unshaded micro_era5 for each unique height to get reference air temperatures ===\n")
  
  for (height in unique_heights) {
    cat(sprintf("\nGetting unshaded air temperature for height %.3f\n", height))
    
    # Define date range
    start_date <- paste0("01/06/", year)
    end_date <- paste0("31/08/", year)
    
    # Save current directory and timezone
    old_dir <- getwd()
    old_tz <- Sys.getenv("TZ")
    Sys.setenv(TZ = "America/Denver")
    
    # Change to the directory with the NetCDF file
    setwd(site_dir)
    
    # For micro_era5, use exactly the prefix it expects
    spatial_arg <- file_prefix
    
    # Run micro_era5 with NO SHADE (using minshade=0, maxshade=1 to avoid error)
    result <- tryCatch({
      micro_era5(
        loc = c(site_lon, site_lat),
        dstart = start_date,
        dfinish = end_date,
        Usrhyt = height,
        minshade = 0,     # NO SHADE
        maxshade = 1,     # Set to 1 to avoid error, but runshade=0 means it won't be used
        runshade = 0,     # This ensures only minshade is used
        RUF = 0.004,
        spatial = spatial_arg,
        weather.elev = site_elev,
        dem.res = 30,
        dem = NA,
        dem2 = NA
      )
    }, error = function(e) {
      cat("Error running micro_era5:", e$message, "\n")
      return(NULL)
    }, finally = {
      setwd(old_dir)
      Sys.setenv(TZ = old_tz)
    })
    
    if (!is.null(result)) {
      # Extract and store the unshaded air temperatures
      metout_data <- as.data.frame(result$metout)
      
      # Create datetime
      date_df <- data.frame(
        DOY = metout_data$DOY,
        TIME = metout_data$TIME,
        YEAR = year
      )
      date_df$date_only <- as.Date(paste(date_df$YEAR, date_df$DOY), format="%Y %j")
      date_df$hour <- date_df$TIME %/% 60
      date_df$minute <- date_df$TIME %% 60
      date_df$datetime <- as.POSIXct(
        paste(date_df$date_only, sprintf("%02d:%02d:00", date_df$hour, date_df$minute)),
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      )
      date_df$dtuse <- with_tz(date_df$datetime, tzone = "America/Denver")
      
      # Store the unshaded air temps with timestamps
      unshaded_air_temps[[as.character(height)]] <- data.frame(
        dtuse = date_df$dtuse,
        Tair_unshaded = metout_data$TALOC,
        stringsAsFactors = FALSE
      )
      
      cat(sprintf("Stored %d unshaded air temperature values for height %.3f\n", 
                  nrow(unshaded_air_temps[[as.character(height)]]), height))
    }
  }
  
  # CHANGE 3: Now process each combination as before, but replace air temperatures
  cat("\n=== Processing shade combinations with unshaded air temperatures ===\n")
  
  all_results <- list()
  
  for (i in 1:nrow(target_combinations)) {
    shade <- target_combinations$shade[i]
    height <- target_combinations$height[i]
    
    combo_id <- sprintf("shade%.1f_height%.3f", shade, height)
    cat(sprintf("\nProcessing combination %d/%d: %s\n", 
                i, nrow(target_combinations), combo_id))
    
    # Convert decimal shade to percentage for NicheMapR
    shade_pct <- shade * 100
    
    # Define date range
    start_date <- paste0("01/06/", year)
    end_date <- paste0("31/08/", year)
    
    # Save current directory and timezone
    old_dir <- getwd()
    old_tz <- Sys.getenv("TZ")
    Sys.setenv(TZ = "America/Denver")
    
    # Change to the directory with the NetCDF file
    setwd(site_dir)
    
    # For micro_era5, use exactly the prefix it expects
    spatial_arg <- file_prefix
    
    cat("Running micro_era5 with spatial arg:", spatial_arg, "\n")
    cat("Current working directory:", getwd(), "\n")
    cat("Files in directory:", paste(list.files(), collapse=", "), "\n")
    
    # Run micro_era5 for this combination (WITH SHADE if shade > 0)
    result <- tryCatch({
      micro_era5(
        loc = c(site_lon, site_lat),
        dstart = start_date,
        dfinish = end_date,
        Usrhyt = height,
        minshade = shade_pct,     # Shade percentage
        maxshade = 100,           # KEEP THIS AT 100 as in original code
        runshade = 0,
        RUF = 0.004,
        spatial = spatial_arg,
        weather.elev = site_elev,
        dem.res = 30,
        dem = NA,
        dem2 = NA, 
        ERR=2
      )
    }, error = function(e) {
      cat("Error running micro_era5:", e$message, "\n")
      return(NULL)
    }, finally = {
      setwd(old_dir)
      Sys.setenv(TZ = old_tz)
    })
    
    if (is.null(result)) {
      cat("Skipping this combination due to error.\n")
      next
    }
    
    # Process results for this combination
    cat("Processing results...\n")
    
    # Get metout and soil data
    metout_data <- as.data.frame(result$metout)
    soil_data <- as.data.frame(result$soil)
    
    # Apply shade factor to solar radiation - EXACTLY AS IN YOUR ORIGINAL CODE
    metout_data$SOLR <- metout_data$SOLR * (1 - shade)
    
    # Create a dataframe with time information
    date_df <- data.frame(
      DOY = metout_data$DOY,
      TIME = metout_data$TIME,
      YEAR = year,
      stringsAsFactors = FALSE
    )
    
    # Create dates
    date_df$date_only <- as.Date(paste(date_df$YEAR, date_df$DOY), format="%Y %j")
    
    # Create times
    date_df$hour <- date_df$TIME %/% 60
    date_df$minute <- date_df$TIME %% 60
    
    # Create full datetime
    date_df$datetime <- as.POSIXct(
      paste(
        date_df$date_only,
        sprintf("%02d:%02d:00", date_df$hour, date_df$minute)
      ),
      format = "%Y-%m-%d %H:%M:%S",
      tz = "UTC"
    )
    
    # Convert to Denver time
    date_df$dtuse <- with_tz(date_df$datetime, tzone = "America/Denver")
    
    # Add date components
    date_df$month <- month(date_df$dtuse)
    date_df$day <- day(date_df$dtuse)
    date_df$hour <- hour(date_df$dtuse)
    
    # Create date without year
    date_df$dt_noyr <- date_df$dtuse
    year(date_df$dt_noyr) <- 2024
    
    # Ensure row counts match
    min_rows <- min(nrow(date_df), nrow(metout_data), nrow(soil_data))
    date_df <- date_df[1:min_rows, ]
    metout_data <- metout_data[1:min_rows, ]
    soil_data <- soil_data[1:min_rows, ]
    
    # CHANGE 4: Replace air temperature with unshaded values
    # Get the unshaded air temps for this height
    unshaded_temps <- unshaded_air_temps[[as.character(height)]]
    
    if (is.null(unshaded_temps) || nrow(unshaded_temps) == 0) {
      cat(sprintf("WARNING: No unshaded temperatures found for height %.3f, using original temperatures\n", height))
      # Keep original temperatures if we don't have unshaded data
    } else {
      # Match by timestamp - being more careful about the merge
      temp_df <- data.frame(
        dtuse = date_df$dtuse, 
        row_num = 1:length(date_df$dtuse),
        stringsAsFactors = FALSE
      )
      
      matched_temps <- merge(
        temp_df,
        unshaded_temps,
        by = "dtuse",
        all.x = TRUE,
        sort = FALSE
      )
      matched_temps <- matched_temps[order(matched_temps$row_num), ]
      
      # Replace the air temperature if we have matching data
      if (nrow(matched_temps) > 0 && !all(is.na(matched_temps$Tair_unshaded))) {
        original_Tair <- metout_data$TALOC
        metout_data$TALOC <- matched_temps$Tair_unshaded
        
        # Handle any NAs by keeping original values
        na_indices <- is.na(metout_data$TALOC)
        if (any(na_indices)) {
          cat(sprintf("WARNING: %d NA values in unshaded temperatures, keeping original values\n", sum(na_indices)))
          metout_data$TALOC[na_indices] <- original_Tair[na_indices]
        }
        
        cat(sprintf("Replaced air temperatures: original range [%.2f, %.2f], new range [%.2f, %.2f]\n",
                    min(original_Tair, na.rm = TRUE), max(original_Tair, na.rm = TRUE),
                    min(metout_data$TALOC, na.rm = TRUE), max(metout_data$TALOC, na.rm = TRUE)))
      } else {
        cat("WARNING: Could not match timestamps, keeping original temperatures\n")
      }
    }
    
    # Create dataframe in simplified format - EXACTLY AS IN YOUR ORIGINAL CODE
    result_df <- data.frame(
      # Core datetime and location
      dtuse = date_df$dtuse,
      site = site_name,
      shade = shade,
      height = height,
      
      # Essential climate variables - USING REPLACED AIR TEMPERATURE
      Tair = metout_data$TALOC,      # Now using unshaded air temperature
      Tsoil = soil_data$D0cm,        # Soil temperature (still affected by shade)
      rad = metout_data$SOLR,        # Solar radiation (adjusted for shade)
      wind = if("VLOC" %in% names(metout_data)) {
        metout_data$VLOC
      } else {
        metout_data$VREF
      },
      
      # Time variables
      doy = metout_data$DOY,
      dt_noyr = date_df$dt_noyr,
      
      # Additional time fields
      month = date_df$month,
      day = date_df$day,
      hour = date_df$hour,
      year = year(date_df$dtuse),
      
      stringsAsFactors = FALSE
    )
    
    # Filter out 5/31 18:00 observations - EXACTLY AS IN YOUR ORIGINAL CODE
    result_df <- result_df %>%
      filter(!(month == 5 & day == 31 & hour == 18))
    
    # Add to results list
    all_results[[combo_id]] <- result_df
    
    cat(sprintf("Finished combination %s with %d observations\n", 
                combo_id, nrow(result_df)))
  }
  
  # Combine all results
  if (length(all_results) > 0) {
    combined_results <- bind_rows(all_results)
    
    # Sort by datetime and combination
    combined_results <- combined_results %>%
      arrange(dtuse, shade, height)
    
    # Add this after creating combined_results but before returning
    combined_results <- combined_results %>%
      group_by(dtuse, site, shade) %>%
      mutate(Tsoil = mean(Tsoil)) %>%
      ungroup()
    
    # Save RDS file
    simplified_file <- file.path(final_site_dir, 
                                 sprintf("simplified_%s_%d.rds", site_name, year))
    saveRDS(combined_results, simplified_file)
    cat(sprintf("\nSaved combined file with %d rows to %s\n", 
                nrow(combined_results), simplified_file))
    
    # Also save as CSV for easy viewing
    csv_file <- file.path(final_site_dir, 
                          sprintf("simplified_%s_%d.csv", site_name, year))
    write.csv(combined_results, csv_file, row.names = FALSE)
    cat(sprintf("Saved CSV version to %s\n", csv_file))
    
    # Return summary
    summary <- combined_results %>%
      group_by(site, shade, height) %>%
      summarize(
        n_obs = n(),
        min_date = min(dtuse),
        max_date = max(dtuse),
        .groups = "drop"
      )
    
    return(list(
      summary = summary,
      simplified = combined_results
    ))
  } else {
    cat("No successful combinations to save.\n")
    return(NULL)
  }
}


# Define site locations
sites <- data.frame(
  name = c("Eldo", "A1", "B1", "C1", "D1"),
  lat = c(39.944, 40.015, 40.023, 40.036, 40.059),
  lon = c(-105.262, -105.377, -105.430, -105.547, -105.617),
  elev = c(1740, 2195, 2591, 3048, 3515)
)

# Define the specific combinations you want
target_combinations <- data.frame(
  shade = c(0.0, 0.0, 0.8, 0.8, 0.3),
  height = c(0.01, 0.15, 0.01, 0.15, 0.03)
)

# Create a local temp directory that we definitely have write access to
local_temp <- "C:/Temp/era_processingJulia09152025"
if (!dir.exists(local_temp)) {
  dir.create(local_temp, recursive = TRUE)
}




# Example usage:
# For a single site and year
process_single_site <- function(site_name, year, nc_dir = "C:/Users/jmsmi/era5_workflow/merged",
                                output_dir = "C:/Users/jmsmi/era5_workflow/micro") {
  site_info <- sites[sites$name == site_name,]
  
  cat(sprintf("\n\n===== Processing %s for year %d =====\n\n", site_name, year))
  
  result <- process_era5_site_with_combinations(
    site_name = site_name,
    site_lon = site_info$lon,
    site_lat = site_info$lat,
    site_elev = site_info$elev,
    year = year,
    target_combinations = target_combinations,
    nc_dir = nc_dir,
    output_dir = output_dir,
    local_temp_dir = local_temp
  )
  
  return(result)
}

# For multiple sites and years
process_multiple_sites <- function(site_names = c("Eldo", "A1", "B1", "C1", "D1"),
                                   years = 2015:2019,
                                   nc_dir = "C:/Users/jmsmi/era5_workflow_new/merged_corrected",
                                   output_dir = "C:/Users/jmsmi/era5_workflow_new/micro2") {
  for (year in years) {
    for (site_name in site_names) {
      process_single_site(site_name, year, nc_dir, output_dir)
    }
  }
}

# Example calls:
#process_single_site("Eldo", 2019)
#process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 1951:1959) #c("Eldo", "A1", "B1", "C1", "D1"), 1950:1959

#process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 2021:2024)

process_multiple_sites(site_names = c("Eldo", "A1", "C1", "D1"), years = 2019)
process_multiple_sites(site_names = c("B1"), years = 2019)
process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 2020)
process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 2021:2024)
#process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 2015:2018)
process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 1950:1953)
process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 1958:1959)
process_multiple_sites(site_names = c("Eldo", "A1", "B1", "C1", "D1"), years = 2015:2018)
# 
# # Specify the current path and filename of the .nc file
# old_filepath <- "C:/Users/jmsmi/era5_workflow_new/merged_corrected/era5_summer_2018_corrected.nc"
# 
# # Specify the new path and filename for the .nc file
# new_filepath <- "C:/Users/jmsmi/era5_workflow_new/merged_corrected/era5_summer_2018.nc"
# 
# # Rename the file
# file.rename(from = old_filepath, to = new_filepath)
# 
# # Optional: Verify if the renaming was successful
# if (file.exists(new_filepath)) {
#   print(paste("File successfully renamed to:", new_filepath))
#   # You might want to also check if the old file no longer exists
#   if (!file.exists(old_filepath)) {
#     print("Old file removed.")
#   }
# } else {
#   print("File renaming failed.")
# }