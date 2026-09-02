"This is the main script that will run functions, organize outputs and create
necessary file directories."

library(usethis)
library(dplyr)
library(tidyr)
library(purrr)
library(trend)

source("streamflow_signatures.R")
source("hydroclimate_signatures.R")
source("stream_temperature_signatures.R")
source("trends.R")


# Example: ignore folders and files
use_git_ignore(c("data/daymet", "notes.txt", "*.csv", ".Rhistory"))

#CREATE FILE DIRECTORIES AND FOLDERS


##-----------streamflow------------------------------------------------------####

#Call up scripts with downloading and signature functions!



streamflow_functions <- list(calculate_percentiles, calculate_FDC, calculate_CV_IQD, 
                             calculate_q_seasonality, calculate_flow_pulses, calculate_flashiness,
                             calculate_frequency_lows, calculate_frequency_highs)


streamflow_output_dfs <- lapply(
  streamflow_functions,
  function(f) f(cleaned_dv_qDepth_annual, save_path = NULL)
)

streamflow_combined_df <- reduce(
  streamflow_output_dfs,
  full_join,
  by = c("monitoring_location_id", "water_year")
)

#check for missing years within temporal coverage at each site
#preserve missing years as NA gap to get accurate trends
library(dplyr)
library(tidyr)

missing_wy_check <- streamflow_combined_df |>
  distinct(monitoring_location_id, water_year) |>
  group_by(monitoring_location_id) |>
  complete(water_year = seq(min(water_year), max(water_year))) |>
  filter(is.na(monitoring_location_id)) |>
  summarise(
    missing_years = list(water_year),
    n_missing = n(),
    .groups = "drop"
  )

missing_wy_check <- streamflow_combined_df |>
  distinct(monitoring_location_id, water_year) |>
  group_by(monitoring_location_id) |>
  summarise(
    start_year= min(water_year),
    end_year= max(water_year),
    missing_years = list(
      setdiff(
        min(water_year):max(water_year),
        water_year
      )
    ),
    n_missing = lengths(missing_years),
    .groups = "drop"
  ) |>
  filter(n_missing > 0)

#probably want to limit to 60% completeness in the first and last decade 
site_years <- streamflow_combined_df |>
  distinct(monitoring_location_id, water_year)

years_to_remove <- site_years |>
  group_by(monitoring_location_id) |>
  group_modify(~{
    
    yrs <- sort(.x$water_year)
    
    first_decade <- min(yrs):(min(yrs) + 9)
    last_decade  <- (max(yrs) - 9):max(yrs)
    
    n_first <- sum(yrs %in% first_decade)
    n_last  <- sum(yrs %in% last_decade)
    
    tibble(
      water_year = c(
        if (n_first <= 5) intersect(yrs, first_decade),
        if (n_last <= 5) intersect(yrs, last_decade)
      )
    )
  }) |>
  ungroup() |>
  distinct()


streamflow_cleaned <- streamflow_combined_df |>
  anti_join(
    years_to_remove,
    by = c("monitoring_location_id", "water_year")
  )

#running trends on streamflow signatures
streamflow_trends<- compute_trends(streamflow_cleaned)
#save output dataframe



##--------------hydroclimate-------------------------------------------------####

hydroclimate_functions <- list(calculate_air_temp, calculate_prcp_seasonality, calculate_prcp_timing, 
                               calculate_rain_snow, calculate_spring_days, calculate_seasonal_prcp, calculate_seasonal_airtemp)


hydroclimate_output_dfs <- lapply(
  hydroclimate_functions,
  function(f) f(clean_daymet_annual, save_path = NULL)
)

hydroclimate_combined_df <- reduce(
  hydroclimate_output_dfs,
  full_join,
  by = c("monitoring_location_id", "water_year")
)


hydroclimate_trends<- compute_trends(hydroclimate_combined_df)


##---------------------stream temperature-------------------------------------####

streamtemp_functions <- list(calculate_tw_20C_dur, calculate_tw_avg, calculate_tw_grow_days, 
                               calculate_tw_fall, calculate_tw_spring, calculate_tw_summer)


streamtemp_output_dfs <- lapply(
  streamtemp_functions,
  function(f) f(cleaned_dv_tw, save_path = NULL)
)



streamtemp_combined_df <- reduce(
  streamtemp_output_dfs,
  full_join,
  by = c("monitoring_location_id", "water_year")
)

missing_wy_check <- streamtemp_combined_df |>
  distinct(monitoring_location_id, water_year) |>
  group_by(monitoring_location_id) |>
  summarise(
    missing_years = list(
      setdiff(
        min(water_year):max(water_year),
        water_year
      )
    ),
    n_missing = lengths(missing_years),
    .groups = "drop"
  ) |>
  filter(n_missing > 0)


site_years <- streamtemp_combined_df |>
  distinct(monitoring_location_id, water_year)

years_to_remove <- site_years |>
  group_by(monitoring_location_id) |>
  group_modify(~{
    
    yrs <- sort(.x$water_year)
    
    first_decade <- min(yrs):(min(yrs) + 9)
    last_decade  <- (max(yrs) - 9):max(yrs)
    
    n_first <- sum(yrs %in% first_decade)
    n_last  <- sum(yrs %in% last_decade)
    
    tibble(
      water_year = c(
        if (n_first <= 5) intersect(yrs, first_decade),
        if (n_last <= 5) intersect(yrs, last_decade)
      )
    )
  }) |>
  ungroup() |>
  distinct()


streamtemp_cleaned <- streamtemp_combined_df |>
  anti_join(
    years_to_remove,
    by = c("monitoring_location_id", "water_year")
  )


streamtemp_trends<- compute_trends(streamtemp_cleaned)

##---------------------------seasonal version--------------------------------####



#Call up scripts with downloading and signature functions!

#can functionalize this to use for hydroclimate, streamflow, and temp

streamflow_functions <- list(calculate_percentiles,
                              calculate_flow_pulses, calculate_flashiness,
                             calculate_frequency_lows, calculate_frequency_highs)


fall_streamflow_output_dfs <- lapply(
  streamflow_functions,
  function(f) f(fall_dv_q, save_path = NULL)
)

fall_streamflow_combined_df <- reduce(
  fall_streamflow_output_dfs,
  full_join,
  by = c("monitoring_location_id", "water_year")
)

#check for missing years within temporal coverage at each site

library(dplyr)
library(tidyr)


#probably want to limit to 60% completeness in the first and last decade 


fall_streamflow_cleaned <- fall_streamflow_combined_df |>
  anti_join(
    years_to_remove,
    by = c("monitoring_location_id", "water_year")
  )

#running trends on streamflow signatures
fall_streamflow_trends<- compute_trends(fall_streamflow_cleaned)
#save output dataframe


