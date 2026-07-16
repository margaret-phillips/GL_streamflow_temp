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


# Example: ignore folders and files
use_git_ignore(c("data/daymet", "notes.txt", "*.csv", ".Rhistory"))

#CREATE FILE DIRECTORIES AND FOLDERS

##-------trends---------------------------------------------------------------####

#move this to a separate script and then source it..?

compute_trends <- function(df, time_col = "water_year") {
  
  n_years_df <- df %>%
    group_by(monitoring_location_id) %>%
    summarise(
      n_years = n_distinct(.data[[time_col]]),
      .groups = "drop"
    )
  
  df_long <- df %>%
    pivot_longer(
      cols = -c(monitoring_location_id, all_of(time_col)),
      names_to = "signature",
      values_to = "value"
    )
  
  results_long <- df_long %>%
    group_by(monitoring_location_id, signature) %>%
    group_modify(~ {
      
      df2 <- .x %>%
        filter(!is.na(value)) %>%
        group_by(.data[[time_col]]) %>%
        summarise(
          value = mean(value),
          .groups = "drop"
        ) %>%
        arrange(.data[[time_col]])
      
      n_years_sig <- n_distinct(df2[[time_col]])
      mean_val <- mean(df2$value, na.rm = TRUE)
      
      #only compute trends if >= 10 years
      if (n_years_sig < 10) {
        return(tibble(
          sen_slope = NA_real_,
          sen_pval = NA_real_,
          mk_tau = NA_real_,
          mk_pval = NA_real_,
          mean_value = mean_val,
          has_trend = FALSE
        ))
      }
      
      sen <- trend::sens.slope(df2$value)
      mk  <- trend::mk.test(df2$value)
      
      tibble(
        sen_slope = sen$estimates,
        sen_pval = sen$p.value,
        mk_tau = cor(df2[[time_col]], df2$value, method = "kendall"),
        mk_pval = mk$p.value,
        mean_value = mean_val,
        has_trend = TRUE
      )
    }) %>%
    ungroup()
  
  
  valid_sites <- results_long %>%
    group_by(monitoring_location_id) %>%
    summarise(any_trend = any(has_trend), .groups = "drop") %>%
    filter(any_trend) %>%
    pull(monitoring_location_id)
  
  
  results_long <- results_long %>%
    filter(monitoring_location_id %in% valid_sites)
  
  
  results_wide <- results_long %>%
    select(-has_trend) %>%
    pivot_wider(
      names_from = signature,
      values_from = c(
        sen_slope,
        sen_pval,
        mk_tau,
        mk_pval,
        mean_value
      ),
      names_sep = "__"
    )
  
  results_final <- results_wide %>%
    left_join(n_years_df, by = "monitoring_location_id")
  
  return(results_final)
}


##-----------streamflow------------------------------------------------------####

#Call up scripts with downloading and signature functions!

#can functionalize this to use for hydroclimate, streamflow, and temp

streamflow_functions <- list(calculate_percentiles, calculate_FDC, calculate_CV_IQD, 
                             calculate_q_seasonality, calculate_flow_pulses, calculate_flashiness)


streamflow_output_dfs <- lapply(
  streamflow_functions,
  function(f) f(cleaned_dv_qDepth_annual, save_path = NULL)
)

streamflow_combined_df <- reduce(
  streamflow_output_dfs,
  full_join,
  by = c("monitoring_location_id", "water_year")
)


#running trends on streamflow signatures
streamflow_trends<- compute_trends(streamflow_combined_df)
#save output dataframe



##--------------hydroclimate-------------------------------------------------####

hydroclimate_functions <- list(calculate_air_temp, calculate_prcp_seasonality, calculate_prcp_timing, 
                               calculate_rain_snow, calculate_spring_days, calculate_seasonal_prcp, calculate_seasonal_airtemp)


hydroclimate_output_dfs <- lapply(
  hydroclimate_functions,
  function(f) f(daymet_annual_clean, save_path = NULL)
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

streamtemp_trends<- compute_trends(streamtemp_combined_df)
#not calling stream temp trends. only 5 sites have 25 years of at least on signature


