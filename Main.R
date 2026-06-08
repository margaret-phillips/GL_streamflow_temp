"This is the main script that will run functions, organize outputs and create
necessary file directories."

library(usethis)
library(dplyr)
library(tidyr)
library(purrr)
library(trend)

source("streamflow_signatures.R")
source("hydroclimate_signatures.R")



# Example: ignore folders and files
use_git_ignore(c("data/daymet", "notes.txt", "*.csv", ".Rhistory"))

#CREATE FILE DIRECTORIES AND FOLDERS

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




#Save the output dataframe with all signatures for site years



#Process the outputs and compute trends throughout the time period

#move this to a separate script and then source it..?


compute_trends<- function(df, time_col = "water_year") {
  
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
        summarise(value = mean(value), .groups = "drop") %>%  # or first(value)
        arrange(.data[[time_col]])
      
      
      sen <- trend::sens.slope(df2$value)
      mk  <- trend::mk.test(df2$value)
      
      tibble(
        sen_slope = sen$estimates,
        sen_conf_low = sen$conf.int[1], #maybe remove? too many variables
        sen_conf_high = sen$conf.int[2], #maybe remove? too many variables
        sen_pval = sen$p.value,
        mk_tau = cor(df2[[time_col]], df2$value, method = "kendall"),
        mk_pval = mk$p.value
      )
    }) %>%
    ungroup() %>%
    
    distinct(monitoring_location_id, signature, .keep_all = TRUE)
  
 
  results_wide <- results_long %>%
    pivot_wider(
      names_from = signature,
      values_from = c(
        sen_slope,
        sen_conf_low,
        sen_conf_high,
        sen_pval,
        mk_tau,
        mk_pval
      ),
      names_sep = "__"
    )

  results_final <- results_wide %>%
    left_join(n_years_df, by = "monitoring_location_id")
  
  return(results_final)
}


#save output dataframe



##--------------hydroclimate-------------------------------------------------####

hydroclimate_functions <- list(calculate_air_temp, calculate_prcp_seasonality, calculate_prcp_timing, 
                               calculate_rain_snow, calculate_spring_days)


hydroclimate_output_dfs <- lapply(
  hydroclimate_functions,
  function(f) f(daymet_annual, save_path = NULL)
)

hydroclimate_combined_df <- reduce(
  hydroclimate_output_dfs,
  full_join,
  by = c("monitoring_location_id", "water_year")
)


hydroclimate_trends<- compute_trends(hydroclimate_combined_df)
