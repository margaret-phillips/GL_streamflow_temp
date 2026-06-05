"This is the main script that will run functions, organize outputs and create
necessary file directories."

source("streamflow_signatures.R")

#CREATE FILE DIRECTORIES AND FOLDERS



#Call up scripts with downloading and signature functions

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



#save output dataframe