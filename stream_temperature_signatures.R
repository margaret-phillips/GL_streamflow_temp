"This script uses water temperature data to calculate stream temperature signatures.
Signatures are computed based on corresponding data availability (e.g. maximum stream
temperature is computed when summer months are available, but average temperature is 
only computed if all months are available)."

##------------magnitude------------------------------------------------------####

#avg max and min stream temp

calculate_tw_avg<- function(tw_data, save_path= NULL){
  
  time_range<- 1:12 #require all months to compute average
  
  output_df<- tw_data %>% 
    filter(month %in% time_range) %>% 
    group_by(monitoring_location_id, water_year) %>% 
    summarise(
      has_full_range = all(time_range %in% month),
      complete_months= all(table(month)>= 28),
      tw_avg= mean(tw, na.rm= TRUE),
      .groups= "drop") %>% 
    filter(has_full_range & complete_months) %>% 
    select(-has_full_range, -complete_months)
 
   if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}


calculate_tw_avg(cleaned_dv_tw, save_path= NULL)

calculate_tw_grow_days <- function(tw_data, save_path = NULL) {
  Tgrow <- 5
  time_range <- 2:11
  
  output_df <- tw_data %>% 
    filter(month %in% time_range) %>% 
    group_by(monitoring_location_id, water_year) %>% 
    summarise(
      has_full_range = all(time_range %in% month),
      complete_months = all(table(month) >= 28),
      growing_deg = sum(pmax(tw - Tgrow, 0)),
      .groups = "drop"
    ) %>% 
    filter(has_full_range & complete_months) %>% 
    select(-has_full_range, -complete_months)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

calculate_tw_grow_days(cleaned_dv_tw, save_path= NULL)

##------------------timing---------------------------------------------------####

#first winter/spring day (using 5 day running avg) above 5C
#not going to calculate this for now since I have growing degree days already

##-----------------------duration--------------------------------------------####


#max consecutive days above 20C

calculate_tw_20C_dur<- function(tw_data, save_path= NULL){
  
  time_range<- 5:9
  T_thresh<- 20 #choosing 20 C as a threshold temp, but can be adjusted
  
  output_df<- tw_data %>% 
    filter(month %in% time_range) %>% 
    group_by(monitoring_location_id, water_year) %>% 
    arrange(time) %>% 
    mutate(
      heat_event = coalesce(tw > T_thresh, FALSE),
      heat_event_num = cumsum(heat_event & !lag(heat_event, default = FALSE)),
      ) %>% 
    summarise(
      has_full_range = all(time_range %in% month),
      complete_months = all(table(month) >= 28),
      total_heat_events = max(heat_event_num),
      max_dur_abv20 = {
        r <- rle(heat_event)
        if (any(r$values)) max(r$lengths[r$values]) else 0
      },
      .groups= "drop"
    ) %>% 
    filter(has_full_range & complete_months) %>% 
    select(-has_full_range, -complete_months)
  
    
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

calculate_tw_20C_dur(cleaned_dv_tw)

##-----------------------------variability-----------------------------------####

#rate of change for spring and Fall (deg C per day)-- (Chu et al., 2010)

calculate_tw_roc_spring <- function(tw_data, save_path = NULL) {
  
  time_range <- 4:6
  
  output_df <- tw_data %>%
    filter(month %in% time_range) %>%
    arrange(monitoring_location_id, water_year, time) %>%
    group_by(monitoring_location_id, water_year) %>%
    mutate(
      delta_tw   = tw - lag(tw),
      delta_time = as.numeric(difftime(time, lag(time), units = "days"))
    ) %>%
    summarise(
      has_full_range = all(time_range %in% month),
      complete_months = all(table(month) >= 28),
      
      # only use true daily steps
      min_spring_tw_roc = min(delta_tw[delta_time == 1], na.rm = TRUE),
      spring_tw_roc     = mean(delta_tw[delta_time == 1], na.rm = TRUE),
      max_spring_tw_roc = max(delta_tw[delta_time == 1], na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    filter(has_full_range & complete_months) %>%
    select(-has_full_range, -complete_months)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

calculate_tw_roc_spring(cleaned_dv_tw)



calculate_tw_roc_fall <- function(tw_data, save_path = NULL) {
  
  time_range <- 9:11
  
  output_df <- tw_data %>%
    filter(month %in% time_range) %>%
    arrange(monitoring_location_id, water_year, time) %>%
    group_by(monitoring_location_id, water_year) %>%
    mutate(
      delta_tw   = tw - lag(tw), #this was creating a problem bc the months cross a water year!
      delta_time = as.numeric(difftime(time, lag(time), units = "days"))
    ) %>%
    summarise(
      has_full_range = all(time_range %in% month),
      complete_months = all(table(month) >= 28),
      
      # only use true daily steps
      min_fall_tw_roc = min(delta_tw[delta_time == 1], na.rm = TRUE), #enforced one day only lags
      fall_tw_roc     = mean(delta_tw[delta_time == 1], na.rm = TRUE),
      max_fall_tw_roc = max(delta_tw[delta_time == 1], na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    filter(has_full_range & complete_months) %>%
    select(-has_full_range, -complete_months)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

