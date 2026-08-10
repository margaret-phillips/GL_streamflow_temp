"This script includes functions to calcualte hydroclimate signatures using daily air 
temperature, SWE, and precipitation data for each site in the daymet dataset. Signatures
are calculated at the annual timescale, with some applying only to specific time periods
within the year. The output signatures dataframes are then processed by a wrapper function
to compute trends, significance, and variability in signatures over the time period."


#required libraries:
library(tidyverse)
library(lubridate)

#read in daymet dataset
daymet_annual<- readRDS("data/daymet/daymet.rds")


clean_daymet_annual<- daymet_annual %>% 
  mutate(monitoring_location_id= paste0("USGS-", site_id),
         water_year = if_else(month >= 10, year + 1, year)) %>% 
  group_by(monitoring_location_id, water_year) %>%
  mutate(wy_doy= row_number()) %>% 
  filter(n() >= 364) #bc of how daymet wraps years and deals with leap years, a water year has btwn 364 and 366 days!

wy_check<- clean_daymet_annual %>% 
  group_by(monitoring_location_id) %>% 
  summarise(complete_yrs= n_distinct(water_year))#verifying that coverage spans complete wy time period!



#NEED TO ADD A FEW SIGNATURES FOR ANNUAL SCALE! and make sure annual daymet has required cols

##-----------air temperature---------------------------------------------------####

#This function calculates minimum, maximum, and average air temp
calculate_air_temp<- function(climate_data, save_path= NULL){
  output_df<- climate_data %>% 
    group_by(monitoring_location_id) %>% 
    mutate(tavg= (tmin + tmax)/2) %>% 
    group_by(monitoring_location_id, water_year) %>% 
    summarise(tmin= min(tmin),
              tmax= max(tmax),
              tavg= mean(tavg))
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}





#This function calculates the last spring day below -2.2C and number of days above 5C
calculate_spring_days<- function(climate_data, save_path= NULL){
  Tgrow<- 5
  Tdamage<- -2.2
  
  output_df<- climate_data %>% 
    group_by(monitoring_location_id, water_year) %>% 
    summarise(
      growing_deg = sum(pmax(tmin - Tgrow, 0), na.rm = TRUE), #cumulative degrees above 5C
      
      last_damage_day = {
        idx <- which(tmin < Tdamage & wy_doy >= 0 & wy_doy <= 180)
        if (length(idx) == 0) NA else max(wy_doy[idx])
      }, #spring doy that corresponds to last value below -2.2C...
      .groups = "drop"
      
    )
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}


##---------------------precipitation -------------------------------------------####

#this function calculates rain, snow, and melt signatures
calculate_rain_snow<- function(climate_data, save_path= NULL){
  swe_period<- list(11, 12, 1, 2, 3, 4)
  
  output_df <- climate_data %>%
    
    group_by(monitoring_location_id, water_year) %>%
    mutate(
      delta_swe = swe - lag(swe),
      acc_swe = ifelse(delta_swe > 0, delta_swe, 0),
      melt = ifelse(delta_swe < 0, -delta_swe, 0),
      rain = ifelse(prcp > 0 & delta_swe <= 0, prcp, 0),
      snow = ifelse(prcp > 0 & delta_swe > 0, prcp, 0),
      melt_event = coalesce(melt > 0, FALSE),
      melt_event_num = cumsum(melt_event & !lag(melt_event, default = FALSE))
    ) %>%
    
    
    mutate(
      total_rain = sum(rain, na.rm = TRUE),
      total_snow = sum(snow, na.rm = TRUE),
      total_prcp = sum(prcp, na.rm = TRUE),
      total_melt_events = max(melt_event_num, na.rm = TRUE),
      max_melt_daily = max(melt, na.rm = TRUE),
      swe_ratio = total_snow / total_prcp,
      ws_swe_ratio= 
        sum(snow[month %in% swe_period], na.rm = TRUE) /
        sum(prcp[month %in% swe_period], na.rm = TRUE)
    ) %>%
    
    group_by(monitoring_location_id, water_year, melt_event_num) %>%
    summarise(
      melt_amt = sum(melt, na.rm = TRUE),
      total_rain = first(total_rain),
      total_snow = first(total_snow),
      total_prcp = first(total_prcp),
      total_melt_events = first(total_melt_events),
      max_melt_daily = first(max_melt_daily),
      swe_ratio = first(swe_ratio),
      ws_swe_ratio= first(ws_swe_ratio),
      .groups = "drop"
    ) %>%
    
    filter(melt_amt > 0) %>%
    
    group_by(monitoring_location_id, water_year) %>%
    summarise(
      total_rain = first(total_rain),
      total_snow = first(total_snow),
      total_prcp = first(total_prcp),
      total_melt_events = first(total_melt_events),
      swe_ratio = first(swe_ratio),
      ws_swe_ratio= first(ws_swe_ratio),
      
      total_melt = sum(melt_amt, na.rm = TRUE),
      avg_melt   = mean(melt_amt, na.rm = TRUE),
      min_melt   = min(melt_amt, na.rm = TRUE),
      max_melt   = max(melt_amt, na.rm = TRUE),
      
      max_melt_daily = first(max_melt_daily),
      n_melt_events = n(),
      .groups = "drop"
    ) %>% 
    select(-total_melt, -total_snow, -total_prcp, -n_melt_events)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}




#this function calculates precipitation timing signatures
calculate_prcp_timing <- function(climate_data, save_path = NULL) {
  #calculate annual totals first: sum up the daily depths for winter-spring of every site year
  annual_totals<- climate_data %>% 
    group_by(monitoring_location_id, water_year) %>% 
    mutate(delta_swe = swe - lag(swe),
           snow = ifelse(prcp > 0 & delta_swe > 0, prcp, 0)
    ) %>% 
    
    summarise(annual_prcp= sum(prcp),
              annual_swe= sum(snow, na.rm= TRUE))
  # then join to daily 
  
  daily_frac <- climate_data %>%
    left_join(annual_totals, by = c("monitoring_location_id", "water_year")) %>%
    arrange(monitoring_location_id, water_year, wy_doy) %>%
    group_by(monitoring_location_id, water_year) %>%
    mutate(
      delta_swe = swe - lag(swe),
      snow = ifelse(prcp > 0 & delta_swe > 0, prcp, 0),
      ros= ifelse(prcp >0 & delta_swe <= 0 & swe >0, prcp, 0), #this is the depth of ROS
      ros_event = coalesce(ros > 0, FALSE),
      ros_event_num = cumsum(ros_event & !lag(ros_event, default = FALSE)),
      frac_annual_swe = snow / annual_swe, 
      frac_annual_prcp= prcp/ annual_prcp,
      cum_frac_swe = cumsum(coalesce(frac_annual_swe, 0)), #had to do this to avoid NA values
      cum_frac_prcp= cumsum(frac_annual_prcp)
    )
  
  #cv and iqd
  prcp_timing <- daily_frac %>%
    summarise(
      doy_25_swe = wy_doy[which(cum_frac_swe >= 0.25)[1]],
      CV_swe = wy_doy[which(cum_frac_swe >= 0.50)[1]],
      doy_75_swe = wy_doy[which(cum_frac_swe >= 0.75)[1]],
      IQD_swe = doy_75_swe - doy_25_swe,
      doy_25_prcp = wy_doy[which(cum_frac_prcp >= 0.25)[1]],
      CV_prcp = wy_doy[which(cum_frac_prcp >= 0.50)[1]],
      doy_75_prcp = wy_doy[which(cum_frac_prcp >= 0.75)[1]],
      IQD_prcp = doy_75_prcp - doy_25_prcp,
      peak_swe_doy= wy_doy[which.max(swe)],
      swe_peak= max(swe),
      snow_cover_duration= sum(swe> 0),
      total_ros= sum(ros),
      num_ros_events= max(ros_event_num),
      avg_ros_depth= mean(ros, na.rm= TRUE),
      max_ros_depth= max(ros, na.rm=TRUE),
      .groups = "drop"
    )
  
  output_df <- annual_totals %>%
    left_join(prcp_timing, by = c("monitoring_location_id", "water_year")) %>%
    select(-doy_25_swe, -doy_75_swe, -doy_25_prcp, -doy_75_prcp) %>% 
    arrange(monitoring_location_id, water_year)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)


}



#this function calculates precipitation seasonality which is the sum of deviations
#of monthly precip from mean montly precip, divided by the total annual precip
#values range from 0 to greater than 1.2
calculate_prcp_seasonality<- function(climate_data, save_path = NULL) {
  
  output_df <- climate_data %>%
    group_by(monitoring_location_id, water_year, month) %>%
    summarise(
      month_mean = mean(prcp, na.rm = TRUE),
      month_total= sum(prcp, na.rm= TRUE),
      .groups = "drop"
    ) %>%
    
    group_by(monitoring_location_id, water_year) %>%
    summarise(
      Ri = sum(month_total, na.rm = TRUE),
      
      seasonality_annual = #could make this more generic, but not hard-coding number of months
        sum(abs(month_mean - Ri/12), na.rm = TRUE) * (1/ Ri), #now this ranges from 0 to 2
      
      .groups = "drop"
    ) %>%
    select(-Ri) %>% 
    
    arrange(monitoring_location_id, water_year)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}





calculate_seasonal_prcp <- function(climate_data, save_path = NULL) {
  
  output_df <- climate_data %>%
    group_by(monitoring_location_id, water_year) %>%
    summarise(
      fall_prcp   = sum(prcp[month %in% 10:12], na.rm = TRUE),
      winter_prcp = sum(prcp[month %in% 1:3],  na.rm = TRUE),
      spring_prcp = sum(prcp[month %in% 4:6],  na.rm = TRUE),
      summer_prcp = sum(prcp[month %in% 7:9],  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(monitoring_location_id, water_year)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}





calculate_seasonal_airtemp <- function(climate_data, save_path = NULL) {
  
  output_df <- climate_data %>%
    group_by(monitoring_location_id, water_year) %>%
    summarise(
      fall_airtemp   = (mean(tmax[month %in% 10:12], na.rm = TRUE) + mean(tmin[month %in% 10:12], na.rm = TRUE)) / 2,
      winter_airtemp = (mean(tmax[month %in% 1:3], na.rm = TRUE) + mean(tmin[month %in% 1:3], na.rm = TRUE)) / 2,
      spring_airtemp = (mean(tmax[month %in% 4:6], na.rm = TRUE) + mean(tmin[month %in% 4:6], na.rm = TRUE)) / 2,
      summer_airtemp = (mean(tmax[month %in% 7:9], na.rm = TRUE) + mean(tmin[month %in% 7:9], na.rm = TRUE)) / 2,
      .groups = "drop"
    ) %>%
    arrange(monitoring_location_id, water_year)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}




##-----------add in a function that computes spi at monthly scale----------#####


