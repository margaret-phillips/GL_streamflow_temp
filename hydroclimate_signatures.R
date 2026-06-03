"This script includes functions to calcualte hydroclimate signatures using daily air 
temperature, SWE, and precipitation data for each site in the daymet dataset. Signatures
are calculated at the annual timescale, with some applying only to specific time periods
within the year. The output signatures dataframes are then processed by a wrapper function
to compute trends, significance, and variability in signatures over the time period."


#required libraries:
library(tidyverse)


#required columns: swe, tmin, tmax, prcp, monitoring_location_id, water_year
daymet_ws<- daymet_ws %>% 
  group_by(monitoring_location_id, year) %>% 
  mutate(doy= row_number())
#NEED TO ADD A FEW SIGNATURES FOR ANNUAL SCALE! and make sure annual daymet has required cols

##-----------air temperature---------------------------------------------------####

#This function calculates minimum, maximum, and average air temp
calculate_air_temp<- function(climate_data, save_path= NULL){
  output_df<- climate_data %>% 
    group_by(monitoring_location_id) %>% 
    mutate(tavg= (tmin + tmax)/2) %>% 
    group_by(monitoring_location_id, year) %>% 
    summarise(tmin= min(tmin),
              tmax= max(tmax),
              tavg= mean(tavg))
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

calculate_air_temp(daymet_ws, save_path= NULL) #calling the fn



#This function calculates the last spring day below -2.2C and number of days above 5C
calculate_spring_days<- function(climate_data, save_path= NULL){
  Tgrow<- 5
  Tdamage<- -2.2
  
  output_df<- climate_data %>% 
    group_by(monitoring_location_id, year) %>% 
    summarise(
      growing_deg = sum(pmax(tmin - Tgrow, 0), na.rm = TRUE), #cumulative degrees above 5C
      last_damage_day = max(doy[tmin < Tdamage], na.rm = TRUE), #doy that corresponds to last value below -2.2C
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
  
  output_df <- climate_data %>%
    
    group_by(monitoring_location_id, year) %>%
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
      swe_ratio = total_snow / total_prcp
    ) %>%
    
    group_by(monitoring_location_id, year, melt_event_num) %>%
    summarise(
      melt_amt = sum(melt, na.rm = TRUE),
      total_rain = first(total_rain),
      total_snow = first(total_snow),
      total_prcp = first(total_prcp),
      total_melt_events = first(total_melt_events),
      max_melt_daily = first(max_melt_daily),
      swe_ratio = first(swe_ratio),
      .groups = "drop"
    ) %>%
    
    filter(melt_amt > 0) %>%
    
    group_by(monitoring_location_id, year) %>%
    summarise(
      total_rain = first(total_rain),
      total_snow = first(total_snow),
      total_prcp = first(total_prcp),
      total_melt_events = first(total_melt_events),
      swe_ratio = first(swe_ratio),
      
      total_melt = sum(melt_amt, na.rm = TRUE),
      avg_melt   = mean(melt_amt, na.rm = TRUE),
      min_melt   = min(melt_amt, na.rm = TRUE),
      max_melt   = max(melt_amt, na.rm = TRUE),
      
      max_melt_daily = first(max_melt_daily),
      n_melt_events = n(),
      .groups = "drop"
    )
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

rain_snow<- calculate_rain_snow(daymet_ws, save_path= NULL) #calling the fn




#this function calculates precipitation timing signatures
calculate_prcp_timing <- function(climate_data, save_path = NULL) {
  #calculate annual totals first: sum up the daily depths for winter-spring of every site year
  annual_totals<- climate_data %>% 
    group_by(monitoring_location_id, year) %>% 
    mutate(delta_swe = swe - lag(swe),
           snow = ifelse(prcp > 0 & delta_swe > 0, prcp, 0)
    ) %>% 
    
    summarise(annual_prcp= sum(prcp),
              annual_swe= sum(snow, na.rm= TRUE))
  # then join to daily 
  
  daily_frac <- climate_data %>%
    left_join(annual_totals, by = c("monitoring_location_id", "year")) %>%
    arrange(monitoring_location_id, year, doy) %>%
    group_by(monitoring_location_id, year) %>%
    mutate(
      delta_swe = swe - lag(swe),
      snow = ifelse(prcp > 0 & delta_swe > 0, prcp, 0),
      frac_annual_swe = snow / annual_swe, 
      frac_annual_prcp= prcp/ annual_prcp,
      cum_frac_swe = cumsum(coalesce(frac_annual_swe, 0)), #had to do this to avoid NA values
      cum_frac_prcp= cumsum(frac_annual_prcp)
    )
  
  #cv and iqd
  prcp_timing <- daily_frac %>%
    summarise(
      doy_25_swe = doy[which(cum_frac_swe >= 0.25)[1]],
      CV_swe = doy[which(cum_frac_swe >= 0.50)[1]],
      doy_75_swe = doy[which(cum_frac_swe >= 0.75)[1]],
      IQD_swe = doy_75_swe - doy_25_swe,
      doy_25_prcp = doy[which(cum_frac_prcp >= 0.25)[1]],
      CV_prcp = doy[which(cum_frac_prcp >= 0.50)[1]],
      doy_75_prcp = doy[which(cum_frac_prcp >= 0.75)[1]],
      IQD_prcp = doy_75_prcp - doy_25_prcp,
      .groups = "drop"
    )
  
  output_df <- annual_totals %>%
    left_join(prcp_timing, by = c("monitoring_location_id", "year")) %>%
    arrange(monitoring_location_id, year)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(list(
    annual_metrics = output_df,
    daily_metrics = daily_frac #maybe don't return this after all?
  ))
}

calculate_prcp_timing(daymet_ws, save_path = NULL)



#this function calculates precipitation seasonality which is the sum of deviations
#of monthly precip from mean montly precip, divided by the total annual precip
#values range from 0 to greater than 1.2
calculate_prcp_seasonality<- function(climate_data, save_path = NULL) {
  
  output_df <- climate_data %>%
    group_by(monitoring_location_id, water_year, month) %>%
    summarise(
      month_mean = mean(prcp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    group_by(monitoring_location_id, water_year) %>%
    summarise(
      Ri = sum(prcp, na.rm = TRUE),
      
      seasonality_annual = #could make this more generic, but not hard-coding number of months
        sum(abs(month_mean - Ri/12), na.rm = TRUE) / (1/ Ri), #now this ranges from 0 to 2
      
      .groups = "drop"
    ) %>%
    
    arrange(monitoring_location_id, water_year)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}