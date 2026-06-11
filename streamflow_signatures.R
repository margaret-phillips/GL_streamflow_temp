"This script includes functions to compute streamflow signatures that characterize
magnitude, timing, frequency, duration, and rate of change. Each function is run
with the annual, complete dataframe of daily streamflow normalized to a flow depth.
When a value of streamflow magnitude is reported, the calculation uses flow depth in mm
which is discharge normalized by drainage area. Otherwise, streamflow in cfs is used.
The output data frame is then processed by a wrapper function that computes trends, 
significance and variability in signatures over the time period."



# required libraries:
library(tidyverse)
library(lubridate)
library(slider)

#loading in saved discharge df:
#required columns: q, daily_depth_mm, monitoring_location_id, water_year, doy (day of year), wy_doy
cleaned_dv_qDepth_annual<- readRDS("data/processed/cleaned_dv_qDepth_annual.rds")


##--------streamflow magnitude-----------------------------------------------####

#this function computes flow quartiles for every site year using daily flow depth in mm
calculate_percentiles<- function(Q_data, save_path){
  output_df<- Q_data %>% 
    group_by(monitoring_location_id, water_year) %>% 
    
    summarise(
      
      q10 = quantile(daily_depth_mm, probs = 0.10, na.rm = TRUE),
      q25 = quantile(daily_depth_mm, probs = 0.25, na.rm = TRUE),
      q50 = quantile(daily_depth_mm, probs = 0.50, na.rm = TRUE),
      q75 = quantile(daily_depth_mm, probs = 0.75, na.rm = TRUE),
      q90 = quantile(daily_depth_mm, probs = 0.90, na.rm = TRUE),
      q95 = quantile(daily_depth_mm, probs = 0.95, na.rm = TRUE),
      q99 = quantile(daily_depth_mm, probs = 0.99, na.rm = TRUE),
      .groups = "drop"
    ) 
  
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

calculate_percentiles(cleaned_dv_qDepth_annual, save_path= NULL) #calling the fn



# this function calculates the slope of the flow duration curve and low, mid, and high range
# slopes for every site year
calculate_FDC <- function(Q_data, save_path,
                        min_n = 10, #should increase this
                        high_range = c(0.00, 0.20),
                        mid_range  = c(0.20, 0.80),
                        low_range  = c(0.90, 1.00)) {
  
  slope_calc <- function(exceedance, flow) { #this is unnecessary since I already filtered to adequate data
    if (length(flow) < min_n) return(NA_real_)
    
    log_flow <- log10(flow + 1e-10)
    model <- try(lm(log_flow ~ exceedance), silent = TRUE)
    
    if (inherits(model, "try-error")) NA_real_
    else coef(model)[["exceedance"]]
  }
  
  output_df<- Q_data %>%
    group_by(monitoring_location_id, water_year) %>%
    summarise(
      n_obs = sum(!is.na(q)),
      
      slopes = {
        flows <- q[!is.na(q)]
        n <- length(flows)
        
        if (n < min_n) {
          list(c(
            overall = NA,
            high = NA,
            mid = NA,
            low = NA
          ))
        } else {
          
          sorted_flows <- sort(flows, decreasing = TRUE)
          exceedance <- (1:n) / (n + 1)
          
          overall <- slope_calc(exceedance, sorted_flows)
          
          high_idx <- exceedance >= high_range[1] &
            exceedance <  high_range[2]
          
          mid_idx  <- exceedance >= mid_range[1] &
            exceedance <  mid_range[2]
          
          low_idx  <- exceedance >= low_range[1] &
            exceedance <= low_range[2]
          
          list(c(
            overall = abs(slope_calc(exceedance, sorted_flows)),
            high    = abs(slope_calc(exceedance[high_idx],
                                 sorted_flows[high_idx])),
            mid     = abs(slope_calc(exceedance[mid_idx],
                                 sorted_flows[mid_idx])),
            low     = abs(slope_calc(exceedance[low_idx],
                                 sorted_flows[low_idx])
          )))
        }
      },
      .groups = "drop"
    ) %>%
    select(-n_obs) %>% 
    tidyr::unnest_wider(slopes, names_sep = "_slope")
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

calculate_FDC(cleaned_dv_qDepth_annual, save_path= NULL) #calling the fn

##------------streamflow timing & duration-----------------------------------####


#CV, IQD, seasonality index

#this function computes flow depth, center of volume, and inter-quartile distance which are measures
#of flow timing and protractedness using daily flow depth in mm.
# add a version that does winter-spring timing??
calculate_CV_IQD <- function(Q_data, save_path = NULL) {
  #calculate annual totals first: sum up the daily depths for winter-spring of every site year
  annual_totals<- Q_data %>% 
    group_by(monitoring_location_id, water_year) %>% 
    summarise(annual_depth_mm= sum(daily_depth_mm))
  # then join to daily 

  daily_frac <- Q_data %>%
    left_join(annual_totals, by = c("monitoring_location_id", "water_year")) %>%
    arrange(monitoring_location_id, water_year, wy_doy) %>%
    group_by(monitoring_location_id, water_year) %>%
    mutate(
      frac_annual_qDepth = daily_depth_mm / annual_depth_mm,
      cum_frac_flow = cumsum(frac_annual_qDepth)
    )
  
  #cv and iqd
  flow_timing <- daily_frac %>%
    summarise(
      doy_25 = wy_doy[which(cum_frac_flow >= 0.25)[1]],
      CV = wy_doy[which(cum_frac_flow >= 0.50)[1]],
      doy_75 = wy_doy[which(cum_frac_flow >= 0.75)[1]],
      IQD = doy_75 - doy_25,
      .groups = "drop"
    )
  
  output_df <- annual_totals %>%
    left_join(flow_timing, by = c("monitoring_location_id", "water_year")) %>%
    arrange(monitoring_location_id, water_year) %>% 
    select(c= -doy_25, -doy_75)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)

}

calculate_CV_IQD(cleaned_dv_qDepth_annual, save_path = NULL) #calling function




#this function computes a seasonality index, which ranges from 0 to 1.2 with 
#higher values indicating higher seasonality (Walsh and Lawler, 1981)
calculate_q_seasonality <- function(Q_data, save_path = NULL) {
  
  output_df <- Q_data %>%
    group_by(monitoring_location_id, water_year, month) %>%
    summarise(
      month_mean = mean(daily_depth_mm, na.rm = TRUE),
      month_total = sum(daily_depth_mm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    group_by(monitoring_location_id, water_year) %>%
    summarise(
      Qa = sum(month_total, na.rm = TRUE),
      seasonality_annual =
        sum(abs(month_mean - (Qa / 12)), na.rm = TRUE) * (1 / Qa),
      .groups = "drop"
    ) %>%
    select(-Qa) %>% 
    
    arrange(monitoring_location_id, water_year)
  
  if (!is.null(save_path)) {
    saveRDS(output_df, save_path)
  }
  
  return(output_df)
}

calculate_q_seasonality(cleaned_dv_qDepth_annual, save_path= NULL)
##-----------------streamflow frequency--------------------------------------####

#flow pulses of different % and sensitivity analysis
calculate_flow_pulses <- function(Q_data, save_path = NULL) {
  
  output_df <- Q_data %>%
    mutate(time = as.Date(time)) %>%
    group_by(monitoring_location_id, water_year) %>%
    arrange(time, .by_group = TRUE) %>%
    
    mutate(
      pct_change = (q - lag(q)) / lag(q), #percent change btwn consecutive timesteps
      
      direction = case_when(
        pct_change >= 0.02  ~  1,   # increase ≥ 2%
        pct_change <= -0.02 ~ -1,   # decrease ≥ 2%
        TRUE ~ 0                  # ignore changes less than 2%
      ),
      
      prev_direction = lag(direction) #Do i need this?
    ) %>%
    
    # counting both increase to decrease and decrease to increase
    summarise(
      two_pct_pulses = sum(
        (prev_direction == 1 & direction == -1) |
          (prev_direction == -1 & direction == 1),
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  if (!is.null(save_path)) saveRDS(output_df, save_path)
  
  return(output_df)
}  

calculate_flow_pulses(Q_data= cleaned_dv_qDepth_annual, save_path=NULL)


##-----------------------streamflow rate of change-----------------------------####

#this function computes a rate of change metric, the Richard Baker flashiness index,
#for every site year using daily discharge
calculate_flashiness <- function(Q_data, save_path = NULL) {
  
  output_df <- Q_data %>%
    mutate(time = as.Date(time)) %>%
    group_by(monitoring_location_id, water_year) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(q_diff = abs(q - dplyr::lag(q))) %>%
    summarise(
      num = sum(q_diff, na.rm = TRUE),
      den = sum(q,  na.rm = TRUE),
      rb_index = if_else(den > 0, num / den, NA_real_),
      .groups = "drop"
    ) %>% 
    select(c= -num, -den) #only want to return rb index, wy, and site id
  
  
  if (!is.null(save_path)) saveRDS(output_df, save_path)
  
  return(output_df)
}

calculate_flashiness(cleaned_dv_qDepth_annual, save_path= NULL) #calling the fn

  
  