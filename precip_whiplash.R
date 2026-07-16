"Isolated dry = precipitation is lower than the 20th percentile and previous timestep  precipitation was greater than 20th percentile 
& change from previous timestep did not exceed 60 percentile points.
Isolated wet = precipitation exceeds the 80th percentile and previous timestep  precipitation was less than 80th percentile 
& change from previous timestep did not exceed 60 percentile points. 
Recurring dry (dry-to-dry) = precipitation is lower than 20th percentile & previous  timestep precipitation was lower than 20th percentile. 
Recurring wet (wet-to-wet) = precipitation is higher than 80th percentile & previous  timestep precipitation was higher than 80th percentile.
Dry-to-wet whiplash (dry-to-wet) = precipitation percentile increased by at least 60  percentile points from previous timestep.
Wet-to-dry whiplash (wet-to-dry) = precipitation percentile decreased by at least 60  percentile points from previous timestep."


#this function determines extreme type at the annual scale based on the above logic for every site (water)year

calculate_prcp_annual_extremes<- function(climate_data, save_path){
  
  annual_prcp<- climate_data %>%
    group_by(monitoring_location_id, water_year) %>% 
    summarise(
      total_annual_prcp= sum(prcp, na.rm= TRUE),
      .groups= "drop") #annual totals for precip at each site
  
  annual_prcp<- annual_prcp %>%
    group_by(monitoring_location_id) %>% 
    mutate(prcp_percentile= percent_rank(total_annual_prcp) * 100,
           percentile_change= prcp_percentile - lag(prcp_percentile)) %>% #change from previous timestep
    ungroup() #this is now percentiles for annual precip throughout the time pd & change from previous
  
  annual_prcp<- annual_prcp %>% 
    group_by(monitoring_location_id) %>% 
    mutate(extreme_class= case_when(
      is.na(percentile_change) ~ "NA",
      prcp_percentile > 80 & lag(prcp_percentile) < 80 & percentile_change <= 60 ~ "isolated wet",
      prcp_percentile < 20 & lag(prcp_percentile) > 20 & percentile_change >= -60 ~ "isolated dry",
      prcp_percentile < 20 & lag(prcp_percentile) < 20 ~ "recurring dry",
      prcp_percentile > 80 & lag(prcp_percentile) > 80 ~ "recurring wet",
      percentile_change >= 60 ~ "dry-to-wet",
      percentile_change <= -60 ~ "wet-to-dry",
      TRUE ~ NA_character_ )) %>% 
    ungroup()
  
  if (!is.null(save_path)) {
    saveRDS(annual_prcp, save_path)
  }
  
  return(annual_prcp)
}

#calculate_prcp_annual_extremes(daymet_annual_clean, save_path = NULL)



#this function calculates extreme type for seasons

calculate_prcp_szn_extremes<- function(climate_data, save_path){
  
  season_order<- c("fall", "winter", "spring", "summer") #aligned with water_year
  
  szn_prcp<- daymet_annual_clean %>% 
    mutate(season= case_when(
      month %in% 1:3 ~ "winter",
      month %in% 4:6 ~ "spring",
      month %in% 7:9 ~ "summer",
      month %in% 10:12 ~ "fall"
    )) %>% 
    group_by(monitoring_location_id, water_year, season) %>% 
    summarise(total_szn_prcp= sum(prcp, na.rm= TRUE),
              .groups= "drop")
  
  szn_prcp<- szn_prcp %>% 
    mutate(season= factor(season, levels= season_order)) %>% 
    arrange(monitoring_location_id, water_year, season) #ordering df so that percentile change is correct (spring after winter, etc.)
  
  szn_prcp<- szn_prcp %>% 
    group_by(monitoring_location_id, season) %>% 
    mutate(szn_prcp_percentile= percent_rank(total_szn_prcp) * 100 #this is comparison for seasons e.g.) spring vs. spring
    ) %>% 
    ungroup()
  
  szn_prcp<- szn_prcp %>% 
    group_by(monitoring_location_id) %>% 
    mutate(szn_percentile_change= szn_prcp_percentile - lag(szn_prcp_percentile)) %>% #this needs to be here to compare consecutive seasons
    mutate(
      szn_extreme_class= case_when(
        is.na(szn_percentile_change) ~ "NA",
        szn_prcp_percentile > 80 & lag(szn_prcp_percentile) < 80 & szn_percentile_change <= 60 ~ "isolated wet",
        szn_prcp_percentile < 20 & lag(szn_prcp_percentile) > 20 & szn_percentile_change >= -60 ~ "isolated dry",
        szn_prcp_percentile < 20 & lag(szn_prcp_percentile) < 20 ~ "recurring dry",
        szn_prcp_percentile > 80 & lag(szn_prcp_percentile) > 80 ~ "recurring wet",
        szn_prcp_percentile >= 60 ~ "dry-to-wet",
        szn_prcp_percentile <= -60 ~ "wet-to-dry",
        TRUE ~ NA_character_ )) %>% 
    ungroup()
  
  if(!is.null(save_path)){
    saveRDS(szn_prcp, save_path)
  }
  
  return(szn_prcp)
}

calculate_prcp_szn_extremes(daymet_annual_clean, save_path= NULL)


#some sort of fn to take extremes output and summarize them:
# 1. number of each event type at a site over time
# 2. spatial plot of linear slopes?
# need to pivot to wide format e.g.) extreme_class becomes recurring dry, recurring wet, etc.?


#summarize the number of each type of extreme at each site over time period--show distribution
summary_annual_prcp<- annual_prcp %>% 
  group_by(monitoring_location_id) %>% 
  summarise(n_drytowet= n_distinct(extreme_class== "dry to wet"),
            n_wettodry= n_distinct(extreme_class== "wet to dry"),
            n_recurringwet= n_distinct(extreme_class== "recurring wet"),
            n_recurringdry= n_distinct(extreme_class== "recurring dry"),
            n_isolatedwet= n_distinct(extreme_class== "isolated wet"),
            n_isolateddry= n_distinct(extreme_class== "isolated dry"))
                                  
