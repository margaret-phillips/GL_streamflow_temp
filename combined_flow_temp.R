"This script analyzes streamflow and stream temperature combined using data downloaded
in the download scripts."

# necessary libraries:




#combine streamflow and temp by mutual days
combined_dv_temp_flow<- inner_join(cleaned_dv_tw, cleaned_dv_qDepth_annual, by= c("time", "monitoring_location_id")) %>% 
  select(monitoring_location_id, time, tw, q, rep_drainage_area, daily_depth_mm, wy_doy.y, doy, month.y, year.y, water_year.y) %>% 
  rename(water_year= water_year.y,
         wy_doy= wy_doy.y,
         month= month.y,
         year= year.y)
  