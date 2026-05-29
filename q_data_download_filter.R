" This script includes the workflow for pulling daily streamflow and stream temp
from NWIS, filtering based on data availability, and interpolation for missing values.
NOTE: downloaded daily values do not include NAs for missing data, so NA values must
be added to the dataset based on missing dates and all columns must be interpolated.
The script also downloads metadata to compute a regulation metric and remove substantially
regulated sites. The end result is two dataframes (one annual and one winter-spring) that
are used to run all functions and compute streamflow and stream temp signatures."

# packages needed:
library(dataRetrieval)
library(tidyverse)
library(zoo)
library(purrr)
library(sf)
library(stringr)
library(imputeTS)

# conversion factors:
ft3_per_ML<- 35314.6667 #conversion factor for ft3 to megaliters
sec_per_day<- 86400 #conversion factor for seconds per day
km2_per_mi2<- 2.58998811 #square km per sq mile
ft2_per_km2<- 1.0764E+7 #square ft per sq kilometer
mm_per_ft<- 304.8 #mm per foot

##---------------requesting data for USGS sites with daymet data-------------####
daymet_sites<- readRDS("daymet_sites_summary.RDS") #these are GL sites with daymet data from 1980 to 2022

#converting site_id to monitoring_location_id for NWIS request
ids<- unique(daymet_sites$site_id)
ids<- paste0("USGS-", ids) # need to be in "USGS-" format for requesting



# downloading daily mean approved values for streamflow and temp from sites
##batching site ids to avoid exceeding request limits
start_date<- "1979-10-01" # CHANGE THIS TO WY: 1979-10-01
end_date<- "2025-09-30"

batch_ids <- function(x, n = 20) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  split(x, ceiling(seq_along(x) / n))
}


drop_geom <- function(x) { #this didn't work for some reason
  if (inherits(x, "sf")) {
    sf::st_drop_geometry(x)
  } else {
    x
  }
}

## this function gets daily values for the start and end date, parameter code, and statistic specified
get_dv <- function(id_batch, pcode) {
  read_waterdata_daily(
    monitoring_location_id = id_batch,
    parameter_code         = pcode,
    statistic_id           = "00003",  # daily mean
    time                   = c(start_date, end_date)
  ) |>
    drop_geom() #this makes it less demanding and you can request metadata afterwards, but it didn't work, so dropped it later
}


safe_get_dv <- purrr::possibly(
  get_dv,
  otherwise = NULL
)


dv_q <- map(
  batch_ids(ids, n = 20),
  ~ safe_get_dv(.x, "00060")
) |>
  purrr::compact() |>
  purrr::list_rbind()

#how does the download process treat NA values
sum(!complete.cases(dv_tw)) #confirmed.. rows with NA values are not downloaded. need to manually add them

#now that data is downloaded, cleaning up columns
dv_q_filtered<- dv_q %>%  #also want to remove geometry for now (will request separately with metadata)
  select(c(-"geometry", -"time_series_id", -"last_modified")) %>%  #Unnecessary columns
  filter(approval_status== "Approved") %>%  #keeping only approved data
  rename(q= value) #renaming to parameter

##--------------------de-duplicate dataframe------------####
dv_q_filtered%>%
  count(monitoring_location_id) %>%
  arrange(desc(n)) %>%
  head()

# Do any site-date pairs duplicate?
dup_check <- dv_q_filtered %>% #temporary dv
  count(monitoring_location_id, time) %>%
  filter(n > 1)
dup_check  # should be empty

#if it's not empty, remove duplicates:

dv_q_filtered <- dv_q_filtered %>%
  distinct(monitoring_location_id, time, .keep_all = TRUE)

#double check that it's empty now:
dup_check <- dv_q_filtered %>% #temporary dv
  count(monitoring_location_id, time) %>%
  filter(n > 1)
dup_check  # should be empty


##-------------------------require 300 days of data per water year, remove site years with long gaps, fill gaps <= 5 days----------####
#at this point, dataframe should be de-duplicated, but doesn't need to be filtered to adequate data yet

dv_q_full <- dv_q_filtered %>%
  mutate(
    time = as.Date(time),
    year = year(time),
    month= month(time),#ADD WATER YEAR HERE
    water_year = ifelse(month >= 10, year + 1, year),
    wy_doy = as.integer(time - ymd(paste0(water_year - 1, "-10-01"))) + 1
  ) %>%
  
  group_by(monitoring_location_id, water_year) %>% #GROUP BY WY
  
  # create complete daily sequence within each site-year
  complete(
    time = seq(min(time), max(time), by = "day")
  ) %>%
  
  #refill metadata columns first (better than interpolating Q first)
  fill(
    parameter_code,
    statistic_id,
    unit_of_measure,
    approval_status,
    qualifier,
    .direction = "downup"
  ) %>%
  
  ungroup()

#calculate coverage and longest NA gap
site_year_summary <- dv_q_full %>%
  group_by(monitoring_location_id, water_year) %>%
  summarize(
    
    #number of non-NA discharge values
    n_obs = sum(!is.na(q)),
    
    #longest consecutive NA run
    max_gap = {
      r <- rle(is.na(q))
      max(c(0, r$lengths[r$values]))
    },
    
    .groups = "drop"
  )


#Keep only good site-years
# at least 300 observed days
# no NA gaps > 5 days

good_site_years <- site_year_summary %>%
  filter(
    n_obs >= 300,
    max_gap <= 5
  )

dv_q_adequate_data <- dv_q_full %>%
  semi_join(
    good_site_years,
    by = c("monitoring_location_id", "water_year")
  )


## interpolate only gaps that are < 5 days long


dv_q_interp <- dv_q_adequate_data %>%
  arrange(monitoring_location_id, time) %>%
  
  group_by(monitoring_location_id, water_year) %>%
  
  mutate(
    
    # flag rows that were originally missing
    interpolated = is.na(q),
    
    # interpolate only gaps <= 5 days
    q = na_interpolation(
      q,
      option = "linear",
      maxgap = 5
    )
  ) %>%
  
  ungroup()


## check interpolated rows

interpolated_rows <- dv_q_interp %>%
  filter(interpolated & !is.na(q))

remaining_gaps <- dv_q_interp %>% #verify that no gaps (data that doesn't meet criteria) remain
  filter(is.na(q))


##-----------------------------add necessary columns to run functions (DOY, month, year)----------####
dv_q_annual <- dv_q_interp %>%
  mutate(
    month = month(time),
    doy = yday(time)
  ) 

##----------------------------------create winter-spring dataframe (Jan- May)----------####
winter_spring<- 1:5

dv_q_ws<- dv_q_annual %>% 
  filter(month%in%winter_spring)

##--------------------------------------now download metadata for sites -------------------####
#at this point, we need metadata in order to:
# 1. remove sites that are substantially regulated (using approach from Dudley et al., 2018 & 2020)
# 2. normalize discharge for watershed area by converting to a flow depth

q_sites<- unique(daymet_sites$monitoring_location_id) #will use this list to request metadata

meta <- read_waterdata_monitoring_location(
  monitoring_location_id = q_sites,
  properties = c("monitoring_location_id",
                 "monitoring_location_name",
                 "state_name",
                 "hydrologic_unit_code",
                 "drainage_area",
                 "contributing_drainage_area",
                 "site_type_code",
                 "geometry")
)

# When a gage sits partway along a flowline rather than at its outlet, the downstream portion of the catchment 
# does not contribute to the gage (https://doi-usgs.github.io/nhdplusTools/articles/drainage_area_estimation.html)
# this is the difference btwn drainage area and contrib drainage area
# not all sites have contrib values and only a few have discrepancies.. will re-visit after regulation filtering

##------------------------------------------compute regulation metric and remove sites that exceed--------------####

# for now, just pulling in dams metadata from gagesii
gagesii_meta<- readxl::read_excel("gagesii_data/gagesII_sept30_2011_conterm.xlsx", sheet= "HydroMod_Dams")

#filtering gagesii meta to GL q sites
gagesii_GL<- gagesii_meta %>% 
  mutate(monitoring_location_id= paste0("USGS-", STAID)) %>% 
  filter(monitoring_location_id%in%q_sites)


regulation_metric <- meta %>%
  select(monitoring_location_id, drainage_area, contributing_drainage_area) %>% #drainage area is in mi2, converting to km2
  mutate(
    rep_drainage_area = pmin(drainage_area, contributing_drainage_area, na.rm = TRUE)* km2_per_mi2 #choosing whichever drainage area is smaller
  )

#computing mean basin annual volume/time

mean_basin_q<- dv_q_annual %>% 
  select(monitoring_location_id, q, water_year) %>% 
  group_by(monitoring_location_id, water_year) %>% 
  mutate(mean_annual= mean(q) * sec_per_day/ ft3_per_ML) %>%  #mean q in ML/ day for every water year
  group_by(monitoring_location_id) %>% 
  summarise(total_mean_q= mean(mean_annual)) #now averaging over the entire period

#joining dam storage column from gagesii and total_mean_q from above
regulation_metric<- regulation_metric %>% 
  left_join(mean_basin_q %>% select(monitoring_location_id, total_mean_q), by= "monitoring_location_id") %>% 
  left_join(gagesii_GL %>% select(monitoring_location_id, STOR_NID_2009), by= "monitoring_location_id") #units: ML/ km2 (megaliter)


# now adding normalized dam storage and flagging those with greater than 180 days
regulation_metric<- regulation_metric %>% 
  mutate(norm_dam_storage= (rep_drainage_area * STOR_NID_2009)/ total_mean_q) %>% #normalized dam storage
  mutate(too_regulated= ifelse(norm_dam_storage>= 180, TRUE, FALSE))

#filtering out sites that exceed normalized dam storage
regulation_filtered<- regulation_metric %>% 
  filter(too_regulated== "FALSE") #this eliminates 8 sites

acceptable_sites<- regulation_filtered$monitoring_location_id

##----------------------------------------------saving filtered and cleaned dataframes for annual and winter-spring--------####

#annual q first--just filtering to sites that meet dam storage criteria (not too regulated):
cleaned_dv_q_annual<- dv_q_annual %>% 
  filter(monitoring_location_id%in%acceptable_sites)

#now meta data from NWIS
cleaned_meta<- meta %>% 
  filter(monitoring_location_id%in%acceptable_sites)

#metadata from gagesii
cleaned_gagesii_GL<- gagesii_GL %>% 
  filter(monitoring_location_id%in%acceptable_sites)


#flow depth per day (annual)
cleaned_dv_qDepth_annual<- cleaned_dv_q_annual %>% 
  left_join(regulation_filtered %>% select(monitoring_location_id, rep_drainage_area), by= "monitoring_location_id") %>%
  mutate(daily_depth_mm= q * sec_per_day/ (rep_drainage_area * ft2_per_km2) * mm_per_ft)

#ws q--just filtering to sites that meet dam storage criteria (not too regulated):
cleaned_dv_q_ws<- dv_q_ws %>% 
  filter(monitoring_location_id%in%acceptable_sites)

#ws flow depth per day
cleaned_dv_qDepth_ws<- cleaned_dv_qDepth_annual %>% 
  filter(month%in%winter_spring)

  
#now save the interpolated complete datasets to be used for running functions (and save versions with flow DEPTH)

#function to save dataframes to output directory as RDS objects for easy re-loading
save_rds <- function(out_dir, ...) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  purrr::iwalk(
    list(...),
    ~ saveRDS(.x, file.path(out_dir, paste0(.y, ".rds")))
  )
}

#calling fn for cleaned and filtered dfs
save_rds(
  out_dir = "data/processed",
  cleaned_dv_q_annual = cleaned_dv_q_annual,
  cleaned_meta = cleaned_meta,
  cleaned_gagesii_GL= cleaned_gagesii_GL,
  cleaned_dv_qDepth_annual = cleaned_dv_qDepth_annual,
  cleaned_dv_q_ws = cleaned_dv_q_ws,
  cleaned_dv_qDepth_ws = cleaned_dv_qDepth_ws
)