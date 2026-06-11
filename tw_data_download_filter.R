"This script includes the workflow for pulling daily stream temp data
from NWIS, filtering based on data availability, and interpolation for missing values.
The end result is a dataframe that can be used to run stream temperature signatures"


#necessary libraries:
library(dataRetrieval)
library(stringr)
library(imputeTS)
library(tidyverse)
library(zoo)
library(purrr)
library(sf)
library(stringr)
library(imputeTS)
library(missMDA)
library(FactoMineR)

##---------------download tw data from NWIS for sites with daymet data---------------------####

daymet_sites<- readRDS("daymet_sites_summary.RDS") #these are GL sites with daymet data from 1980 to 2022

#converting site_id to monitoring_location_id for NWIS request
ids<- unique(daymet_sites$site_id)
ids<- paste0("USGS-", ids) # need to be in "USGS-" format for requesting



# downloading daily mean approved values for temp from sites
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

dv_tw <- map(
  batch_ids(ids, n = 20),
  ~ safe_get_dv(.x, "00010")
) |>
  purrr::compact() |>
  purrr::list_rbind()

##--------------------de-duplicate dataframe------------####
#now that data is downloaded, cleaning up columns
dv_tw_filtered<- dv_tw %>%  #also want to remove geometry for now (will request separately with metadata)
  select(c(-"time_series_id", -"last_modified")) %>%  #Unnecessary columns
  filter(approval_status== "Approved") %>%  #keeping only approved data
  rename(tw= value) #renaming to parameter

dv_tw_filtered%>%
  count(monitoring_location_id) %>%
  arrange(desc(n)) %>%
  head()

# Do any site-date pairs duplicate?
dup_check <- dv_tw_filtered %>% #temporary dv
  count(monitoring_location_id, time) %>%
  filter(n > 1)
dup_check  # should be empty

#if it's not empty, remove duplicates:

dv_tw_filtered <- dv_tw_filtered %>%
  distinct(monitoring_location_id, time, .keep_all = TRUE)

#double check that it's empty now:
dup_check <- dv_tw_filtered %>% #temporary dv
  count(monitoring_location_id, time) %>%
  filter(n > 1)
dup_check  # should be empty



##-----------------------------re-doing completeness criteria w monthly scale-----------------####
#group by month and require 20 or more days per month and consecutive gaps of 5 days or less

#at this point, dataframe should be de-duplicated, but doesn't need to be filtered to adequate data yet

dv_tw_full <- dv_tw_filtered %>%
  mutate(
    time = as.Date(time),
    year = year(time),
    month= month(time),
    water_year = ifelse(month >= 10, year + 1, year),
    wy_doy = as.integer(time - ymd(paste0(water_year - 1, "-10-01"))) + 1
  ) %>%
  
  group_by(monitoring_location_id, water_year) %>% #GROUP BY WY
  
  # create complete daily sequence within each site-year
  
  complete(
    time = seq(
      floor_date(min(time), "month"),
      ceiling_date(max(time), "month") - days(1), #need to do this so that gaps at beginning or end of month are counted!
      by = "day"
    )
  ) %>%
  
  #refill metadata columns first (better than interpolating tw first)
  fill(
    parameter_code,
    statistic_id,
    unit_of_measure,
    approval_status,
    qualifier,
    .direction = "downup"
  ) %>%
  
  
  mutate(
    year = year(time),
    month = month(time),
    water_year = ifelse(month >= 10, year + 1, year),
    wy_doy = as.integer(time - ymd(paste0(water_year - 1, "-10-01"))) + 1
  ) %>% #needed to re-derive this once filled!
  
  
  ungroup()

#calculate coverage and longest NA gap
site_year_summary_month <- dv_tw_full %>%
  filter(water_year<= 2022) %>% 
  group_by(monitoring_location_id, water_year, month) %>%
  summarize(
    
    #number of non-NA discharge values
    n_obs = sum(!is.na(tw)),
    
    #longest consecutive NA run
    max_gap = {
      r <- rle(is.na(tw))
      max(c(0, r$lengths[r$values]))
    },
    
    .groups = "drop"
  )
#need to filter this df and semi-join so that acceptable months at a site are kept..
adequate_tw_data<- site_year_summary_month %>% 
  filter(n_obs>= 20 & max_gap <= 5)

#need to join to create a dataframe that meets completeness criteria for interpolation, etc below!!
dv_tw_adequate_data <- dv_tw_full %>%
  semi_join(
    adequate_tw_data,
    by = c("monitoring_location_id", "water_year", "month") #including month in join since filter criteria used it
  )


#OR: several dataframes? one annual, one summer.. to optimize available data


##-------------------------------filter to sites with at least 11 yrs of data-----------------####

dv_tw_adequate_data <- dv_tw_adequate_data %>%
  group_by(monitoring_location_id) %>%
  filter(n_distinct(water_year[water_year <= 2022]) >= 11)

tw_sites<- dv_tw_adequate_data %>% 
  summarise(n_sites= n_distinct(monitoring_location_id))

excluded_sites<- setdiff(ids, sites)

!tw_sites$monitoring_location_id %in% excluded_sites #one site needs to be eliminated sadly

dv_tw_adequate_data<- dv_tw_adequate_data %>% 
  filter(!monitoring_location_id%in% excluded_sites)

#re-echecking:
tw_sites<- dv_tw_adequate_data %>% 
  summarise(n_sites= n_distinct(monitoring_location_id)) #confirmed that it dropped.

##-------------------------------------fill gaps <= 5 days using linear interpolation-----####

#handling gaps with linear interpolation since gaps are limited to 5 consecutive days or less. 
#also have a version saved that uses PCA imputation, but prob overkill for such short gaps..

dv_tw_interp <- dv_tw_adequate_data %>%
  arrange(monitoring_location_id, time) %>%
  
  group_by(monitoring_location_id, water_year, month) %>%
  
  mutate(
    
    # flag rows that were originally missing
    interpolated = is.na(tw),
    
    # interpolate only gaps <= 5 days
    tw = na_interpolation(
      tw,
      option = "linear",
      maxgap = 5
    )
  ) %>%
  
  ungroup()


cleaned_dv_tw<- dv_tw_interp #renaming before saving!


##----------------------------------------saving df to processed folder---------####

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
  cleaned_dv_tw = cleaned_dv_tw)
  

##--------------------------------------sanity checks------------------------####

monthly_completeness <- cleaned_dv_tw %>%
  mutate(
    year = year(time),
    month = month(time)
  ) %>%
  group_by(monitoring_location_id, year, month) %>%
  summarise(
    n_obs = n(),
    expected_days = days_in_month(min(time)),
    complete = n_obs == expected_days,
    .groups = "drop"
  )

#this confirms that I'm not missing any data per month



