# this script includes the workflow for pulling daily streamflow and stream temp
# from NWIS, filtering based on data availability, and interpolation for missing values
# NOTE: downloaded daily values do not include NAs for missing data, so NA values must
# be added to the dataset based on missing dates and all columns must be interpolated

# packages needed:
library(dataRetrieval)
library(tidyverse)
library(zoo)
library(purrr)
library(sf)
library(stringr)
library(imputeTS)


##---------------------------requesting data for sites with daymet data-------------####
daymet_sites<- readRDS("daymet_sites_summary.RDS") #these are GL sites with daymet data from 1980 to 2022

#converting site_id to monitoring_location_id for NWIS request
ids<- unique(daymet_sites$site_id)
ids<- paste0("USGS-", ids) # need to be in "USGS-" format for requesting



# downloading daily mean approved values for streamflow and temp from sites
##batching site ids to avoid exceeding request limits
start_date<- "1980-01-01"
end_date<- "2025-05-31"

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

##-----------------------------------------------de-duplicate dataframe------------####
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


##----------------------------require 300 days of data per year, remove site years with long gaps, fill gaps <= 5 days----------####
#at this point, dataframe should be de-duplicated, but doesn't need to be filtered to adequate data yet

dv_q_full <- dv_q_filtered %>%
  mutate(
    time = as.Date(time),
    year = year(time)
  ) %>%
  
  group_by(monitoring_location_id, year) %>%
  
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
  group_by(monitoring_location_id, year) %>%
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
    by = c("monitoring_location_id", "year")
  )


## interpolate only gaps that are < 5 days long


dv_q_interp <- dv_q_adequate_data %>%
  arrange(monitoring_location_id, time) %>%
  
  group_by(monitoring_location_id, year) %>%
  
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

#################################################################################

#add necessary columns to run functions (DOY, month, year)

#now save the interpolated complete dataset to be used for running functions