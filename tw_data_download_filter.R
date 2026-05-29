"This script includes the workflow for pulling daily stream temp data
from NWIS, filtering based on data availability, and interpolation for missing values.
The end result is a dataframe that can be used to run stream temperature signatures"


#necessary libraries:
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

##-------------------------require 300 days of data per water year, remove site years with long gaps-----------####
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
site_year_summary <- dv_tw_full %>%
  group_by(monitoring_location_id, water_year) %>%
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


#Keep only good site-years
# at least 300 observed days
# no NA gaps > 5 days

good_site_years <- site_year_summary %>%
  filter(
    n_obs >= 250,
    max_gap <= 30 #some gages are operated seasonally... toss out months with inadequate data? or require completeness march- nov?
  )

dv_tw_adequate_data <- dv_tw_full %>%
  semi_join(
    good_site_years,
    by = c("monitoring_location_id", "water_year")
  )


##-------------------------------filter to sites with at least 11 yrs of data-----------------####

dv_tw_adequate_data <- dv_tw_adequate_data %>%
  group_by(monitoring_location_id) %>%
  filter(n_distinct(water_year[water_year <= 2022]) >= 11)

dv_tw_adequate_data %>% 
  summarise(n_sites= n_distinct(monitoring_location_id))


##-------------------------------------fill gaps <= 5 days using either linear interpolation or PCA imputation-----####

#this function interpolates gaps with linear interpolation or with PCA depending on
#the amount of correlated data available
hybrid_impute_tw <- function(df,
                             maxgap = 7, #originally set this to 5--can iterate on this number
                             min_sites = 3,
                             min_overlap = 0.3,
                             ncp = 2) {
  
  df <- df %>%
    arrange(monitoring_location_id, time) %>%
    mutate(time = as.Date(time))
  

# Linear interpolation first

  df_lin <- df %>%
    group_by(monitoring_location_id) %>%
    arrange(time) %>%
    mutate(
      tw_linear = na.approx(tw, maxgap = maxgap, na.rm = FALSE),
      method = ifelse(is.na(tw) & !is.na(tw_linear), "linear", NA)
    ) %>%
    ungroup()
  

#process each year to see whether there are enough sites for PCA
  df_out <- df_lin %>%
    group_split(year) %>%
    lapply(function(d) {
      
      # pivot to wide
      wide <- d %>%
        select(time, monitoring_location_id, tw_linear) %>%
        pivot_wider(names_from = monitoring_location_id,
                    values_from = tw_linear) %>%
        arrange(time)
      
      time_index <- wide$time
      mat <- wide %>% select(-time)
      
      # diagnostics
      n_sites <- ncol(mat)
      overlap <- mean(!is.na(mat))
      
      # skip PCA if not enough structure
      if (n_sites < min_sites || overlap < min_overlap) {
        return(d %>% mutate(tw_filled = tw_linear))
      }
      
      # run PCA imputation safely
      imputed_mat <- tryCatch({
        imputePCA(mat, ncp = ncp, maxiter = 100)$completeObs
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(imputed_mat)) {
        return(d %>% mutate(tw_filled = tw_linear))
      }
      
      # back to long
      df_pca <- imputed_mat %>%
        as.data.frame() %>%
        mutate(time = time_index) %>%
        pivot_longer(-time,
                     names_to = "monitoring_location_id",
                     values_to = "tw_pca")
      
      # merge back
      d2 <- d %>%
        left_join(df_pca,
                  by = c("time", "monitoring_location_id")) %>%
        mutate(
          tw_filled = case_when(
            !is.na(tw) ~ tw,
            is.na(tw) & !is.na(tw_linear) ~ tw_linear,
            is.na(tw) & is.na(tw_linear) & !is.na(tw_pca) ~ tw_pca,
            TRUE ~ NA_real_
          ),
          method = case_when(
            !is.na(tw) ~ "observed",
            is.na(tw) & !is.na(tw_linear) ~ "linear",
            is.na(tw) & is.na(tw_linear) & !is.na(tw_pca) ~ "pca",
            TRUE ~ "missing"
          )
        )
      
      return(d2)
    }) %>%
    bind_rows()
  
  return(df_out)
}

df_tw_interp <- hybrid_impute_tw(dv_tw_adequate_data) #calling the fn on the filtered tw dataset

#also need to filter to sites with at least 11 yrs data

df_tw_interp %>% 
  summarise(n_sites= n_distinct(monitoring_location_id))

