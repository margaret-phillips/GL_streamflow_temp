#this script does bootstrapping for trends and trend statistics

library(boot)

#MIGHT want to resample residuals instead and then add them back to the trend estimate to get uncertainty

#first need to write a general function that calculates mann-kendall and sen's slope statistics to be called by bootstrap
compute_trend_metric <- function(
    df,
    metric,
    time_col = "water_year",
    min_years = 10
) {
  
  df %>%
    group_by(monitoring_location_id) %>%
    group_modify(~ {
      
      df2 <- .x %>%
        filter(!is.na(.data[[metric]])) %>%
        group_by(.data[[time_col]]) %>%
        summarise(
          value = mean(.data[[metric]], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(.data[[time_col]])
      
      n_years <- n_distinct(df2[[time_col]])
      mean_val <- mean(df2$value, na.rm = TRUE)
      
      if (n_years < min_years) {
        return(
          tibble(
            metric = metric,
            sen_slope = NA_real_,
            sen_pval = NA_real_,
            mk_tau = NA_real_,
            mk_pval = NA_real_,
            mean_value = mean_val,
            n_years = n_years
          )
        )
      }
      
      sen <- trend::sens.slope(df2$value)
      mk  <- trend::mk.test(df2$value)
      
      tibble(
        metric = metric,
        sen_slope = unname(sen$estimates),
        sen_pval = sen$p.value,
        mk_tau = cor(
          df2[[time_col]],
          df2$value,
          method = "kendall"
        ),
        mk_pval = mk$p.value,
        mean_value = mean_val,
        n_years = n_years
      )
    }) %>%
    ungroup()
}

bootstrap_trends <- function(
    df,
    metric,
    n_boot = 100,
    time_col = "water_year"
) {
  
  purrr::map_dfr(seq_len(n_boot), function(i) {
    
    boot_df <- df %>%
      group_by(
        monitoring_location_id,
        .data[[time_col]]
      ) %>%
      slice_sample(
        n = n(),
        replace = TRUE
      ) %>%
      ungroup()
    
    compute_trend_metric(
      boot_df,
      metric = metric,
      time_col = time_col
    ) %>%
      mutate(iteration = i)
  })
}

#calling bootstrapping code for one signature
boot_results <- bootstrap_trends(
  df = streamflow_combined_df,
  metric = "q75",
  n_boot = 100
)






##version that parallelizes the bootstrapping!

library(dplyr)
library(purrr)
library(furrr)
library(trend)

#filtering streamflow sigs to chosen natural flow regime components
signatures_df<- streamflow_cleaned |> #this version excludes years that violate the temporal coverage rule
  dplyr::select(monitoring_location_id, water_year, annual_rb, annual_low_frequency, annual_high_frequency,
         q75, q10, WSIQD, ws_depth_mm, two_pct_pulses)

#filtering hydroclimate signatures:
hydro_sigs_df<- hydroclimate_combined_df |> 
  dplyr::select(monitoring_location_id, water_year, growing_deg, summer_airtemp, fall_airtemp,
                total_ros, annual_prcp, ws_swe_ratio, winter_prcp)

bootstrap_trend_metric_parallel <- function(
    df,
    metric,
    time_col = "water_year",
    n_boot = 500,
    min_years = 10,
    workers = parallelly::availableCores() - 1
) {

  
  site_list <- split(
    df,
    df$monitoring_location_id
  )
  
  future_map_dfr(
    site_list,
    function(site_dat) {
      
      site_id <- unique(site_dat$monitoring_location_id)
      
      dat <- site_dat %>%
        dplyr::select(
          monitoring_location_id,
          dplyr::all_of(time_col),
          dplyr::all_of(metric)
        ) %>%
        filter(!is.na(.data[[metric]])) %>%
        arrange(.data[[time_col]])
      
      n_years <- nrow(dat)
      
      if (n_years < min_years) {
        
        return(
          tibble(
            monitoring_location_id = site_id,
            metric = metric,
            iteration = seq_len(n_boot),
            sen_slope = NA_real_,
            mk_tau = NA_real_,
            mk_pval = NA_real_,
            mean_value = mean(dat[[metric]], na.rm = TRUE),
            n_years = n_years
          )
        )
      }
      
      years <- dat[[time_col]]
      values <- dat[[metric]]
      
      # Robust trend estimate
      sen0 <- trend::sens.slope(values)
      
      slope0 <- as.numeric(sen0$estimates)
      
      # Robust intercept
      intercept0 <- median(
        values - slope0 * years,
        na.rm = TRUE
      )
      
      # Robust fitted trend
      fitted0 <- intercept0 + slope0 * years
      
      # Residuals from robust trend
      residuals0 <- values - fitted0
      
      map_dfr(
        seq_len(n_boot),
        function(iter) {
          
          # Resample residuals
          boot_residuals <- sample(
            residuals0,
            size = length(residuals0),
            replace = TRUE
          )
          
          # Generate bootstrap series
          boot_values <- fitted0 + boot_residuals
          
          sen <- trend::sens.slope(boot_values)
          mk <- trend::mk.test(boot_values)
          
          tibble(
            monitoring_location_id = site_id,
            metric = metric,
            iteration = iter,
            sen_slope = as.numeric(sen$estimates),
            mk_tau = cor(
              years,
              boot_values,
              method = "kendall"
            ),
            mk_pval = mk$p.value,
            mean_value = mean(boot_values),
            n_years = n_years
          )
        }
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
}


#calling bootstrapping fn for one metric
boot_results <- bootstrap_trend_metric_parallel(
  signatures_df,
  metric = "q75",
  n_boot = 500
)

#calling bootstrapping code for all signatures in the df
metrics <- setdiff(
  names(hydro_sigs_df),
  c("monitoring_location_id", "water_year")
)



#removing worker plan from function and just establishing it globally once:
future::plan(
  future::multisession,
  workers = 8
)

hydro_all_boot <- purrr::map_dfr(
  metrics,
  ~ bootstrap_trend_metric_parallel(
    hydro_sigs_df,
    metric = .x,
    n_boot = 500
  )
)

trend_summary <- boot_results %>%
  group_by(
    monitoring_location_id,
    metric
  ) %>%
  summarise(
    sen_median = median(sen_slope, na.rm = TRUE),
    sen_lcl = quantile(
      sen_slope,
      0.025,
      na.rm = TRUE
    ),
    sen_ucl = quantile(
      sen_slope,
      0.975,
      na.rm = TRUE
    ),
    
    tau_median = median(
      mk_tau,
      na.rm = TRUE
    ),
    tau_lcl = quantile(
      mk_tau,
      0.025,
      na.rm = TRUE
    ),
    tau_ucl = quantile(
      mk_tau,
      0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

hydro_boot_trend_summary <- hydro_all_boot %>%
  group_by(
    monitoring_location_id,
    metric
  ) %>%
  summarise(
    sen_median = median(sen_slope, na.rm = TRUE),
    sen_lcl = quantile(
      sen_slope,
      0.025,
      na.rm = TRUE
    ),
    sen_ucl = quantile(
      sen_slope,
      0.975,
      na.rm = TRUE
    ),
    
    tau_median = median(
      mk_tau,
      na.rm = TRUE
    ),
    tau_lcl = quantile(
      mk_tau,
      0.025,
      na.rm = TRUE
    ),
    tau_ucl = quantile(
      mk_tau,
      0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )





ggplot(
  boot_trend_summary,
  aes(
    y = reorder(monitoring_location_id, sen_median),
    x = sen_median
  )
) +
  geom_pointrange(
    aes(
      xmin = sen_lcl,
      xmax = sen_ucl
    )
  ) +
  geom_vline(
    xintercept = 0,
    color = "red",
    linetype = 2
  ) +
  theme_bw()

boot_trend_summary |> 
  summarize(
    sig_inc= sum(tau_median> 0 & tau_lcl<0 & tau_ucl>0)
  )

library(forecast)
gghistogram(x = boot_trend_summary$tau_median, add.kde= TRUE)


ggplot(boot_trend_summary, aes(x = tau_median)) +
  geom_density(
    alpha = 0.3
  ) +
  labs(
    x = "median mk tau of q75 trend",
    y = "Density"
    
  ) +
  geom_vline(xintercept = 0, size = .3, color= "blue", linetype= "dashed")+
  theme_minimal()


ggplot(boot_trend_summary, aes(x= tau_median))+
  geom_histogram(binwidth = 1, boundary= -0.5)+
  labs(
    x= "Number of extreme events at a site",
    y= "number of sites"
  )+
  theme_minimal()

##--------plotted summaries for both clusters---------------------------------####

boot_tau_summary <- all_boot %>%

  group_by(
    monitoring_location_id,
    metric
  ) %>%
  summarise(
    tau_median = median(mk_tau, na.rm = TRUE),
    tau_lcl = quantile(mk_tau, 0.025, na.rm = TRUE),
    tau_ucl = quantile(mk_tau, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

boot_tau_summary <- boot_tau_summary %>%
  mutate(
    metric = recode(
      metric,
      annual_high_frequency = "high flow frequency",
      annual_low_frequency = "low flow frequency",
      annual_rb = "flashiness index",
      q10 = "q10",
      q75= "q75",
      two_pct_pulses= "two percent pulses",
      ws_depth_mm= "winter-spring total flow depth",
      WSIQD= "winter-spring flow protractedness"
    )
  )

boot_tau_summary <- boot_tau_summary %>%
  left_join(kmeans_climate_plotting |> dplyr::select(monitoring_location_id, cluster), by= "monitoring_location_id")

ci_band <- boot_tau_summary %>%
  group_by(metric, cluster) %>%
  summarise(
    ci_lcl = median(tau_lcl, na.rm = TRUE),
    ci_ucl = median(tau_ucl, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(boot_tau_summary, aes(x = tau_median)) +
  
  geom_rect(
    data = ci_band,
    inherit.aes = FALSE,
    aes(
      xmin = ci_lcl,
      xmax = ci_ucl,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "skyblue",
    alpha = 0.3
  ) +
  
  geom_density(

    alpha = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "red"
  ) +
  
  facet_wrap(
    vars(metric, cluster), scales = "free_y") +
  
  theme_minimal() +
  
  labs(
    x = "Median Mann-Kendall Tau",
    y = "Density"
  )

#same thing, but for sen's slope!---------------------------------------------####
boot_sen_summary <- all_boot %>%
  
  group_by(
    monitoring_location_id,
    metric
  ) %>%
  summarise(
    sen_median = median(sen_slope, na.rm = TRUE),
    sen_lcl = quantile(sen_slope, 0.025, na.rm = TRUE),
    sen_ucl = quantile(sen_slope, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

boot_sen_summary <- boot_sen_summary %>%
  mutate(
    metric = recode(
      metric,
      annual_high_frequency = "high flow frequency",
      annual_low_frequency = "low flow frequency",
      annual_rb = "flashiness index",
      q10 = "q10",
      q75= "q75",
      two_pct_pulses= "two percent pulses",
      ws_depth_mm= "winter-spring total flow depth",
      WSIQD= "winter-spring flow protractedness"
    )
  )

boot_sen_summary <- boot_sen_summary %>%
  left_join(kmeans_climate_plotting |> dplyr::select(monitoring_location_id, cluster), by= "monitoring_location_id")

ci_band <-  boot_sen_summary |> 
  group_by(metric, cluster) %>%
  summarise(
    ci_lcl = median(sen_lcl, na.rm = TRUE),
    ci_ucl = median(sen_ucl, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(boot_sen_summary, aes(x = sen_median)) +
  
  geom_rect(
    data = ci_band,
    inherit.aes = FALSE,
    aes(
      xmin = ci_lcl,
      xmax = ci_ucl,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "skyblue",
    alpha = 0.3
  ) +
  
  geom_density(
    
    alpha = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "red"
  ) +
  
  facet_wrap(
    vars(metric, cluster), scales = "free_y") +
  
  theme_minimal() +
  
  labs(
    x = "Median Sen's slope",
    y = "Density"
  )


##--------version without clusters, but with conf intervals for each site-------####

boot_tau_summary <- all_boot %>%
  
  group_by(
    monitoring_location_id,
    metric
  ) %>%
  summarise(
    tau_median = median(mk_tau, na.rm = TRUE),
    tau_lcl = quantile(mk_tau, 0.05, na.rm = TRUE), #0.025 and 0.975 for 95% conf interval
    tau_ucl = quantile(mk_tau, 0.95, na.rm = TRUE),
    .groups = "drop"
  )



boot_tau_summary <- boot_tau_summary %>%
  mutate(
    metric = recode(
      metric,
      annual_high_frequency = "high flow frequency",
      annual_low_frequency = "low flow frequency",
      annual_rb = "flashiness index",
      q10 = "q10",
      q75= "q75",
      two_pct_pulses= "two percent pulses",
      ws_depth_mm= "winter-spring total flow depth",
      WSIQD= "winter-spring flow protractedness"
    )
  )

ci_band <- boot_tau_summary  %>%
  group_by(metric) %>%
  summarise(
    ci_lcl = median(tau_lcl, na.rm = TRUE),
    ci_ucl = median(tau_ucl, na.rm = TRUE),
    .groups = "drop"
  )


flow_boot<- ggplot(boot_tau_summary , aes(x = tau_median)) +
  
  geom_rect(
    data = ci_band,
    inherit.aes = FALSE,
    aes(
      xmin = ci_lcl,
      xmax = ci_ucl,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "skyblue",
    alpha = 0.3
  ) +
  
  geom_density(
    
    alpha = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "red"
  ) +
  
  facet_wrap( ~ metric, scales = "free_y",
              labeller= labeller(metric= label_wrap_gen(width= 20))) +
  
  theme_minimal() +
  
  labs(
    x = "Median Mann-Kendall Tau",
    y = "Density"
  )

ggsave(plot= flow_boot, filename = "streamflow_boot.png", dpi= 600, width= 4.5, height= 3.7, bg= "white", units= "in")
##----------version for hydroclimate---------------------------------------######


boot_hydro_tau_summary <- hydro_all_boot %>%
  
  group_by(
    monitoring_location_id,
    metric
  ) %>%
  summarise(
    tau_median = median(mk_tau, na.rm = TRUE),
    tau_lcl = quantile(mk_tau, 0.05, na.rm = TRUE), #0.025 and 0.975 for 95% conf interval
    tau_ucl = quantile(mk_tau, 0.95, na.rm = TRUE),
    .groups = "drop"
  )



boot_hydro_tau_summary <- boot_hydro_tau_summary %>%
  mutate(
    metric = recode(
      metric,
      annual_prcp = "annual precip",
      fall_airtemp = "avg. fall air temp",
      growing_deg = "growing degree days",
      summer_airtemp = "avg. summer air temp",
      total_ros= "total rain-on-snow",
      winter_prcp= "winter precip",
      ws_swe_ratio= "winter-spring swe (snow) ratio"
    )
  )

ci_band <- boot_hydro_tau_summary  %>%
  group_by(metric) %>%
  summarise(
    ci_lcl = median(tau_lcl, na.rm = TRUE),
    ci_ucl = median(tau_ucl, na.rm = TRUE),
    .groups = "drop"
  )


hydro_boot<- ggplot(boot_hydro_tau_summary , aes(x = tau_median)) +
  
  geom_rect(
    data = ci_band,
    inherit.aes = FALSE,
    aes(
      xmin = ci_lcl,
      xmax = ci_ucl,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "skyblue",
    alpha = 0.3
  ) +
  
  geom_density(
    
    alpha = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "red"
  ) +
  
  facet_wrap( ~ metric, scales = "free_y",
              labeller= labeller(metric= label_wrap_gen(width= 20))) +
  
  theme_minimal() +
  
  labs(
    x = "Median Mann-Kendall Tau",
    y = "Density"
  )

ggsave(plot= hydro_boot, filename = "hydroclimate_boot.png", dpi= 600, width= 4.5, height= 3.7, bg= "white", units= "in")