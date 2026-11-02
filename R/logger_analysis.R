source(here::here("R", "functions.R"))
install_and_load_all_packages()
install.packages("myClim")
library(myClim)

tms <- mc_read_files(here::here("Data", "Loggers Data"), dataformat_name = "TOMST", silent = TRUE)
tms_p <- mc_prep_crop(tms, start = as.POSIXct("2024-03-08", tz = "UTC"),
                             end   = as.POSIXct("2024-11-28", tz = "UTC"))
mc_plot_raster(tms_p, sensors = c("TMS_moist"))
mc_plot_raster(tms_p, sensors = c("TMS_T3"))
tms_day <- mc_agg(tms_p, fun = "mean", period = "day", min_coverage = 0.95)
tms_calc <- mc_calc_vwc(tms_day)
tms_long <- mc_reshape_long(tms_calc)

tms_long |>
  dplyr::filter(sensor_name == "TMS_moist_mean") |>
  ggplot(aes(x = datetime, y = value)) +
  geom_line(color = "steelblue") +
  facet_wrap(~ locality_id) +
  labs(x = "Date", y = "Soil Moisture") +
  theme_bw()

tms_long |>
  dplyr::filter(sensor_name == "TMS_T1_mean") |>
  ggplot(aes(x = datetime, y = value)) +
  geom_line(color = "steelblue") +
  facet_wrap(~ locality_id) +
  labs(x = "Date", y = "Soil Temperature") +
  theme_bw()

