library(pacman)

pacman::p_load(
  dplyr,
  ggplot2,
  tidyr,
  readr,
  stringr,
  lubridate,
  forcats,
  zoo
)

#load data
climate_data<-read.csv("temp_precip2.csv")

#filter out for mombasa from 2020 to 2024
mombasa_climate<-climate_data|>
  filter(county=="Mombasa" & year>=2020 & year<=2024)|>
  mutate(
    month_name = case_when(
      is.na(month_name) | month_name == "" ~ "Dec",
      TRUE ~ month_name
    ))

#plot the total_precipitation against month_name
ggplot(mombasa_climate, aes(x=fct_relevel(month_name, month.abb), y=total_precip, group=year, color=factor(year)))+
  geom_line(size=1)+
  geom_point(size=2)+
  labs(
    title="Total Monthly Precipitation in Mombasa (2020-2024)",
    x="Month",
    y="Total Precipitation (mm)",
    color="Year"
  )+
  theme_minimal()

#Modifying the mombasa climate data to suit the Gillepsie model
rainfall_df <- mombasa_climate %>%
  arrange(year, month) %>%
  mutate(
    month_index = row_number()
  ) %>%
  select(month_index, year, month, total_precip)

#Extract rain data as a vector
rain <- rainfall_df$total_precip

#Standardize the rainfall by using the mean
rain_scaled <- rain / mean(rain)

#smooth rainfall data without creating noise
rain_smooth <- zoo::rollmean(rain_scaled, k = 3, fill = "extend")


temp_df<-mombasa_climate|>
  arrange(year, month) %>%
  mutate(
    month_index = row_number()
  ) %>%
  select(month_index, year, month, mean_temp)
