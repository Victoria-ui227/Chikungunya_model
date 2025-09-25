#Before building the Chikungunya model, are we able to understand the data and what it says

library(pacman)

pacman::p_load(tidyverse, ggplot2, lubridate, data.table,sf, viridis, plotly)

chik_data<-fread("data/CHIKV_data_KHIS.csv")

fever_per_county_data<-fread("data/KHIS_data_fever_county.csv")

fever_data<-fread("data/KHIS_data_fever.csv")

kenya_counties<-st_read("Kenya/gadm41_KEN_1.shp")

kenya_counties<-kenya_counties|>
  mutate(county=NAME_1)
kenya_outline<-kenya_counties|>
  mutate(geometry=st_union(geometry))

# Analyzing Chik data -----------------------------------------------------

#data cleaning

chik_data<-chik_data|>
  rename(chik_cases = `MOH 705A Rev 2020_ Chikungunya`)

chik_data <- chik_data |>
  separate(periodname, into = c("month", "year"), sep = " ")

chik_data$month = factor(chik_data$month, levels = month.name)

chik_data <- chik_data |>
  rename(county = organisationunitname) |>
  mutate(
    county = case_when(
      county == "Elgeyo Marakwet County" ~ "Elgeyo-Marakwet",
      TRUE ~ str_remove(county, " County$")
    )
  )

  
#from the data given what is the monthly distribution of chik cases
monthly_data<-chik_data|>
  group_by(month,year)|>
  summarise(total = sum(chik_cases, na.rm = TRUE), .groups = "drop")


monthly_data |>
  ggplot(aes(x = month, y = total, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  facet_wrap(~ year, ncol = 2, scales = "free_y") +
  labs(
    title = "Monthly Cases per Year",
    x = "Month",
    y = "Total Cases"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # tilt month labels
  )

#which months per year have the highest cases
top_months <- monthly_data %>%
  group_by(year) %>%
  filter(total == max(total, na.rm = TRUE)) %>%
  ungroup()

top_months

monthly_data %>%
  ggplot(aes(x = month, y = total)) +
  geom_col(show.legend = FALSE, alpha = 0.7) +
  # highlight the max months
  geom_col(data = top_months, aes(x = month, y = total),
           fill = "tomato", show.legend = FALSE) +
  geom_text(data = top_months,
            aes(label = total, y = total + 0.5),
            color = "black", size = 3, fontface = "bold") +
  facet_wrap(~ year, ncol = 2, scales = "free_y") +
  labs(
    title = "Monthly Cases per Year (Highlighted Highest Months)",
    x = "Month",
    y = "Total Cases"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#which year has the highest cases
yearly_cases <- monthly_data |>
  group_by(year) |>
  summarise(total = sum(total, na.rm = TRUE), .groups = "drop")

yearly_cases |>
  ggplot(aes(x = factor(year), y = total, fill="orange")) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = total),
            vjust = -0.5,   # move labels above the bars
            size = 4) +
  labs(
    title = "Total Cases per Year",
    x = "Year",
    y = "Total Cases"
  ) +
  theme_minimal()

  

#which counties have the highest cases for all years
county_cases<-chik_data|>
  group_by(county)|>
  summarise(total = sum(chik_cases, na.rm = TRUE), .groups = "drop")


county_plot <- ggplot(county_cases, aes(x = reorder(county, total), y = total)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = total),
            vjust = -0.2,   # use hjust instead of vjust for coord_flip
            size = 4) +
  labs(
    title = "Chikungunya Cases per County (2020–2025)",
    x = "County",
    y = "Total Cases"
  ) +
  theme_minimal()


county_plot
map_data<-kenya_counties|>
  left_join(county_cases, by="county")

map_plot<-ggplot(map_data)+
  geom_sf(aes(fill=total,text = paste0("County: ", county, "<br>Cases: ", total)),size=0.2)+
  scale_fill_viridis_c(
    option="C",
    direction = -1,
    na.value="grey90",
    name=NULL
  )+
  labs(title="Distribution of Chikungunya cases from 2020-2025 per county")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text =element_blank(),
        axis.ticks = element_blank())
ggplotly(map_plot, tooltip = "text")

#distribution of cases per year per county
cases_per_county<-chik_data|>
  group_by(county, year)|>
  summarise(total = sum(chik_cases, na.rm = TRUE), .groups = "drop")

top_counties <- cases_per_county %>%
  group_by(year) %>%
  slice_max(order_by = total, n = 3, with_ties = FALSE) %>%
  ungroup()


cases_per_county %>%
  ggplot(aes(x = county, y = total)) +
  geom_col(data = top_counties, aes(x = county, y = total),
           fill = "tomato", show.legend = FALSE) +
  geom_text(data = top_counties,
            aes(label = total, y = total +1),
            color = "black", size = 3) +
  geom_col(show.legend = FALSE, alpha = 0.7) +
  facet_wrap(~ year, ncol = 2, scales = "free_y") +
  labs(
    title = "County Cases per Year",
    x = "county",
    y = "Total Cases"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Analyzing the county and subcounty data ---------------------------------
fever_data<-fever_data|>
  rename(haemorrhagic_fevers_total=`IDSR Other Viral Haemorrhagic Fevers Total`,
         over_5_cases_fever_less_7_days=`MOH 515 Rev 2020_Number of over 5 years cases with fever for less than 7 days`,
         under_5_fever_less_7_days=`MOH 515 Rev 2020_Number of under 5 years cases with fever for less than 7 days`,
         neurological_disorder=`MOH_744A_No. of Neurological disorders`,
         location=organisationunitname)

fever_data <- fever_data |>
  separate(periodname, into = c("month", "year"), sep = " ")

fever_data$month = factor(fever_data$month, levels = month.name)

#how many people over the age of 5 had fever for less than 7 days
sum(fever_data$over_5_cases_fever_less_7_days, na.rm = TRUE)

#how many people in the country under five had fever for more than 5 days

sum(fever_data$under_5_fever_less_7_days, na.rm = TRUE)

#what are the fever cases per month
monthly_fever<-fever_data|>
  group_by(month, year)|>
  summarise(total_under_5=sum(under_5_fever_less_7_days, na.rm=TRUE),
            total_over_5=sum(over_5_cases_fever_less_7_days, na.rm = TRUE),
            .groups = "drop")

sub_county_cases<-fever_data|>
  filter(str_detect(str_trim(location), regex("Sub\\W*County\\s*$", ignore_case = TRUE)))|>
  group_by(location)|>
  summarise(total_under_5=sum(under_5_fever_less_7_days, na.rm=TRUE),
            total_over_5=sum(over_5_cases_fever_less_7_days, na.rm = TRUE),
            total_neurological_disorders=sum(neurological_disorder, na.rm=TRUE),
            .groups = "drop")|>
  mutate(total_cases=total_under_5+total_over_5)

county_cases<-fever_data|>
  filter(
    str_detect(str_trim(location), regex("County\\s*$", ignore_case = TRUE)) &
      !str_detect(str_trim(location), regex("Sub\\W*County\\s*$", ignore_case = TRUE))
  )|>
  group_by(year,location)|>
  summarise(total_under_5=sum(under_5_fever_less_7_days, na.rm=TRUE),
            total_over_5=sum(over_5_cases_fever_less_7_days, na.rm = TRUE),
            total_neurological_disorders=sum(neurological_disorder, na.rm=TRUE),
            .groups = "drop")|>
  mutate(total_fevers=total_under_5+total_over_5)


#what are the counties that contributed to Chikungunya cases from May 2021 to May 2025

contributing_counties<-chik_data|>
  filter(periodcode >= 202105 & periodcode <= 202505)|>
  group_by(year,county)|>
  summarise(infected=sum(chik_cases, na.rm=TRUE), .groups = "drop")

#what is the population of under 5 and over 5 in these counties that have had Chikungunya in the past 4 years

pop_data<-fread("Population per county,sex,age group.csv")

#population of under 5 and over 5 in each county in each year
over_under5_pop <- pop_data |>
  group_by(year, county) |>
  summarise(
    under_5 = sum(population[Age %in% c("0", "0-4", "1-4")], na.rm = TRUE),
    over_5  = sum(population[!Age %in% c("0", "0-4", "1-4")], na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    county = str_replace_all(county, "-\\s+", "-"),  # fix space after hyphen
    county = case_when(
      county == "Taita-Taveta" ~ "Taita Taveta",
      county == "Nairobi City" ~ "Nairobi",
      TRUE ~ county
    )
  )


contributing_counties<-contributing_counties|>
  mutate(year=as.integer(year))|>
  left_join(over_under5_pop, by=c("year","county"))


#the current county cases do not have data from 2021 that is contributing to our intial conditions and therefore we have to use the fever_per_county data

fever_per_county_data<-fever_per_county_data|>
  rename(over_5_cases_fever_less_7_days=`MOH 515 Rev 2020_Number of over 5 years cases with fever for less than 7 days`,
         under_5_fever_less_7_days=`MOH 515 Rev 2020_Number of under 5 years cases with fever for less than 7 days`,
         location=organisationunitname)
fever_per_county_data <- fever_per_county_data |>
  mutate(
    year  = floor(periodcode / 100),          # extract year
    month_num = periodcode %% 100,           # extract month number
    month = month.name[month_num]            # convert to month name
  ) |>
  select(periodcode, year, month, everything()) 

fever_per_county_data$month = factor(fever_per_county_data$month, levels = month.name)

#calculate the number of fever cases in each year and each county
fever_county <- fever_per_county_data |>
  filter(
    str_detect(str_trim(location), regex("County\\s*$", ignore_case = TRUE))
  ) |>
  group_by(year, location) |>
  summarise(
    fever_under_5 = sum(under_5_fever_less_7_days, na.rm = TRUE),
    fever_over_5  = sum(over_5_cases_fever_less_7_days, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    # remove "County" at the end
    county = str_remove(location, regex("\\s*County$", ignore_case = TRUE)),
    # clean spaces and standardize Elgeyo Marakwet
    county = str_replace_all(county, "\\s+", " "), 
    county = case_when(
      county == "Elgeyo Marakwet" ~ "Elgeyo-Marakwet",
      TRUE ~ county
    )
  ) |>
  select(year, county, fever_under_5, fever_over_5)


#join the fever_county to the contributing counties
contributing_counties<-contributing_counties|>
  left_join(fever_county|>mutate(year=as.integer(year)),
            by=c("year","county"))
  

#now with all the data we have, generate the intial conditions
intial_conditions<-contributing_counties|>
  filter(year==2021)|>
  summarise(total_infected=sum(infected, na.rm=TRUE),
            total_under5=sum(under_5, na.rm=TRUE),
            total_over5=sum(over_5, na.rm=TRUE),
            total_fever_under5=sum(fever_under_5, na.rm=TRUE),
            total_fever_over5=sum(fever_over_5, na.rm=TRUE),
            .groups="drop" )
