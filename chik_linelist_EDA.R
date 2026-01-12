library(pacman)

pacman::p_load(tidyverse, ggplot2, lubridate, data.table,sf, viridis, plotly, readxl, janitor, gt, scales,fitdistrplus, MASS)

kenya_sf<-st_read("Kenya/gadm41_KEN_0.shp")
kenya_counties<-st_read("Kenya/gadm41_KEN_1.shp")
kenya_subcounties<-st_read("Kenya/gadm41_KEN_2.shp")

kenya_counties<-kenya_counties|>
  mutate(county=NAME_1)
kenya_outline <- kenya_counties |>
  group_by(NAME_1) |>
  mutate(geometry = st_union(geometry))

# EDA of the Chikungunya Line List ----------------------------------------

chik_linelist <- read_excel("Chikungunya linelistJuly_2025.xlsx")

chik_linelist <- chik_linelist|>
  clean_names()

#rename the columns and also the column contents
#rename the important columns
chik_linelist<-chik_linelist|>
  rename(status=in_patien_t_outpatient)
names(chik_linelist) <- sub("_yes_n[o0]$", "", names(chik_linelist))



chik_linelist <- chik_linelist |>
  mutate(
    status = str_to_lower(status), # lowercase
    status=str_trim(status),                                          # remove spaces
    status=case_when(
      str_detect(status, "in") ~ "inpatient",                      # any variation of "in"
      str_detect(status, "out|op") ~ "outpatient",      # any variation of "out", "op", "outp"
      TRUE ~ NA_character_                                           # anything else → NA
    )
  )

#group the age_years into category for easier analysis

chik_linelist <- chik_linelist |>
  mutate(
    age_years = as.numeric(age_years),  # convert to numeric
    age_cat = case_when(
      age_years < 5 ~ "0-4",
      age_years >= 5 & age_years < 10 ~ "5-9",
      age_years >= 10 & age_years < 15 ~ "10-14",
      age_years >= 15 & age_years < 20 ~ "15-19",
      age_years >= 20 & age_years < 25 ~ "20-24",
      age_years >= 25 & age_years < 30 ~ "25-29",
      age_years >= 30 & age_years < 35 ~ "30-34",
      age_years >= 35 & age_years < 40 ~ "35-39",
      age_years >= 40 & age_years < 45 ~ "40-44",
      age_years >= 45 & age_years < 50 ~ "45-49",
      age_years >= 50 & age_years < 55 ~ "50-54",
      age_years >= 55 & age_years < 60 ~ "55-59",
      age_years >= 60 & age_years < 65 ~ "60-64",
      age_years >= 65 ~ "65+",
      TRUE ~ NA_character_
    )
  )

chik_linelist<-chik_linelist|>
  mutate(county=case_when(
    sub_county=="Likoni"~"Mombasa",
    TRUE~ county
  ))

#create the age-sex pyramid
# Prepare data
pyramid_data <- chik_linelist |>
  mutate(
    age_group2 = if_else(age_years < 60, "Under 60", "Above 60")
  ) |>
  select(age_years, age_group2, sex) |>
  group_by(sex, age_group2) |>
  summarise(total = n(), .groups = "drop") |>
  mutate(
    total = ifelse(tolower(sex) == "male", -total, total)  # males on the left
  )

median_age <- median(chik_linelist$age_years, na.rm = TRUE)
median_age

ggplot(pyramid_data|>
         filter(!is.na(age_group2)), aes(x = age_group2, y = total, fill = sex)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = abs, name = "Number of Cases") +
  scale_fill_manual(values = c("Male" = "#0072B2", "Female" = "#D55E00")) +
  labs(
    title = "Age–Sex Distribution for Chikungunya Cases",
    x = "Age Group",
    fill = "Sex"
  ) +
  annotate(
    "text",
    x = "Under 60",
    y = 0,
    label = paste("Median age =", round(median_age, 1), "years"),
    size = 4,
    fontface = "bold",
    hjust = -1
  ) +
  theme_minimal()


ggplot(pyramid_data, aes(x = age_group2, y = total, fill = sex)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = abs, name = "Number of Cases") +
  scale_fill_manual(values = c("Male" = "#0072B2", "Female" = "#D55E00")) +
  labs(
    title = "Age–Sex Distribution for Chikungunya Cases",
    x = "Age Category",
    fill = "Sex"
  ) +
  annotate(
    "text",
    x = "Under 60",
    y = 0,
    label = paste("Median age =", round(median_age, 1), "years"),
    color = "black", size = 4, fontface = "bold", hjust = -1
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# ggplot(pyramid_data, aes(x = age_cat, y = total, fill = sex)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   coord_flip() +
#   scale_y_continuous(labels = abs, name = "Number of Cases") +
#   scale_fill_manual(values = c("Male" = "#0072B2", "Female" = "#D55E00")) +
#   labs(
#     title = "Age–Sex Pyramid for Chikungunya Cases",
#     x = "Age Group",
#     fill = "Sex"
#   ) +
#   annotate(
#     "text", x = "30-34", y = 0, 
#     label = paste("Median age =", round(median_age, 1), "years"),
#     color = "black", size = 4, fontface = "bold", hjust = -1
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
#   )

#how many cases per county
chik_linelist<-chik_linelist|>
  mutate(county=case_when(
    county=="MOMBASA"~"Mombasa",
    county=="mombasa"~"Mombasa",
    county=="msa"~"Tana River",
    TRUE~county
  ))

county_cases<-chik_linelist|>
  group_by(county)|>
  summarise(total=n(), .groups = "drop")

#how many cases per sub_county
chik_linelist<-chik_linelist|>
  mutate(sub_county=case_when(
    sub_county=="NYALI"~"Nyali",
    sub_county=="likoni"~"Likoni",
    sub_county=="likon"~"Likoni",
    TRUE~sub_county
  ))

sub_county_cases<-chik_linelist|>
  group_by(sub_county)|>
  summarise(total=n(), .groups = "drop")

sub_county_cases |>
  filter(!sub_county=="")|>
  ggplot(aes(x = total, y = reorder(sub_county, total))) +
  geom_col(fill = "maroon") +
  labs(
    title = "Total Cases by Subcounty",
    x = "Number of Cases",
    y = "Subcounty"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


# Status, sex, age_cat distributions
chik_linelist |>
  count(status) |>
  mutate(pct = n / sum(n)) |>
  gt()

chik_linelist |>
  count(sex) |>
  mutate(pct = n / sum(n)) |>
  gt()

chik_linelist |>
  count(age_cat) |>
  mutate(pct = n / sum(n)) |>
  arrange(desc(pct)) |>
  gt()

# Symptom prevalence table
symptom_cols <- c("fever","rash","pain_behind_the_eyes","vomiting","headache",
                  "joint_pains_arthritis","rapid_breathing","nausea",
                  "swollen_glands","severe_abdominal","chills")

symptom_prev <- chik_linelist |>
  pivot_longer(all_of(symptom_cols), names_to = "symptom", values_to = "present") |>
  count(symptom, present) |>
  group_by(symptom) |>
  mutate(pct = n / sum(n)) |>
  filter(present == "yes") |>
  arrange(desc(pct))

symptom_prev |> gt()

# Bar plot of symptom prevalence
symptom_prev |>
  ggplot(aes(x = reorder(symptom, pct), y = pct)) +
  geom_col(fill = "#2C7FB8") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Symptom", y = "Prevalence", title = "Symptom prevalence") +
  theme_minimal()

# Prevalence of each symptom by status
by_status <- chik_linelist |>
  pivot_longer(all_of(symptom_cols), names_to = "symptom", values_to = "present") |>
  filter(!is.na(status)) |>
  count(status, symptom, present) |>
  group_by(status, symptom) |>
  mutate(pct = n / sum(n)) |>
  filter(present == "yes")

by_status |>
  ggplot(aes(x = reorder(symptom, pct), y = pct, fill = status)) +
  geom_col(position = position_dodge(width = 0.7)) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Symptom", y = "Prevalence", title = "Symptom prevalence by status") +
  theme_minimal()

by_status_sex <- chik_linelist |>
  pivot_longer(all_of(symptom_cols), names_to = "symptom", values_to = "present") |>
  filter(!is.na(status), !is.na(sex)) |>
  count(status, sex, symptom, present) |>
  group_by(status, symptom) |>
  mutate(pct = n / sum(n)) |>
  ungroup() |>
  filter(present == "yes")

by_status_sex |>
  ggplot(aes(x = fct_reorder(symptom, pct, .fun = sum), y = pct, fill = sex)) +
  geom_col(position = position_stack()) +
  coord_flip() +
  #facet_wrap(~ status, nrow = 1) +
  scale_fill_manual(values = c("Male" = "gold", "Female" = "steelblue")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Symptom",
    y = "Prevalence among 'yes' responses",
    fill = "Sex",
    title = "Symptom Prevalence by Status and Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Status by sex and age_cat
chik_linelist |>
  filter(!is.na(sex), !is.na(status))|>
  count(sex, status) |>
  group_by(sex) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = sex, y = pct, fill = status)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Status distribution by sex", x = "Sex", y = "Proportion") +
  theme_minimal()

chik_linelist |>
  filter(!is.na(status), status == "inpatient") |>
  count(age_cat) |>
  mutate(
    age_cat = factor(
      age_cat,
      levels = c(
        "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
        "30-34", "35-39", "40-44", "45-49", "50-54",
        "55-59", "60-64", "65+"
      ),
      ordered = TRUE
    )
  ) |>
  ggplot(aes(x = age_cat, y = n, fill = "inpatient")) +
  geom_col(show.legend = FALSE) +
  labs(
    title = "Number of Inpatients by Age Category",
    x = "Age Category",
    y = "Number of Inpatients"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

chik_linelist |>
  filter(!is.na(status), !status == "inpatient") |>
  count(age_cat) |>
  mutate(
    age_cat = factor(
      age_cat,
      levels = c(
        "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
        "30-34", "35-39", "40-44", "45-49", "50-54",
        "55-59", "60-64", "65+"
      ),
      ordered = TRUE
    )
  ) |>
  ggplot(aes(x = age_cat, y = n, fill = "inpatient")) +
  geom_col(show.legend = FALSE) +
  labs(
    title = "Number of Outpatient by Age Category",
    x = "Age Category",
    y = "Number of Outpatients"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#when were the patients diagnosed
#convert the confirmation date to date

# Define a function to clean date formats
clean_date_column <- function(date_col) {
  date_col <- as.character(date_col)
  
  cleaned <- case_when(
    # Excel serial number (e.g. 45814)
    grepl("^[0-9]+$", date_col) ~ as.character(as.Date(as.numeric(date_col), origin = "1899-12-30")),
    
    # Text month (e.g. "JUNE 2, 2025" or "June 2, 2025")
    grepl("[A-Za-z]", date_col) ~ as.character(parse_date_time(date_col, orders = c("b d, Y", "B d, Y"))),
    
    # Dot-separated (e.g. "3.5.2025")
    grepl("\\.", date_col) ~ as.character(dmy(date_col, quiet = TRUE)),
    
    # Slash-separated (e.g. "29/5/2025")
    grepl("/", date_col) ~ as.character(dmy(date_col, quiet = TRUE)),
    
    TRUE ~ NA_character_
  )
  
  as.Date(cleaned)
}

# Apply the function to all relevant columns
chik_linelist <- chik_linelist %>%
  mutate(
    date_seen_clean = clean_date_column(date_seen_ddmmyyyy),
    onset_date_clean = clean_date_column(date_of_onset_of_illness_ddmmyyyy),
    confirmation_date_clean = clean_date_column(confirmation_date)
  )

chik_linelist <- chik_linelist %>%
  mutate(
    diff_days = as.numeric(date_seen_clean - onset_date_clean)
  )

#get the monthly data 
monthly_data<-chik_linelist|>
  group_by(date_seen_clean)|>
  summarise(total=n(), .groups="drop")

ggplot(monthly_data, aes(x = date_seen_clean, y = total)) +
  geom_line() +
  labs(title = "Monthly Trend", x = "Month", y = "Value")

ggplot(monthly_data, aes(x = factor(format(date_seen_clean, "%b"), levels = month.abb), y = total)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Monthly Totals",
    x = "Month",
    y = "Total"
  ) +
  theme_minimal()


#Create a map to show clustering in the mombasa subcounties
subcounty_data<-chik_linelist|>
  filter(!is.na(sub_county))|>
  group_by(county,sub_county)|>
  summarise(total_cases=n(), .groups="drop")

subcounty_map_data<-kenya_subcounties|>
  left_join(subcounty_data, by=c("NAME_2"="sub_county"))|>
  filter(NAME_1 =="Mombasa")

map_plot <- ggplot(subcounty_map_data) +
  geom_sf(
    aes(fill = total_cases),
    color = "white",
    size = 0.2
  ) +
  geom_sf_text(
    data = subcounty_map_data,
    aes(label = NAME_2),
    color = "black",
    size = 3,
    fontface = "bold"
  )+
  scale_fill_viridis_c(option = "plasma", name = "Total Cases", direction=-1) +
  labs(
    title = "Chikungunya Cases by Subcounty - Mombasa County"
  ) +
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

map_plot

#calculating the suspected and confirmed cases
sus_conf_cases<- chik_linelist|>
  mutate(
    case_status = case_when(
      tolower(specimen_taken) == "yes" ~ "Confirmed",
      TRUE ~ "Suspected"
    )
  ) |>
  group_by(case_status) |>
  summarise(total_cases = n(), .groups="drop")

sus_conf_cases|>
  ggplot(aes(x=case_status, y=total_cases, fill=case_status))+
  geom_col()+
  scale_fill_manual(values = c("Confirmed" = "gold", "Suspected" = "steelblue"))
  labs(
    title="Confirmed vs Suspected cases",
    x="Case Status",
    y="Total Cases"
  )+
  theme_minimal()
  
  #for the confirmed cases and suspected, how many were below 60 and above 60 years
  age_summary <- chik_linelist |>
    mutate(
      case_type = case_when(
        tolower(specimen_taken) == "yes" ~ "Confirmed",
        tolower(specimen_taken) == "no" | specimen_taken == "" | is.na(specimen_taken) ~ "Suspected",
        TRUE ~ "Unknown"
      ),
      age_group2 = if_else(age_years < 60, "Under 60", "Above 60")
    ) |>
    filter(case_type != "Unknown") |>
    group_by(case_type, age_group2) |>
    summarise(total_cases = n(), .groups = "drop")
  
  ggplot(age_summary|>
           filter(!is.na(age_group2)), aes(x = age_group2, y = total_cases, fill = case_type)) +
    geom_col() +
    scale_fill_manual(values = c("Confirmed" = "gold", "Suspected" = "steelblue"))+
    labs(
      title = "Distribution of Chikungunya Cases by Age Group",
      x = "Age Group",
      y = "Number of Cases",
      fill = "Case Type"
    ) +
    theme_minimal(base_size = 14)

  #what are the symptoms prevalent in either the suspected or the confirmed cases
  symptom_by_case_type <- chik_linelist |>
    mutate(
      case_type = case_when(
        tolower(specimen_taken) == "yes" ~ "Confirmed",
        tolower(specimen_taken) == "no" | specimen_taken == "" | is.na(specimen_taken) ~ "Suspected",
        TRUE ~ "Unknown"
      )
    ) |>
    filter(case_type != "Unknown") |>
    pivot_longer(
      all_of(symptom_cols),
      names_to = "symptom",
      values_to = "present"
    ) |>
    count(case_type, age_cat, symptom, present) |>
    group_by(case_type, age_cat, symptom) |>
    mutate(pct = n / sum(n)) |>
    filter(present == "yes")
  
  
  #visualizing symptom prevalence by case type
ggplot(symptom_by_case_type, aes(x = fct_reorder(symptom, pct), y = pct, fill = case_type)) +
    geom_col(position = position_dodge(width = 0.7)) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = c("Confirmed" = "gold", "Suspected" = "steelblue"))+
    labs(
      x = "Symptom",
      y = "Prevalence",
      fill = "Case Type",
      title = "Symptom Prevalence by Case Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )  

mombasa_pop <- pop_data |>
  filter(county == "Mombasa" & year==2025) |>
  mutate(
    # Extract the lower bound of the age group
    age_lower = as.numeric(str_extract(Age, "^\\d+")),
    
    # Create the new category
    age_group2 = if_else(age_lower < 60, "Under 60", "Above 60")
  )|>
  group_by(age_group2) |>
  summarise(total_population = sum(population, na.rm = TRUE), .groups = "drop")

#from the chiklinelist, find the proportion above 60 and below 60
chik_age_dist <- chik_linelist |>
  mutate(
    age_group2 = if_else(age_years < 60, "Under 60", "Above 60")
  ) |>
  group_by(age_group2, status) |>
  summarise(total_cases = n(), .groups = "drop")|>
  mutate(proportion = round(total_cases / sum(total_cases),3))

symptom_by_case_type1 <- chik_linelist |>
  mutate(
    case_type = case_when(
      tolower(specimen_taken) == "yes" ~ "Confirmed",
      tolower(specimen_taken) == "no" | specimen_taken == "" | is.na(specimen_taken) ~ "Suspected",
      TRUE ~ "Unknown"
    )
  ) |>
  group_by(case_type, age_cat, other_symptoms_specify)|>
  summarise(cases=n(), na.rm=TRUE,
            .groups="drop")


clean_dates <- function(x) {
  
  # Step 1: Excel numeric dates
  x_numeric <- suppressWarnings(as.numeric(x))
  
  excel_dates <- dplyr::if_else(
    !is.na(x_numeric),
    as.Date(x_numeric, origin = "1899-12-30"),
    as.Date(NA)     # must be Date NA
  )
  
  # Step 2: Parse character formats
  char_dates_raw <- parse_date_time(
    x,
    orders = c("d.m.y", "d/m/y", "B d, Y"),
    exact = FALSE
  )
  
  char_dates <- as.Date(char_dates_raw)
  
  # Step 3: Combine
  final <- dplyr::coalesce(excel_dates, char_dates)
  
  return(final)
}



chik_linelist <- chik_linelist |>
  mutate(
    date_seen_clean = clean_dates(date_seen_ddmmyyyy),
    onset_date_clean = clean_dates(date_of_onset_of_illness_ddmmyyyy),
    confirmation_date_clean = clean_dates(confirmation_date)
  )


#using the date_seen_clean, count the number of cases per day
daily_cases <- chik_linelist |>
  group_by(date_seen_clean) |>
  summarise(daily_total = n(), .groups = "drop")

#find the diffeerence between date seen and onset date
chik_linelist <- chik_linelist |>
  mutate(
    diff_days = as.numeric(confirmation_date_clean - date_seen_clean)
  )

#if the patient was inpatient and had joint_pain_arthritis, they are severe, the rest are mild. group them into >60 and below 60 and find the proprtion of mild and severe
chik_linelist <- chik_linelist |>
  mutate(
    severity = case_when(
      status == "inpatient" & joint_pains_arthritis == "yes" ~ "Severe",
      TRUE ~ "Mild"
    ),
    age_group2 = if_else(age_years < 60, "Under 60", "Above 60")
  )

#now find the proportion of mild and severe cases by age group
severity_by_age <- chik_linelist |>
  group_by(age_group2, severity) |>
  summarise(total_cases = n(), .groups = "drop") |>
  mutate(proportion = round(total_cases / sum(total_cases),3))|>
  gt()

severity_by_age
chik_clean <- chik_linelist |>
  filter(!is.na(confirmation_date_clean))
#plot the epi curve using confirmation_date_clean
ggplot(chik_clean, aes(x = confirmation_date_clean)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "steelblue") +
  labs(
    title = "Epi Curve: Chikungunya Outbreak",
    x = "Date of Confirmation",
    y = "Number of Cases"
  ) +
  theme_minimal()
#find the week with the highest number of cases
peak_week <- chik_linelist |>
  mutate(
    confirmation_date_clean = as.Date(confirmation_date_clean),
    week = isoweek(confirmation_date_clean),
    year = year(confirmation_date_clean)
  ) |>
  filter(!is.na(confirmation_date_clean)) |>
  count(year, week, name = "cases") |>
  slice_max(cases)

peak_week

#find the peak week according to the date_seen
peak_week_seen <- chik_linelist |>
  mutate(
    date_seen_clean = as.Date(date_seen_clean),
    week = isoweek(date_seen_clean),
    year = year(date_seen_clean)
  ) |>
  filter(!is.na(date_seen_clean)) |>
  count(year, week, name = "cases") 
peak_week_seen

peak_week_subcounty <- chik_linelist |>
  mutate(
    date_seen_clean = as.Date(date_seen_clean),
    week = isoweek(date_seen_clean),
    year = year(date_seen_clean)
  ) |>
  filter(
    year == peak_week$year,
    week == peak_week$week
  ) |>
  count(subcounty, name = "cases") |>
  arrange(desc(cases))

peak_week_subcounty


#plot the epi curve using seen_date_clean
ggplot(chik_clean, aes(x = date_seen_clean)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "steelblue") +
  labs(
    title = "Epi Curve: Chikungunya Outbreak",
    x = "Date of Confirmation",
    y = "Number of Cases"
  ) +
  theme_minimal()

ggplot(data = chik_linelist, aes(x = age_years)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Age Distribution", x = "Age", y = "Count") +
  theme_minimal()

#filtered out the dates after 2025-05-02 since the cases were more uniform after that date
daily_cases <- chik_linelist |>
  mutate(date_seen_clean = as.Date(date_seen_clean)) |>               # Ensure date
  filter(!is.na(date_seen_clean)) |>                                  # Remove NA dates
  filter(date_seen_clean >= as.Date("2025-05-02")) |>                 # Keep only dates >= 2025-05-02
  group_by(date_seen_clean) |>                                        # Group by date
  summarise(cases = n(), .groups = "drop") |>                         # Count cases per date
  arrange(date_seen_clean) |>                                         # Order by date
  complete(date_seen_clean = seq(min(date_seen_clean, na.rm = TRUE), 
                                 max(date_seen_clean, na.rm = TRUE), 
                                 by = "1 day"), 
           fill = list(cases = 0)) |>                                 # Fill missing dates
  mutate(day = as.integer(date_seen_clean - min(date_seen_clean, na.rm = TRUE)) + 1)  # Day number



# Plot
ggplot(daily_cases, aes(x = day, y = cases)) +
  geom_line() +
  labs(
    title = "Cases Over Time",
    x = "Day",
    y = "Number of Cases"
  ) +
  theme_minimal()

#Find the best distribution that represents our data
cases <- daily_cases$cases
cases <- cases[!is.na(cases)]

#since the daily cases are discrete, find the best discrete distribution
#fit poisson
fit_pois <- fitdist(cases, "pois")
summary(fit_pois)
#fit negative binomial
fit_nb <- fitdist(cases, "nbinom")
summary(fit_nb)
#fit geometric
fit_geom <- fitdist(cases, "geom")
summary(fit_geom)

#fit norm
fit_norm <- fitdist(cases, "norm")
summary(fit_norm)
#compare the distributions using AIC
aic_values <- data.frame(
  Distribution = c("Poisson", "Negative Binomial", "Geometric", "Normal"),
  AIC = c(fit_pois$aic, fit_nb$aic, fit_geom$aic, fit_norm$aic)
)
aic_values <- aic_values |>
  arrange(AIC)

gofstat(list(fit_pois, fit_nb, fit_geom, fit_norm))


