#########################################################################
# QBIO7008- ANALYSING THE EFFECT OF ENV. VARIABLES ON FISH GROWTH RATE  #
#########################################################################

## This project aims to investigate the complex relationships between hydrological
## and environmental factors, and the growth rates of three key lotic fish species, 
## Golden Perch (Macquaria ambigua), Bony Bream (Nematalosa erebi), and 
## Common Carp (Cyprinus carpio), within South and West Queensland's dryland river systems.
## Utilizing otolith-derived incremental growth data, along with river flow metrics
## and annual temperature readings, this research delves into how spatial and 
## temporal variations in streamflow and thermal regimes influence fish growth. 


#########################################################################

## Installing and loading required packages

# Equatiomatic has been removed from the CRAN repository. 
# Install specific version of equatiomatic from archive. First download onto PC.
# install.packages("C:\\Users\\hawwa\\Downloads\\equatiomatic_0.3.1.tar.gz", repos = NULL, type="source")
library(equatiomatic) 

# Loading all other packages
library(car)
library(readxl)
library(loo)
library(caret)
library(tidyverse)
library(DHARMa)
library(caret)
library(broom.mixed)
library(corrplot)
library(ggplot2)
library(gridExtra)
library(brms)
library(lubridate)
library(data.table)
library(ggrepel)
library(glmmTMB)
library(mgcViz)

#library(LDdiag) 
#also removed from CRAN repository. devtools required to install it from the archives
#devtools::install_version("LDdiag","0.1")
library(LDdiag)



#########################################################################
## Importing and cleaning data from La Trobe U's Microsoft Power BI Solution
#########################################################################

## Importing data downloaded from Microsoft Power BI
d1 <- read_excel("./Data/Bonybream.xlsx")
d2 <- read_excel("./Data/Carp.xlsx")
d3 <- read_excel("./Data/Goldenperch.xlsx")
Gauges <- read_excel("./Data/Gauges.xlsx")


# Checking the general structure of the datasets
head(d1)
head(d2)
head(d3)
head(Gauges)

## Merging species datasets to form a combined dataset for all 3 species
comb_data <- rbind(d1, d2, d3)

## Incorporating gauge ID data from Power BI into the species dataset (comb_data)

# Checking the number of unique GaugeID entries
length(unique(Gauges$GaugeID))
# 10 entries. So 10 gauges for 11 sites
# There are 2 gauges for Site 7. 
# Excluding records with GaugeID '422209A' , the second gauge
Gauges <- Gauges %>%
  filter(!(GaugeID == "422209A" & SiteID == 7))

# Checking for columns which can be used in merge
head(Gauges)
head(comb_data)

# Merging the comb_data with Gauges dataset based on SiteID
full_data <- comb_data %>%
  left_join(Gauges, by = "SiteID")

# Checking the structure and the first few rows of the merged dataset
str(full_data)
head(full_data)

# Modifying FishID to remove the last three characters (year identifier)
# In case the variable is needed for comparison between years
full_data <- full_data %>%
  mutate(FishID_normalized = str_sub(FishID, 1, -4))

# Adding year in which each fish would have been 2 years old, 1 year old and spawned
full_data <- full_data %>%
  mutate(
    Year_2yo = Year - Age + 2,
    Year_1yo = Year - Age + 1,
    Year_Spawning = Year - Age
  )

# Adding a Species column based on FishID_normalized prefixes
full_data <- full_data %>%
  mutate(Species = case_when(
    grepl("BB", substr(FishID_normalized, 3, 4), ignore.case = TRUE) ~ "BonyBream",
    grepl("GP", substr(FishID_normalized, 3, 4), ignore.case = TRUE) ~ "GoldenPerch",
    grepl("CC", substr(FishID_normalized, 3, 4), ignore.case = TRUE) ~ "CommonCarp",
    TRUE ~ "Unknown"  # A fallback option if no match is found
  ))


#########################################################################
## Incorporating environmental data from WMIP into full_data
#########################################################################

# Loading required packages and setting directory to folder containing WMIP data
library(readr)
library(dplyr)
library(lubridate)

WMIP_directory <- "./Data/WMIP/"

# Listing all filenames within directory
file_names <- list.files(path=WMIP_directory, pattern = "\\.csv$", full.names = TRUE)

# Creating a list to store data frames
data_frames <- list()

# Loop that goes through each file in WMIP directory
for (file_path in file_names) {
  gauge_id <- gsub(".*/|_flow\\.csv$", "", file_path)
  
  # Reading CSV files, skipping the first 3 rows, and reading all columns as character
  df <- read_csv(file_path, skip = 3, col_types = cols(.default = col_character()))
  
  # Renaming specific columns based on their positions
  colnames(df)[1] <- "Date"
  colnames(df)[2] <- "Mean_Level_m"
  colnames(df)[4] <- "Min_Level_m"
  colnames(df)[6] <- "Max_Level_m"
  colnames(df)[8] <- "Mean_Discharge_Ml"
  colnames(df)[10] <- "Min__Discharge_Ml"
  colnames(df)[12] <- "Max__Discharge_Ml"
  colnames(df)[14] <- "Max_Vol_Ml"
  
  # Converting necessary columns to appropriate types
  df <- df %>%
    mutate(
      Mean_Level_m = as.numeric(Mean_Level_m),
      Min_Level_m = as.numeric(Min_Level_m),
      Max_Level_m = as.numeric(Max_Level_m),
      Mean_Discharge_Ml = as.numeric(Mean_Discharge_Ml),
      Min__Discharge_Ml = as.numeric(Min__Discharge_Ml),
      Max__Discharge_Ml = as.numeric(Max__Discharge_Ml),
      Max_Vol_Ml = as.numeric(Max_Vol_Ml),
      GaugeID = gauge_id
    ) %>%
    dplyr::select(Date, Mean_Level_m, Min_Level_m, Max_Level_m, Mean_Discharge_Ml, Min__Discharge_Ml, Max__Discharge_Ml, Max_Vol_Ml, GaugeID)
  
  data_frames[[gauge_id]] <- df
}


# Combine all data frames into one
comb_flow <- bind_rows(data_frames)

# Remove rows where all the data are NA
comb_flow <- comb_flow %>% filter(!if_all(everything(), is.na))
# Remove rows where Date column starts with "Data use licence:" or "Glossary of terms"
comb_flow <- comb_flow %>%
  filter(
    !str_detect(Date, "^Data use licence:"),
    !str_detect(Date, "^Glossary of terms")
  )

# Removinge rows where Date is NA
comb_flow <- comb_flow %>% filter(!is.na(Date))

# Remove "00:00" from the Date column and parse the date
# Also changing the name as a safety measure, in case changes need to be reversed on the spot
comb_flow1 <- comb_flow %>%
  mutate(
    Date = str_remove(Date, " 00:00$"),  # Remove the "00:00" part
    Date = parse_date_time(Date, orders = c("Ymd", "Ymd HMS", "dmy", "dmy HMS"))
  )

# Adding an 'Year' column, extracted from the Date column
comb_flow1 <- comb_flow1 %>%
  mutate(Year = year(Date))

# Display the combined data frame to verify the result
head(comb_flow1)

# Quick check to make sure data from all years have been incorporated
count_2021 <- comb_flow1 %>% filter(Year == 2021) %>% nrow()
print(count_2021)
# 8 gauges * 365 days = 2920 records for each year
# looks good


#########################################################################
## Calculating the longest no flow period and total no. of zero flow days

# Loading libraries
library(dplyr)
library(data.table)


# Summarising full_data dataset to make unique rows for each GaugeID
cease_to_flow <- full_data %>% 
  dplyr::select(GaugeID, Cease_to_flow_level) %>%
  distinct()

# Merging cease_to_flow information with `comb_flow1`
comb_flow1 <- comb_flow1 %>%
  left_join(cease_to_flow, by = "GaugeID")

# Calculating the number of zero flow days per year for each gauge
zero_flow_days <- comb_flow1 %>%
  group_by(GaugeID, Year) %>%
  summarize(
    ZeroFlowDays = sum(Mean_Level_m <= Cease_to_flow_level, na.rm = TRUE),
    .groups = 'drop'
  )

# Converting dataset to a data.table
setDT(comb_flow1)

# Adding a logical column to identify no-flow days
comb_flow1[, `:=`(No_flow = Mean_Level_m <= Cease_to_flow_level)]

# Calculating the longest no-flow spell by gauge and year
longest_no_flow_spell <- comb_flow1[, {
  no_flow_days <- .SD[No_flow == TRUE, .(Date)]
  if (nrow(no_flow_days) > 0) {
    no_flow_days[, diff_date := c(1, diff(Date))]
    no_flow_days[, spell_id := cumsum(diff_date > 1)]
    spell_lengths <- no_flow_days[, .(length = .N), by = spell_id][, max(length, na.rm = TRUE)]
  } else {
    spell_lengths = 0
  }
  .(LongestZeroFlow = as.double(spell_lengths))
}, by = .(GaugeID, Year)]

# Calculating zero flow days and longest no flow spell for Year_1yo
zero_flow_days_1yo <- comb_flow1 %>%
  group_by(GaugeID, Year) %>%
  summarize(
    ZeroFlowDays = sum(Mean_Level_m <= Cease_to_flow_level, na.rm = TRUE),
    .groups = 'drop'
  )

longest_no_flow_spell_1yo <- comb_flow1[, {
  no_flow_days <- .SD[No_flow == TRUE, .(Date)]
  if (nrow(no_flow_days) > 0) {
    no_flow_days[, diff_date := c(1, diff(Date))]
    no_flow_days[, spell_id := cumsum(diff_date > 1)]
    spell_lengths <- no_flow_days[, .(length = .N), by = spell_id][, max(length, na.rm = TRUE)]
  } else {
    spell_lengths = 0
  }
  .(LongestZeroFlow = as.double(spell_lengths))
}, by = .(GaugeID, Year)]

# Calculating zero flow days and longest no flow spell for Year_Spawning
zero_flow_days_spawning <- comb_flow1 %>%
  group_by(GaugeID, Year) %>%
  summarize(
    ZeroFlowDays = sum(Mean_Level_m <= Cease_to_flow_level, na.rm = TRUE),
    .groups = 'drop'
  )

longest_no_flow_spell_spawning <- comb_flow1[, {
  no_flow_days <- .SD[No_flow == TRUE, .(Date)]
  if (nrow(no_flow_days) > 0) {
    no_flow_days[, diff_date := c(1, diff(Date))]
    no_flow_days[, spell_id := cumsum(diff_date > 1)]
    spell_lengths <- no_flow_days[, .(length = .N), by = spell_id][, max(length, na.rm = TRUE)]
  } else {
    spell_lengths = 0
  }
  .(LongestZeroFlow = as.double(spell_lengths))
}, by = .(GaugeID, Year)]

# Merging zero flow days and longest no flow spell into full_data
full_data <- full_data %>%
  left_join(zero_flow_days, by = c("GaugeID" = "GaugeID", "Year_2yo" = "Year"), suffix = c("", "_2yo")) %>%
  left_join(longest_no_flow_spell, by = c("GaugeID" = "GaugeID", "Year_2yo" = "Year"), suffix = c("", "_2yo")) %>%
  left_join(zero_flow_days_1yo, by = c("GaugeID" = "GaugeID", "Year_1yo" = "Year"), suffix = c("", "_1yo")) %>%
  left_join(longest_no_flow_spell_1yo, by = c("GaugeID" = "GaugeID", "Year_1yo" = "Year"), suffix = c("", "_1yo")) %>%
  left_join(zero_flow_days_spawning, by = c("GaugeID" = "GaugeID", "Year_Spawning" = "Year"), suffix = c("", "_spawning")) %>%
  left_join(longest_no_flow_spell_spawning, by = c("GaugeID" = "GaugeID", "Year_Spawning" = "Year"), suffix = c("", "_spawning"))


# Summarizing Mean_Level_m and Mean_Discharge_Ml by year and gauge
comb_flow_summary <- comb_flow1 %>%
  group_by(GaugeID, Year) %>%
  summarize(
    mean_level = mean(Mean_Level_m, na.rm = TRUE),
    mean_discharge = mean(Mean_Discharge_Ml, na.rm = TRUE),
    .groups = 'drop'
  )

# Incorporating the summarized level and discharge data (at 2yo)
full_data <- full_data %>%
  left_join(comb_flow_summary, by = c("GaugeID" = "GaugeID", "Year_2yo" = "Year"), suffix = c("", "_2yo"))


# Display the first few rows of the updated full_data
head(full_data)







#########################################################################
## Incorporating environmental data from SILO into full_data
#########################################################################

# Loading packages and setting directory to folder containing SILO data
library(dplyr)
SILO_directory <- "./Data/SILO/"

# Listing all CSV files in the directory
file_names <- list.files(path = SILO_directory, pattern = "\\.csv$", full.names = TRUE)

# Creating a list to store data frames
data_frames <- list()

# Loop that goes through each file in SILO directory
for (file_path in file_names) {
  # Extracting gauge ID from the filename
  gauge_id <- gsub(".*/|\\.csv$", "", file_path)
  
  # Reading CSV files
  df <- read.csv(file_path)
  
  # Adding gauge ID column
  df$GaugeID <- gauge_id
  
  # Storing the dataframe in the list
  data_frames[[gauge_id]] <- df
}

# Combining all gauge data frames into one, excluding the 'metadata' column
comb_gauge <- bind_rows(data_frames) %>%
  dplyr::select(-metadata)
# Checking the resultant dataset
head(comb_gauge)
str(comb_gauge)

# Extracting data year from the YYYY.MM.DD column and mean temp col.
comb_gauge <- comb_gauge %>%
  mutate(
    Date = ymd(YYYY.MM.DD),
    Year = year(Date),
    mean_temp = (max_temp + min_temp) / 2
  )

# Calculate annual averages for evaporation and rain in comb_gauge
comb_gauge1 <- comb_gauge %>%
  group_by(GaugeID, Year) %>%
  summarize(
    mean_daily_rain = mean(daily_rain, na.rm = TRUE),
    mean_evap_pan = mean(evap_pan, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate annual averages for mean_temp in comb_gauge
mean_temp_annual <- comb_gauge %>%
  group_by(GaugeID, Year) %>%
  summarize(
    mean_temp = mean(mean_temp, na.rm = TRUE),
    .groups = 'drop'
  )

# Incorporating mean_temp for different life stages into the full_data dataset
full_data <- full_data %>%
  left_join(mean_temp_annual, by = c("GaugeID" = "GaugeID", "Year_2yo" = "Year"), suffix = c("", "_2yo")) %>%
  left_join(mean_temp_annual, by = c("GaugeID" = "GaugeID", "Year_1yo" = "Year"), suffix = c("", "_1yo")) %>%
  left_join(mean_temp_annual, by = c("GaugeID" = "GaugeID", "Year_Spawning" = "Year"), suffix = c("", "_spawning")) %>%
  left_join(comb_gauge1, by = c("GaugeID" = "GaugeID", "Year_2yo" = "Year"), suffix = c("", "_2yo"))

# Display the first few rows of the updated full_data
head(full_data)



# Rearranging, renaming and selecting variables from full_data for clarity
fishData <- full_data %>%
  dplyr::select(
    FishID_normalized,
    Year_2yo,
    Year_1yo,
    Year_Spawning,
    SiteID,
    GaugeID,
    River,
    Age,
    Annual_Growth_Rate,
    Species,
    mean_StreamLvl_m = mean_level,
    mean_StreamFlow_megalitre = mean_discharge,
    MeanTemp_degC_2yo = mean_temp,
    MeanTemp_degC_1yo = mean_temp_1yo,
    MeanTemp_degC_spawning = mean_temp_spawning,
    Catchment_area_sq_kms = `Catchment area_sq.kms`,
    Stream_Distance_from_station_to_mouth_km = `Stream Distance from station to mouth_km`,
    LongestZeroFlow_2yo = LongestZeroFlow,
    LongestZeroFlow_1yo,
    LongestZeroFlow_spawning,
    ZeroFlowDays_2yo = ZeroFlowDays,
    ZeroFlowDays_1yo,
    ZeroFlowDays_spawning,
    AvgEvap_mm_2yo = mean_evap_pan,
    AvgDailyRain_mm_2yo = mean_daily_rain
  )

head(fishData)

## Filtering data to only include fish in the correct age range (at least 2yo in 2021)
fishData <- fishData %>%
  filter(Year_1yo != 2021)

# Checking the first few rows of the filtered dataset
head(fishData)



#########################################################################
## Basic Data Visualisation Prior to Modeling
#########################################################################

#########################################################################
## Creating a table of data on sites for report
unique_sites <- full_data %>%
  select(SiteID, GaugeID, River, Control, Latitude, Longitude) %>%
  distinct(SiteID, .keep_all = TRUE) %>%
  arrange(SiteID)  # Sorting by SiteID in ascending order

# saving table as a CSV
write.csv(unique_sites, "Unique_Site_Information.csv", row.names = FALSE)


#########################################################################
## Visualising correlation matrix for pertinent variables

# Computing the correlation matrix for a subset of variables (not including species column)
corr_matrix <- fishData %>%
  dplyr::select(9, 11:ncol(.)) %>%
  cor(use = "complete.obs")

# Creating a mixed correlation plot

corrplot.mixed(
  corr_matrix,
  tl.cex = 0.5,  # Size of text labels
  number.cex = 0.7,  # Size of correlation coefficients
  upper = "ellipse",
  color = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
  tl.col = "black"
)


# Tabulating the corrmatrix for ease of reference
corr_data <- as.data.frame(corr_matrix) %>%
  rownames_to_column(var = "Variable1") %>%
  pivot_longer(cols = -Variable1, names_to = "Variable2", values_to = "Correlation") %>%
  filter(Variable1 != Variable2) 

#########################################################################
## Scatterplot to look at spread of data and 'relationships'

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(GGally)

# Assigning colours to species
species_colors <- c("BonyBream" = "blue", "GoldenPerch" = "salmon", "CommonCarp" = "darkgrey")


# Ensuring the dataset 'fishData' contains the 'Species' column used for coloring
plot_data <- fishData %>%
  dplyr::select(9:ncol(.))

# Check names and adjust formula as needed based on actual names in 'plot_data'
names(plot_data)

# Removing NA records
plot_data <- plot_data %>%
  na.omit()
plot_data <- plot_data %>%
  filter(!is.na(Species))

# Creating the scatter plot matrix
pairs(~ Annual_Growth_Rate + mean_StreamLvl_m + mean_StreamFlow_megalitre +
        MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning +
        Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km +
        LongestZeroFlow_2yo + LongestZeroFlow_1yo + LongestZeroFlow_spawning +
        ZeroFlowDays_2yo + ZeroFlowDays_1yo + ZeroFlowDays_spawning + 
        AvgEvap_mm_2yo + AvgDailyRain_mm_2yo,
      data = plot_data,
      col = species_colors[plot_data$Species],
      pch = 16)

# Scatterplot matrix split into different sections for ease of visual examination
pairs(~ Annual_Growth_Rate + LongestZeroFlow_2yo + LongestZeroFlow_1yo + LongestZeroFlow_spawning +
        ZeroFlowDays_2yo + ZeroFlowDays_1yo + ZeroFlowDays_spawning,
      data = plot_data,
      col = species_colors[plot_data$Species],
      pch = 16,
      main = "Scatter Plot Matrix for Stream Level and Flow Variables")

pairs(~ Annual_Growth_Rate + MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning,
      data = plot_data,
      col = species_colors[plot_data$Species],
      pch = 16,
      main = "Scatter Plot Matrix for Temperature and Zero Flow Variables")

pairs(~ Annual_Growth_Rate + ZeroFlowDays_2019 + ZeroFlowDays_2020 + ZeroFlowDays_2021 +
        AvgEvap_pan + AvgDailyRain,
      data = plot_data,
      col = species_colors[plot_data$Species],
      pch = 16,
      main = "Scatter Plot Matrix for Zero Flow Days, Evaporation, and Rain Variables")


#########################################################################
## Scatter plot with trend lines

# Loading the necessary libraries
library(ggplot2)
library(gridExtra) # For arranging multiple plots

# Creating scatter plots with smoothing lines for temperature variables
p1 <- ggplot(plot_data, aes(x = MeanTemp_degC_2yo, y = Annual_Growth_Rate, color = Species)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = species_colors) +
  theme_classic() +
  labs(x = "Mean Temp. 2yo (°C)", y = "Annual Growth Rate") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") # Hide legend

p2 <- ggplot(plot_data, aes(x = MeanTemp_degC_1yo, y = Annual_Growth_Rate, color = Species)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = species_colors) +
  theme_classic() +
  labs(x = "Mean Temp. 1yo (°C)", y = NULL) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") # Hide legend

p3 <- ggplot(plot_data, aes(x = MeanTemp_degC_spawning, y = Annual_Growth_Rate, color = Species)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = species_colors) +
  theme_classic() +
  labs(x = "Mean Temp. Spawning (°C)", y = NULL) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") # Hide legend


# Arranging the temp plots in a grid
grid.arrange(p1, p2, p3, nrow = 1)


# Create scatter plots with smoothing lines for longest zero flow variables
p1 <- ggplot(plot_data, aes(x = LongestZeroFlow_2yo, y = Annual_Growth_Rate, color = Species)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = species_colors) +
  theme_classic() +
  labs(x = "Longest Zero Flow (2yo)", y = "Annual Growth Rate") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") # Hide legend

p2 <- ggplot(plot_data, aes(x = LongestZeroFlow_1yo, y = Annual_Growth_Rate, color = Species)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = species_colors) +
  theme_classic() +
  labs(x = "Longest Zero Flow (1yo)", y = NULL) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") # Hide legend

p3 <- ggplot(plot_data, aes(x = LongestZeroFlow_spawning, y = Annual_Growth_Rate, color = Species)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = species_colors) +
  theme_classic() +
  labs(x = "Longest Zero Flow (Spawning)", y = NULL) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") # Hide legend

# Arranging the zero flow plots in a grid
grid.arrange(p1, p2, p3, nrow = 1)

#########################################################################
## Exploring the structure of independent variables

## Creating histograms and box plots for continuous predictor variables

# List of continuous variables to plot
continuous_vars <- c("Annual_Growth_Rate", "mean_StreamLvl_m", "mean_StreamFlow_megalitre",
                     "MeanTemp_degC_2yo", "MeanTemp_degC_1yo", "MeanTemp_degC_spawning",
                     "Catchment_area_sq_kms", "Stream_Distance_from_station_to_mouth_km",
                     "LongestZeroFlow_2yo", "LongestZeroFlow_1yo", "LongestZeroFlow_spawning",
                     "ZeroFlowDays_2yo", "ZeroFlowDays_1yo", "ZeroFlowDays_spawning", "AvgEvap_mm_2yo", "AvgDailyRain_mm_2yo")

# Plot histograms
for (var in continuous_vars) {
  p <- ggplot(fishData, aes_string(x = var)) +
    geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = paste("Distribution of", var), x = var, y = "Frequency") +
    theme_minimal()
  print(p)
}


# boxplots

for (var in continuous_vars) {
  p <- ggplot(fishData, aes_string(x = "Species", y = var, fill = "Species")) +
    geom_boxplot() +
    scale_fill_manual(values = species_colors) +
    stat_summary(fun.y = median, geom = "text", aes(label = round(..y.., 2)), 
                 position = position_dodge(width = 0.75), vjust = -0.5, color = "black") +
    labs(title = paste("Boxplot of", var, "by Species"), x = "Species", y = var) +
    theme_minimal()
  print(p)
}

#########################################################################
## Bar graph for Species

# Mapping for species names
species_names <- c("CommonCarp" = "Common Carp", "GoldenPerch" = "Golden Perch", "BonyBream" = "Bony Bream")

# Creating bar plot for Species
p <- ggplot(fishData, aes(x = Species, fill = Species)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  scale_fill_manual(values = species_colors, labels = species_names) +
  scale_x_discrete(labels = species_names) +
  labs(x = "Species", y = "Count") +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank())

print(p)

#Visualising the distribution of zeros across sites and species
all_site_ids <- sort(unique(fishData$SiteID))

ggplot(fishData %>% filter(Annual_Growth_Rate == 0), aes(x = factor(SiteID), fill = Species)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(breaks = all_site_ids) +
  labs(title = "", 
       x = "Site ID", 
       y = "Count of Zero Growth Rates") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##filtering out zero counts in Annual growth rate, based on this result
fishData_Zrem <- fishData %>% filter(Annual_Growth_Rate > 0)


# Plot temperature vs. longest zero flow to visually inspect relationships
ggplot(fishData, aes(x = MeanTemp_degC_2yo, y = LongestZeroFlow_2yo)) +
  geom_point() +
  labs(title = "Mean Temperature vs. Longest Zero Flow (Year 2yo)", x = "Mean Temperature (2yo)", y = "Longest Zero Flow (2yo)")

#########################################################################
## Model Fitting
#########################################################################


#########################################################################
## Model 1: Basic model
m1_lm_all <- lm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre +
                    MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning +
                    Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km +
                    LongestZeroFlow_2yo + LongestZeroFlow_1yo + LongestZeroFlow_spawning +
                    ZeroFlowDays_2yo + ZeroFlowDays_1yo + ZeroFlowDays_spawning + 
                    AvgEvap_mm_2yo + AvgDailyRain_mm_2yo + Species, data = fishData)

# This initial base model uses the full suite of predictor variables
summary(m1_lm_all)
# Testing for Homoscedasticity and Normality
plot(simulateResiduals(m1_lm_all))
# testing for outliers
testOutliers(simulateResiduals(m1_lm_all), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m1_lm_all))
# Extracting model equation
extract_eq(m1_lm_all)


## Removing highly correlated and irrelevant predictors from base model

m1_lm_s <- lm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre +
                MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning +
                LongestZeroFlow_2yo + Species, data = fishData)

# This initial base model uses the full suite of predictor variables
summary(m1_lm_s)
# Testing for Homoscedasticity and Normality
plot(simulateResiduals(m1_lm_s))
# testing for outliers
testOutliers(simulateResiduals(m1_lm_s), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m1_lm_s))
# Extracting model equation
extract_eq(m1_lm_s)


#########################################################################
# Model 2: log transforming variables that are skewed 

# Check skewness of variables
# Calculate skewness excluding NA values
skewness_data_clean <- sapply(fishData[, c("Annual_Growth_Rate", "mean_StreamLvl_m", "mean_StreamFlow_megalitre",
                                           "MeanTemp_degC_2yo", "MeanTemp_degC_1yo", "MeanTemp_degC_spawning",
                                           "LongestZeroFlow_2yo")], 
                              function(x) skewness(na.omit(x)))

print(skewness_data_clean)

# Updtaing fishData to include log-transformations
fishData <- fishData %>%
  mutate(log_mean_StreamFlow_megalitre = log1p(mean_StreamFlow_megalitre),
         log_mean_StreamLvl_m = log1p(mean_StreamLvl_m))


# Fitting the model with log-transformed mean_StreamFlow_megalitre
m2_transformed <- lm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                       MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                       LongestZeroFlow_2yo + Species, data = fishData)

summary(m2_transformed)

# Diagnostics
# Testing for Homoscedasticity and Normality
plot(simulateResiduals(m2_transformed))
# testing for outliers
testOutliers(simulateResiduals(m2_transformed), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m2_transformed))
# esxtracting equation
extract_eq(m2_transformed)
#No significant improvement from the un-transformed model


#########################################################################
# Model 3: Introducing interaction terms
# Because of limited success of m1/m2 in explaining variance (low R-squared)

m3_interaction <- lm(Annual_Growth_Rate ~ Species * (log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                                                       MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                                                       LongestZeroFlow_2yo), data = fishData)
summary(m3_interaction)

# Diagnostics
# Testing for Homoscedasticity and Normality
plot(simulateResiduals(m3_interaction))
# testing for outliers
testOutliers(simulateResiduals(m3_interaction), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m3_interaction))
# Goodness of fit
pregibon.lm(m3_interaction)
# Extracting equation
extract_eq(m3_interaction)


#########################################################################
## Model 4: Species specific lms
# Bony Bream Model
m4_lm_bb <- lm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                     MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                     LongestZeroFlow_2yo, data = fishData, subset = Species == "BonyBream")

# Common Carp Model
m4_lm_cc <- lm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                      MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                      LongestZeroFlow_2yo, 
                    data = fishData, subset = Species == "GoldenPerch")

# Golden Perch Model
m4_lm_gp <- lm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                 MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                 LongestZeroFlow_2yo, 
               data = fishData, subset = Species == "CommonCarp")

# Summarising the models
summary(m4_lm_bb)
summary(m4_lm_cc)
summary(m4_lm_gp)

# Diagnostics
# Testing for Homoscedasticity and Normality
plot(simulateResiduals(m4_lm_bb))
plot(simulateResiduals(m4_lm_cc))
plot(simulateResiduals(m4_lm_gp))

# Assumptions of Homoscedasticity and Normality have not been met in lms


#########################################################################
# Model 5: General Linear Model with standard identity link

m5_glm_all <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                     MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                     LongestZeroFlow_2yo + Species, 
              family = gaussian(), data = fishData)

summary(m5_glm_all)

#Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m5_glm_all))
## Goodness of fit testing via pregibon
pregibon.glm(m5_glm_all)
# extract equation
extract_eq(m5_glm_all)


## Base GLM without zero counts for annual growth rate
# Repeating code used for filtering process above (for reference)
fishData_Zrem <- fishData %>% filter(Annual_Growth_Rate > 0)

# Model:
m5_glm_Zrem <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                    MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                    LongestZeroFlow_2yo + Species, 
                  family = gaussian(), data = fishData_Zrem)

summary(m5_glm_Zrem)

#Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m5_glm_Zrem))
## Goodness of fit testing via pregibon
pregibon.glm(m5_glm_all)
# extract equation
extract_eq(m5_glm_Zrem)


# Species specific GLMs: Bony Bream

m5_glm_bb <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                   MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                   LongestZeroFlow_2yo, 
                   family = gaussian(), data = fishData_Zrem, subset = Species == "BonyBream")

summary(m5_glm_bb)

#Diagnostics
# Testing for Homoscedasticity and Normality
plot(simulateResiduals(m5_glm_bb))
## Goodness of fit testing via pregibon
pregibon.glm(m5_glm_bb)


# Species specific GLMs: Common Carp

m5_glm_cc <- glm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre + 
                   MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                   LongestZeroFlow_2yo, 
                 family = gaussian(), data = fishData_Zrem, subset = Species == "CommonCarp")

# Summarize the model
summary(m5_glm_cc)

# Diagnostics
plot(simulateResiduals(m5_glm_cc))
pregibon.glm(m5_glm_cc)


# Fit GLM for Golden Perch
m5_glm_gp <- glm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre + 
                   MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                   LongestZeroFlow_2yo, 
                 family = gaussian(), data = fishData_Zrem, subset = Species == "GoldenPerch")

# Summarize the model
summary(m5_glm_gp)

# Diagnostics
plot(simulateResiduals(m5_glm_gp))
pregibon.glm(m5_glm_gp)


#########################################################################
# Model 6: General Linear Model with interaction terms

m6_glm_int <- glm(Annual_Growth_Rate ~ Species * (log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                                                    MeanTemp_degC_2yo + MeanTemp_degC_1yo + MeanTemp_degC_spawning + 
                                                    LongestZeroFlow_2yo), 
                  family = gaussian(), data = fishData)

summary(m6_glm_int)

#Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m6_glm_int))
## Goodness of fit testing via pregibon
pregibon.glm(m6_glm_int)


#####################################################################
## Model 7: GLM with polynomial terms

# Updated model with polynomial terms for temperature and longest zero flow variables
# Ensure no NA values
fishData_NArem <- na.omit(fishData)
fishData_clean <- na.omit(fishData_Zrem)

m7_glm_poly_full <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                          poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                          poly(LongestZeroFlow_2yo, 2) + 
                          Species, 
                        family = gaussian(), data = fishData_NArem)
summary(m7_glm_poly_full)


## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_full))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_full)
# Extracting Equation
extract_eq(m7_glm_poly_full)



## Model with both NA and 0 removed (did not pass Pregibon link test)
# so even though it is a better fit, it does not meet model specs
m7_glm_poly_clean <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                           poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                           poly(LongestZeroFlow_2yo, 2) + 
                           Species, 
                         family = gaussian(), data = fishData_clean)

summary(m7_glm_poly_clean)

## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_clean))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_clean)
# Extracting equation
extract_eq(m7_glm_poly_clean)


###############################################################
## Species-specific GLMs with polynomial terms

# Polynomial term GLMs for species: Bony Bream
m7_glm_poly_bb <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                        poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                        poly(LongestZeroFlow_2yo, 2), 
                        family = gaussian(), data = fishData_clean, subset = Species == "BonyBream")
summary(m7_glm_poly_bb)

## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_bb))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_bb)
# Extracting equation
extract_eq(m7_glm_poly_bb)


# Same model without zeros removed (does not pass pregibon link test)

m7_glm_poly_bb2 <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                        poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                        poly(LongestZeroFlow_2yo, 2), 
                      family = gaussian(), data = fishData_NArem, subset = Species == "BonyBream")

## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_bb2))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_bb2)
# Extracting equation
extract_eq(m7_glm_poly_bb2)


# Polynomial term GLMs for species: Common Carp (does not pass)
m7_glm_poly_cc <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                        poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                        poly(LongestZeroFlow_2yo, 2), 
                      family = gaussian(), data = fishData_clean, subset = Species == "CommonCarp")

## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_cc))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_cc)
# Extracting equation
extract_eq(m7_glm_poly_cc)

# Same model without zeros removed (passes)

m7_glm_poly_cc2 <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                         poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                         poly(LongestZeroFlow_2yo, 2), 
                       family = gaussian(), data = fishData_NArem, subset = Species == "CommonCarp")
summary(m7_glm_poly_cc2)


## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_cc2))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_cc2)
# Extracting equation
extract_eq(m7_glm_poly_cc2)


# Polynomial term GLMs for species: Golden Perch (passes pregibon link test)
m7_glm_poly_gp <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                        poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                        poly(LongestZeroFlow_2yo, 2), 
                      family = gaussian(), data = fishData_clean, subset = Species == "GoldenPerch")
summary(m7_glm_poly_gp)

## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_gp))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_gp)
# Extracting equation
extract_eq(m7_glm_poly_gp)

# Same model without zeros removed (passes pregibon link test)

m7_glm_poly_gp2 <- glm(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                         poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                         poly(LongestZeroFlow_2yo, 2), 
                       family = gaussian(), data = fishData_NArem, subset = Species == "GoldenPerch")

summary(m7_glm_poly_gp2)

## Diagnostics
# Testing for Homoscedasticity
plot(simulateResiduals(m7_glm_poly_gp2))
# Goodness of fit testing via pregibon
pregibon.glm(m7_glm_poly_gp2)
# Extracting equation
extract_eq(m7_glm_poly_gp2)



#########################################################################
## Model 8: Generalized Linear Mixed-Effects Models

library(lme4)

# GLMM with random intercepts for SiteID
m8_glmm <- glmer(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                  poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                  poly(LongestZeroFlow_2yo, 2) + Species + (1 | SiteID), data = fishData_NArem, family = gaussian())

summary(m8_glmm)

# Diagnostics
plot(simulateResiduals(m8_glmm))

# Same model with NA and 0 removed
m8_glmm2 <- glmer(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                   poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                   poly(LongestZeroFlow_2yo, 2) + Species + (1 | SiteID), data = fishData_clean, family = gaussian())

summary(m8_glmm2)

# Diagnostics
plot(simulateResiduals(m8_glmm2))




# Species specific glmms

# Bony Bream
m8_glmm_bb <- glmer(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                      poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                      poly(LongestZeroFlow_2yo, 2) + (1 | SiteID), 
                    data = fishData_NArem, subset = Species == "BonyBream", family = gaussian())
summary(m8_glmm_bb)

## Diagnostics for Bony Bream
plot(simulateResiduals(m8_glmm_bb))

# Common Carp
m8_glmm_cc <- glmer(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                      poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                      poly(LongestZeroFlow_2yo, 2) + (1 | SiteID), 
                    data = fishData_NArem, subset = Species == "CommonCarp", family = gaussian())
summary(m8_glmm_cc)

## Diagnostics for Common Carp
plot(simulateResiduals(m8_glmm_cc))

# Golden Perch
m8_glmm_gp <- glmer(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                      poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + 
                      poly(LongestZeroFlow_2yo, 2) + (1 | SiteID), 
                    data = fishData_NArem, subset = Species == "GoldenPerch", family = gaussian())
summary(m8_glmm_gp)

## Diagnostics for Golden Perch
plot(simulateResiduals(m8_glmm_gp))


#########################################################################
## Model 9: General Additive Models

# Fitting GAMs with reduced complexity for smooth terms
m9_gam_bb <- gam(Annual_Growth_Rate ~ s(log_mean_StreamLvl_m, k = 5) + s(log_mean_StreamFlow_megalitre, k = 5) + 
                   s(MeanTemp_degC_2yo, k = 5) + s(MeanTemp_degC_1yo, k = 5) + s(MeanTemp_degC_spawning, k = 5) + 
                   s(LongestZeroFlow_2yo, k = 5), data = fishData, subset = Species == "BonyBream")

# Summarize the new model
summary(m9_gam_bb)
plot(simulateResiduals(m9_gam_bb))

##################################################################
## Model 10: PCA

# Selecting the full suite of predictors
predictors <- fishData[, c("mean_StreamLvl_m", 
                           "mean_StreamFlow_megalitre", "MeanTemp_degC_2yo", 
                           "MeanTemp_degC_1yo", "MeanTemp_degC_spawning", 
                           "Catchment_area_sq_kms", 
                           "Stream_Distance_from_station_to_mouth_km", 
                           "LongestZeroFlow_2yo", "LongestZeroFlow_1yo", 
                           "LongestZeroFlow_spawning", "ZeroFlowDays_2yo", 
                           "ZeroFlowDays_1yo", "ZeroFlowDays_spawning", 
                           "AvgEvap_mm_2yo", "AvgDailyRain_mm_2yo")]

# Identifying columns with zero variance
zero_var_cols <- sapply(predictors, function(x) var(x, na.rm = TRUE) == 0)
print(zero_var_cols)

# Removing zero variance columns
predictors_clean <- predictors[, !zero_var_cols]

# Checking for infinite values
infinite_cols <- sapply(predictors_clean, function(x) any(is.infinite(x)))
print(infinite_cols)

# Removing columns with infinite values
predictors_clean <- predictors_clean[, !infinite_cols]

# Checking for NA values
na_cols <- sapply(predictors_clean, function(x) any(is.na(x)))
print(na_cols)

# Remove rows with NA values in the predictors from the entire dataset
fishData_clean2 <- fishData[complete.cases(predictors_clean), ]

# Selecting the cleaned predictors
predictors_clean2 <- fishData_clean2[, names(predictors_clean)]

# Applying PCA to the cleaned predictors
pca <- prcomp(predictors_clean2, scale. = TRUE)

# Check the proportion of variance explained by each principal component
summary(pca)

# Use the first few principal components that explain most of the variance
# Here, we assume the first three principal components are sufficient
fishData_pca <- data.frame(pca$x[, 1:3]) 
names(fishData_pca) <- c("PC1", "PC2", "PC3")
fishData_pca$Annual_Growth_Rate <- fishData_clean2$Annual_Growth_Rate
fishData_pca$Species <- fishData_clean2$Species

# Fitting the model with principal components
m10_pca <- glm(Annual_Growth_Rate ~ PC1 + PC2 + PC3 + Species, data = fishData_pca)
summary(m10_pca)

# Diagnostics
simulate_res <- simulateResiduals(m10_pca)
plot(simulate_res)
pregibon.glm(m10_pca)

# Scree plot
screeplot(pca, type = "lines", main = "Scree Plot")


# Biplot
# Extract the loadings (eigenvectors)
loadings <- as.data.frame(pca$rotation)
loadings$Variable <- rownames(loadings)
loadings$PC1 <- loadings$PC1 * 5  # Scale the loadings for better visualization
loadings$PC2 <- loadings$PC2 * 5
# Extracting the scores (principal component scores)
scores <- as.data.frame(pca$x)
scores$Observation <- rownames(scores)
# Creating biplot
ggplot() +
  geom_point(data = scores, aes(x = PC1, y = PC2), color = 'blue', alpha = 0.6) +
  geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = Observation), size = 3, color = 'blue', max.overlaps = 10) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = loadings, aes(x = PC1, y = PC2, label = Variable), size = 3, color = "red", max.overlaps = 10) +
  ggtitle("Biplot of Principal Components") +
  theme_minimal() +
  theme(text = element_text(size = 10))


# Loadings plot
library(ggplot2)
loadings <- as.data.frame(pca$rotation)
loadings$Variable <- rownames(loadings)

# Plot for the first two principal components
ggplot(loadings, aes(x = PC1, y = PC2, label = Variable)) +
  geom_point() +
  geom_text(vjust = -0.5, hjust = 0.5) +
  ggtitle("Loadings Plot") +
  xlab("PC1") +
  ylab("PC2")

#########################################################################
# Model 11: Bayesian GLM

#Checking where g++ and make are located
Sys.which("g++")
Sys.which("make")

# Filter data for each species
bony_bream_data <- fishData_clean %>% filter(Species == "BonyBream")
common_carp_data <- fishData_clean %>% filter(Species == "CommonCarp")
golden_perch_data <- fishData_clean %>% filter(Species == "GoldenPerch")

# Model for bony bream

# Define the formula for Bony Bream
formula_bb <- bf(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                   poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + 
                   poly(MeanTemp_degC_spawning, 2) + poly(LongestZeroFlow_2yo, 2))

# Fit the Bayesian model for Bony Bream
m11_bayesian_bb <- brm(formula = formula_bb, 
                data = bony_bream_data, 
                family = gaussian(), 
                chains = 4, 
                iter = 2000, 
                control = list(adapt_delta = 0.95))

# Summarize the model
summary(m11_bayesian_bb)



## Common Carp

# Define the formula for Common Carp
formula_cc <- bf(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                   poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + 
                   poly(MeanTemp_degC_spawning, 2) + poly(LongestZeroFlow_2yo, 2))

# Fit the Bayesian model for Common Carp
m11_bayesian_cc <- brm(formula = formula_cc, 
                data = common_carp_data, 
                family = gaussian(), 
                chains = 4, 
                iter = 2000, 
                control = list(adapt_delta = 0.95))

# Summarize the model
summary(m11_bayesian_cc)



# Golden Perch

# Define the formula for Golden Perch
formula_gp <- bf(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                   poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + 
                   poly(MeanTemp_degC_spawning, 2) + poly(LongestZeroFlow_2yo, 2))

# Fit the Bayesian model for Golden Perch
m11_bayesian_gp <- brm(formula = formula_gp, 
                data = golden_perch_data, 
                family = gaussian(), 
                chains = 4, 
                iter = 2000, 
                control = list(adapt_delta = 0.95))

# Summarize the model
summary(m11_bayesian_gp)



#########################################################################
# Model Comparison 

# Load necessary libraries
library(caret)
library(brms)

# Function for cross-validation
k_fold_cv <- function(model, data, k = 10) {
  folds <- createFolds(data$Annual_Growth_Rate, k = k, list = TRUE, returnTrain = FALSE)
  cv_results <- lapply(folds, function(test_indices) {
    train_indices <- setdiff(1:nrow(data), test_indices)
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    # Fit model
    fit <- update(model, data = train_data)
    
    # Predict and compute error metrics
    predictions <- predict(fit, newdata = test_data)
    actuals <- test_data$Annual_Growth_Rate
    
    # Compute RMSE and MAE
    rmse <- sqrt(mean((predictions - actuals)^2))
    mae <- mean(abs(predictions - actuals))
    
    return(c(RMSE = rmse, MAE = mae))
  })
  
  # Combine results
  cv_results <- do.call(rbind, cv_results)
  return(cv_results)
}

# Bayesian model using brms
m7_bayesian_cv <- function(train_data, test_data) {
  fit <- brm(
    Annual_Growth_Rate ~ PC1 + PC2 + PC3 + Species,
    data = train_data,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "b"), 
      prior(student_t(3, 0, 10), class = "Intercept")
    ),
    iter = 2000, chains = 2, cores = 2,
    control = list(adapt_delta = 0.95)
  )
  predictions <- as.numeric(predict(fit, newdata = test_data))
  return(predictions)
}

# GLM models
m6_pca_cv <- k_fold_cv(m6_pca, fishData_pca, k = 10)
m4_glm_species_cv <- k_fold_cv(m4_glm_species, fishData, k = 10)

# Define custom cross-validation function for Bayesian model
k_fold_cv_bayesian <- function(model_func, data, k = 10) {
  folds <- createFolds(data$Annual_Growth_Rate, k = k, list = TRUE, returnTrain = FALSE)
  cv_results <- lapply(folds, function(test_indices) {
    train_indices <- setdiff(1:nrow(data), test_indices)
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    # Fit and predict using Bayesian model function
    predictions <- model_func(train_data, test_data)
    actuals <- test_data$Annual_Growth_Rate
    
    # Compute RMSE and MAE
    rmse <- sqrt(mean((predictions - actuals)^2))
    mae <- mean(abs(predictions - actuals))
    
    return(c(RMSE = rmse, MAE = mae))
  })
  
  # Combine results
  cv_results <- do.call(rbind, cv_results)
  return(cv_results)
}

# Run cross-validation for Bayesian model
m7_bayesian_cv_results <- k_fold_cv_bayesian(m7_bayesian_cv, fishData_pca, k = 10)

# Aggregate results and compare
glm_pca_results <- colMeans(m6_pca_cv)
glm_species_results <- colMeans(m4_glm_species_cv)
bayesian_results <- colMeans(m7_bayesian_cv_results)

comparison_results <- rbind(
  GLM_PCA = glm_pca_results,
  GLM_Species = glm_species_results,
  Bayesian = bayesian_results
)

print(comparison_results)



#########################################################################
#Visualising the Bayesian GLMs
#########################################################################
plot(m11_bayesian_bb)
plot(m11_bayesian_cc)
plot(m11_bayesian_gp)

pp_check(m11_bayesian_bb, type = "dens_overlay")
pp_check(m11_bayesian_cc, type = "dens_overlay")
pp_check(m11_bayesian_gp, type = "dens_overlay")



#########################################################################
#Visualising the GLMs
#########################################################################
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to get predictions and confidence intervals
get_predictions <- function(model, newdata) {
  preds <- predict(model, newdata = newdata, type = "response", se.fit = TRUE)
  newdata$fit <- preds$fit
  newdata$upr <- preds$fit + 1.96 * preds$se.fit
  newdata$lwr <- preds$fit - 1.96 * preds$se.fit
  return(newdata)
}

# Example datasets for each species
bony_bream_data <- fishData_clean %>% filter(Species == "BonyBream")
common_carp_data <- fishData_clean %>% filter(Species == "CommonCarp")
golden_perch_data <- fishData_clean %>% filter(Species == "GoldenPerch")

# Get predictions for GLMs
preds_glm_bb <- get_predictions(m7_glm_poly_bb, bony_bream_data)
preds_glm_cc <- get_predictions(m7_glm_poly_cc2, common_carp_data)
preds_glm_gp <- get_predictions(m7_glm_poly_gp, golden_perch_data)

# Combine all data for GLM plotting
combined_glm_data <- bind_rows(
  preds_glm_bb %>% mutate(Species = "BonyBream"),
  preds_glm_cc %>% mutate(Species = "CommonCarp"),
  preds_glm_gp %>% mutate(Species = "GoldenPerch")
)

# Summarize predictions
glm_summary <- combined_glm_data %>%
  group_by(Species) %>%
  summarise(
    fit = mean(fit),
    lwr = mean(lwr),
    upr = mean(upr)
  )

# Ensure Species is a factor with the correct levels
glm_summary$Species <- factor(glm_summary$Species, levels = c("BonyBream", "CommonCarp", "GoldenPerch"))

# Assigning colours to species
species_colors <- c("BonyBream" = "blue", "CommonCarp" = "darkgrey", "GoldenPerch" = "salmon")

# Plot GLM predictions
p_glm <- ggplot(glm_summary, aes(x = Species, y = fit, color = Species)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  scale_color_manual(values = species_colors) +
  labs(x = "Species", y = "Annual Growth Rate", title = "GLM Predicted Annual Growth Rate by Species with 95% CI") +
  theme_minimal()

print(p_glm)


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Generate prediction data
temp_seq <- seq(min(fishData_clean$MeanTemp_degC_2yo, na.rm = TRUE), max(fishData_clean$MeanTemp_degC_2yo, na.rm = TRUE), length.out = 100)
predict_data <- data.frame(MeanTemp_degC_2yo = temp_seq)

# Adding other predictors as constants (mean or typical value)
predict_data <- predict_data %>%
  mutate(log_mean_StreamLvl_m = mean(fishData_clean$log_mean_StreamLvl_m, na.rm = TRUE),
         log_mean_StreamFlow_megalitre = mean(fishData_clean$log_mean_StreamFlow_megalitre, na.rm = TRUE),
         MeanTemp_degC_1yo = mean(fishData_clean$MeanTemp_degC_1yo, na.rm = TRUE),
         MeanTemp_degC_spawning = mean(fishData_clean$MeanTemp_degC_spawning, na.rm = TRUE),
         LongestZeroFlow_2yo = mean(fishData_clean$LongestZeroFlow_2yo, na.rm = TRUE))

# Add polynomial terms
predict_data$polyMeanTemp_degC_2yo_1 <- predict_data$MeanTemp_degC_2yo
predict_data$polyMeanTemp_degC_2yo_2 <- I(predict_data$MeanTemp_degC_2yo^2)

# Model prediction
predictions <- predict(m7_glm_poly_bb, newdata = predict_data, se = TRUE)
predict_data$fit <- predictions$fit
predict_data$fit_up <- predict_data$fit + 1.96 * predictions$se.fit
predict_data$fit_low <- predict_data$fit - 1.96 * predictions$se.fit

# Plot
ggplot(predict_data, aes(x = MeanTemp_degC_2yo, y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2) +
  labs(x = "Mean Temperature at Year 2 (°C)", y = "Predicted Annual Growth Rate",
       title = "Effect of Mean Temperature at Year 2 on Annual Growth Rate") +
  theme_minimal()

#############################################################
library(dplyr)

# Assuming mean values for other predictors
predict_data_bb <- expand.grid(
  log_mean_StreamLvl_m = stream_lvl_seq,
  log_mean_StreamFlow_megalitre = stream_flow_seq
)

# Add constant predictors at their mean values
predict_data_bb$MeanTemp_degC_2yo <- mean(fishData_clean$MeanTemp_degC_2yo, na.rm = TRUE)
predict_data_bb$MeanTemp_degC_1yo <- mean(fishData_clean$MeanTemp_degC_1yo, na.rm = TRUE)
predict_data_bb$MeanTemp_degC_spawning <- mean(fishData_clean$MeanTemp_degC_spawning, na.rm = TRUE)
predict_data_bb$LongestZeroFlow_2yo <- mean(fishData_clean$LongestZeroFlow_2yo, na.rm = TRUE)

# Create a vector with sufficient range for polynomial calculation
temp_range_2yo <- seq(min(fishData_clean$MeanTemp_degC_2yo, na.rm = TRUE), max(fishData_clean$MeanTemp_degC_2yo, na.rm = TRUE), length.out = 100)

# Pre-calculate polynomial terms from a range of values
poly_temp_2yo <- poly(temp_range_2yo, 2)  # Creates matrix of polynomial terms

# Assign these polynomial terms back to your prediction data
predict_data_bb$polyMeanTemp_degC_2yo_1 <- rep(poly_temp_2yo[,1], each = nrow(predict_data_bb)/length(temp_range_2yo))
predict_data_bb$polyMeanTemp_degC_2yo_2 <- rep(poly_temp_2yo[,2], each = nrow(predict_data_bb)/length(temp_range_2yo))

# Repeat for other temperature variables as needed using their specific ranges if polynomial terms are required

# Now calculate predictions
predictions_bb <- predict(m7_glm_poly_bb, newdata = predict_data_bb, type = "response", se = TRUE)
predict_data_bb$fit <- predictions_bb$fit
predict_data_bb$fit_up <- predict_data_bb$fit + 1.96 * predictions_bb$se.fit
predict_data_bb$fit_low <- predict_data_bb$fit - 1.96 * predictions_bb$se.fit

# Plot for log_mean_StreamLvl_m
ggplot(predict_data_bb, aes(x = log_mean_StreamLvl_m, y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2, fill = "blue") +
  labs(x = "Log Mean Stream Level (m)", y = "Predicted Annual Growth Rate",
       title = "Effect of Log Mean Stream Level on Annual Growth Rate in Bony Bream") +
  theme_minimal()

# Plot for log_mean_StreamFlow_megalitre
ggplot(predict_data_bb, aes(x = log_mean_StreamFlow_megalitre, y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2, fill = "red") +
  labs(x = "Log Mean Stream Flow (Megalitre)", y = "Predicted Annual Growth Rate",
       title = "Effect of Log Mean Stream Flow on Annual Growth Rate in Bony Bream") +
  theme_minimal()

# Generate a sequence of values for the predictor
stream_flow_seq <- seq(min(fishData_clean$log_mean_StreamFlow_megalitre, na.rm = TRUE), 
                       max(fishData_clean$log_mean_StreamFlow_megalitre, na.rm = TRUE), 
                       length.out = 100)

# Create data frame for prediction
predict_data_bb <- data.frame(log_mean_StreamFlow_megalitre = stream_flow_seq)

# Adding other predictors as constants (mean values)
predict_data_bb$log_mean_StreamLvl_m <- mean(fishData_clean$log_mean_StreamLvl_m, na.rm = TRUE)
predict_data_bb$MeanTemp_degC_2yo <- mean(fishData_clean$MeanTemp_degC_2yo, na.rm = TRUE)
predict_data_bb$MeanTemp_degC_1yo <- mean(fishData_clean$MeanTemp_degC_1yo, na.rm = TRUE)
predict_data_bb$MeanTemp_degC_spawning <- mean(fishData_clean$MeanTemp_degC_spawning, na.rm = TRUE)
predict_data_bb$LongestZeroFlow_2yo <- mean(fishData_clean$LongestZeroFlow_2yo, na.rm = TRUE)

# Calculate polynomial terms if used in the model
predict_data_bb$polyMeanTemp_degC_2yo_1 <- poly(predict_data_bb$MeanTemp_degC_2yo, 2)[,1]  # First degree term of poly for 2nd year temperature
predict_data_bb$polyMeanTemp_degC_2yo_2 <- poly(predict_data_bb$MeanTemp_degC_2yo, 2)[,2]  # Second degree term
predict_data_bb$polyMeanTemp_degC_1yo_1 <- poly(predict_data_bb$MeanTemp_degC_1yo, 2)[,1]
predict_data_bb$polyMeanTemp_degC_1yo_2 <- poly(predict_data_bb$MeanTemp_degC_1yo, 2)[,2]
predict_data_bb$polyMeanTemp_degC_spawning_1 <- poly(predict_data_bb$MeanTemp_degC_spawning, 2)[,1]
predict_data_bb$polyMeanTemp_degC_spawning_2 <- poly(predict_data_bb$MeanTemp_degC_spawning, 2)[,2]
predict_data_bb$polyLongestZeroFlow_2yo_1 <- poly(predict_data_bb$LongestZeroFlow_2yo, 2)[,1]
predict_data_bb$polyLongestZeroFlow_2yo_2 <- poly(predict_data_bb$LongestZeroFlow_2yo, 2)[,2]

# Calculate predictions
predictions_bb <- predict(m7_glm_poly_bb, newdata = predict_data_bb, type = "response", se = TRUE)
predict_data_bb$fit <- predictions_bb$fit
predict_data_bb$fit_up <- predict_data_bb$fit + 1.96 * predictions_bb$se.fit
predict_data_bb$fit_low <- predict_data_bb$fit - 1.96 * predictions_bb$se.fit

# Plot the predictions
ggplot(predict_data_bb, aes(x = log_mean_StreamFlow_megalitre, y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2) +
  labs(x = "Log Mean Stream Flow (Megalitre)", y = "Predicted Annual Growth Rate",
       title = "Effect of Log Mean Stream Flow on Annual Growth Rate in Bony Bream") +
  theme_minimal()




#################################################################
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to generate prediction data and plot for a specific variable
generate_plot <- function(variable, model, data, xlabel, title) {
  var_seq <- seq(min(data[[variable]], na.rm = TRUE), max(data[[variable]], na.rm = TRUE), length.out = 100)
  predict_data <- data.frame(var_seq)
  names(predict_data) <- variable
  
  # Adding other predictors as constants (mean or typical value)
  other_vars <- setdiff(names(data), variable)
  for (var in other_vars) {
    predict_data[[var]] <- mean(data[[var]], na.rm = TRUE)
  }
  
  # Add polynomial terms if applicable
  poly_terms <- names(coef(model))[grepl(paste0("poly\\(", variable), names(coef(model)))]
  for (term in poly_terms) {
    degree <- as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", term))
    predict_data[[term]] <- I(predict_data[[variable]]^degree)
  }
  
  # Model prediction
  predictions <- predict(model, newdata = predict_data, se = TRUE)
  predict_data$fit <- predictions$fit
  predict_data$fit_up <- predict_data$fit + 1.96 * predictions$se.fit
  predict_data$fit_low <- predict_data$fit - 1.96 * predictions$se.fit
  
  # Plot
  ggplot(predict_data, aes_string(x = variable, y = "fit")) +
    geom_line() +
    geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2) +
    labs(x = xlabel, y = "Predicted Annual Growth Rate", title = title) +
    theme_minimal()
}


# Generate and save plots
p1 <- generate_plot("log_mean_StreamLvl_m", m7_glm_poly_bb, fishData_clean, "Log Mean Stream Level (m)", "Effect of Log Mean Stream Level on Annual Growth Rate")
p2 <- generate_plot("log_mean_StreamFlow_megalitre", m7_glm_poly_bb, fishData_clean, "Log Mean Stream Flow (megalitre)", "Effect of Log Mean Stream Flow on Annual Growth Rate")

# Display plots
print(p1)
print(p2)




# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to generate plot
generate_plot <- function(variable, model, data, x_label, plot_title) {
  temp_seq <- seq(min(data[[variable]], na.rm = TRUE), max(data[[variable]], na.rm = TRUE), length.out = 100)
  predict_data <- data.frame(temp_seq)
  colnames(predict_data) <- variable
  
  # Adding other predictors as constants (mean or typical value)
  predict_data <- predict_data %>%
    mutate(log_mean_StreamLvl_m = mean(data$log_mean_StreamLvl_m, na.rm = TRUE),
           log_mean_StreamFlow_megalitre = mean(data$log_mean_StreamFlow_megalitre, na.rm = TRUE),
           MeanTemp_degC_2yo = mean(data$MeanTemp_degC_2yo, na.rm = TRUE),
           MeanTemp_degC_1yo = mean(data$MeanTemp_degC_1yo, na.rm = TRUE),
           LongestZeroFlow_2yo = mean(data$LongestZeroFlow_2yo, na.rm = TRUE))
  
  # Apply polynomial terms
  predict_data[[paste0("poly(", variable, ", 2)1")]] <- poly(predict_data[[variable]], 2)[,1]
  predict_data[[paste0("poly(", variable, ", 2)2")]] <- poly(predict_data[[variable]], 2)[,2]
  
  # Model prediction
  predictions <- predict(model, newdata = predict_data, se = TRUE)
  predict_data$fit <- predictions$fit
  predict_data$fit_up <- predict_data$fit + 1.96 * predictions$se.fit
  predict_data$fit_low <- predict_data$fit - 1.96 * predictions$se.fit
  
  # Plot
  ggplot(predict_data, aes_string(x = variable, y = "fit")) +
    geom_line() +
    geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2) +
    labs(x = x_label, y = "Predicted Annual Growth Rate", title = plot_title) +
    theme_minimal()
}

# Plot for poly(MeanTemp_degC_spawning, 2)1
p1 <- generate_plot("MeanTemp_degC_spawning", m7_glm_poly_cc2, fishData_NArem, "Polynomial Term for Mean Temperature at Spawning (°C)", "Effect of Polynomial Term for Mean Temperature at Spawning on Annual Growth Rate")
print(p1)

# Plot for log_mean_StreamLvl_m
p2 <- generate_plot("log_mean_StreamLvl_m", m7_glm_poly_cc2, fishData_NArem, "Log Mean Stream Level (m)", "Effect of Log Mean Stream Level on Annual Growth Rate")
print(p2)





# Generate plot for poly(MeanTemp_degC_spawning, 2)1
p1 <- generate_plot("MeanTemp_degC_spawning", m7_glm_poly_gp, fishData_clean, "Polynomial Term for Mean Temperature at Spawning (°C)", "Effect of Polynomial Term for Mean Temperature at Spawning on Annual Growth Rate")
print(p1)

# Generate plot for log_mean_StreamFlow_megalitre
p2 <- generate_plot("log_mean_StreamFlow_megalitre", m7_glm_poly_gp, fishData_clean, "Log Mean Stream Flow (megalitre)", "Effect of Log Mean Stream Flow on Annual Growth Rate")
print(p2)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to generate plot
generate_plot <- function(variable, model, data, x_label, plot_title) {
  temp_seq <- seq(min(data[[variable]], na.rm = TRUE), max(data[[variable]], na.rm = TRUE), length.out = 100)
  predict_data <- data.frame(temp_seq)
  colnames(predict_data) <- variable
  
  # Adding other predictors as constants (mean or typical value)
  predict_data <- predict_data %>%
    mutate(log_mean_StreamLvl_m = mean(data$log_mean_StreamLvl_m, na.rm = TRUE),
           log_mean_StreamFlow_megalitre = mean(data$log_mean_StreamFlow_megalitre, na.rm = TRUE),
           MeanTemp_degC_2yo = mean(data$MeanTemp_degC_2yo, na.rm = TRUE),
           MeanTemp_degC_1yo = mean(data$MeanTemp_degC_1yo, na.rm = TRUE),
           LongestZeroFlow_2yo = mean(data$LongestZeroFlow_2yo, na.rm = TRUE))
  
  # Apply polynomial terms
  predict_data[[paste0("poly(", variable, ", 2)1")]] <- poly(predict_data[[variable]], 2)[,1]
  predict_data[[paste0("poly(", variable, ", 2)2")]] <- poly(predict_data[[variable]], 2)[,2]
  
  # Model prediction
  predictions <- predict(model, newdata = predict_data, se = TRUE)
  predict_data$fit <- predictions$fit
  predict_data$fit_up <- predict_data$fit + 1.96 * predictions$se.fit
  predict_data$fit_low <- predict_data$fit - 1.96 * predictions$se.fit
  
  # Plot
  ggplot(predict_data, aes_string(x = variable, y = "fit")) +
    geom_line() +
    geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2) +
    labs(x = x_label, y = "Predicted Annual Growth Rate", title = plot_title) +
    theme_minimal()
}



# Check the number of unique values
unique_values <- length(unique(fishData_clean$log_mean_StreamFlow_megalitre))
print(unique_values)

generate_plot <- function(variable, model, data, x_label, plot_title) {
  temp_seq <- seq(min(data[[variable]], na.rm = TRUE), max(data[[variable]], na.rm = TRUE), length.out = 100)
  predict_data <- data.frame(temp_seq)
  colnames(predict_data) <- variable
  
  # Adding other predictors as constants (mean or typical value)
  predict_data <- predict_data %>%
    mutate(log_mean_StreamLvl_m = mean(data$log_mean_StreamLvl_m, na.rm = TRUE),
           log_mean_StreamFlow_megalitre = mean(data$log_mean_StreamFlow_megalitre, na.rm = TRUE),
           MeanTemp_degC_2yo = mean(data$MeanTemp_degC_2yo, na.rm = TRUE),
           MeanTemp_degC_1yo = mean(data$MeanTemp_degC_1yo, na.rm = TRUE),
           LongestZeroFlow_2yo = mean(data$LongestZeroFlow_2yo, na.rm = TRUE))
  
  # Apply polynomial terms if possible
  if(length(unique(data[[variable]])) >= 3) {
    predict_data[[paste0("poly(", variable, ", 2)1")]] <- poly(predict_data[[variable]], 2)[,1]
    predict_data[[paste0("poly(", variable, ", 2)2")]] <- poly(predict_data[[variable]], 2)[,2]
  }
  
  # Model prediction
  predictions <- predict(model, newdata = predict_data, se = TRUE)
  predict_data$fit <- predictions$fit
  predict_data$fit_up <- predict_data$fit + 1.96 * predictions$se.fit
  predict_data$fit_low <- predict_data$fit - 1.96 * predictions$se.fit
  
  # Plot
  ggplot(predict_data, aes_string(x = variable, y = "fit")) +
    geom_line() +
    geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.2) +
    labs(x = x_label, y = "Predicted Annual Growth Rate", title = plot_title) +
    theme_minimal()
}

# Generate plot for poly(MeanTemp_degC_spawning, 2)1
p1 <- generate_plot("MeanTemp_degC_spawning", m7_glm_poly_gp, fishData_clean, "Polynomial Term for Mean Temperature at Spawning (°C)", "Effect of Polynomial Term for Mean Temperature at Spawning on Annual Growth Rate")
print(p1)

# Generate plot for log_mean_StreamFlow_megalitre
p2 <- generate_plot("log_mean_StreamFlow_megalitre", m7_glm_poly_gp, fishData_clean, "Log Mean Stream Flow (megalitre)", "Effect of Log Mean Stream Flow on Annual Growth Rate")
print(p2)







######################################################################
## Model Validation

# Load necessary libraries
library(caret)
library(brms)

# Function for cross-validation of GLM models
k_fold_cv_glm <- function(model, data, k = 10) {
  folds <- createFolds(data$Annual_Growth_Rate, k = k, list = TRUE, returnTrain = FALSE)
  cv_results <- lapply(folds, function(test_indices) {
    train_indices <- setdiff(1:nrow(data), test_indices)
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    # Fit model
    fit <- update(model, data = train_data)
    
    # Predict and compute error metrics
    predictions <- predict(fit, newdata = test_data)
    actuals <- test_data$Annual_Growth_Rate
    
    # Compute RMSE and MAE
    rmse <- sqrt(mean((predictions - actuals)^2))
    mae <- mean(abs(predictions - actuals))
    
    return(c(RMSE = rmse, MAE = mae))
  })
  
  # Combine results
  cv_results <- do.call(rbind, cv_results)
  return(cv_results)
}

# Function for cross-validation of Bayesian models
k_fold_cv_bayesian <- function(formula, data, k = 10, family = gaussian(), iter = 2000, chains = 4) {
  folds <- createFolds(data$Annual_Growth_Rate, k = k, list = TRUE, returnTrain = FALSE)
  cv_results <- lapply(folds, function(test_indices) {
    train_indices <- setdiff(1:nrow(data), test_indices)
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    # Fit model
    fit <- brm(formula, data = train_data, family = family, iter = iter, chains = chains, control = list(adapt_delta = 0.95))
    
    # Predict and compute error metrics
    predictions <- as.numeric(predict(fit, newdata = test_data))
    actuals <- test_data$Annual_Growth_Rate
    
    # Compute RMSE and MAE
    rmse <- sqrt(mean((predictions - actuals)^2))
    mae <- mean(abs(predictions - actuals))
    
    return(c(RMSE = rmse, MAE = mae))
  })
  
  # Combine results
  cv_results <- do.call(rbind, cv_results)
  return(cv_results)
}

# Define formulas for Bayesian models
formula_bb <- bf(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + poly(LongestZeroFlow_2yo, 2))
formula_cc <- bf(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + poly(LongestZeroFlow_2yo, 2))
formula_gp <- bf(Annual_Growth_Rate ~ log_mean_StreamLvl_m + log_mean_StreamFlow_megalitre + poly(MeanTemp_degC_2yo, 2) + poly(MeanTemp_degC_1yo, 2) + poly(MeanTemp_degC_spawning, 2) + poly(LongestZeroFlow_2yo, 2))

# Define datasets for each species
bony_bream_data <- fishData_clean %>% filter(Species == "BonyBream")
common_carp_data <- fishData_clean %>% filter(Species == "CommonCarp")
golden_perch_data <- fishData_clean %>% filter(Species == "GoldenPerch")

# Perform k-fold cross-validation for each model

# Bony Bream GLM and Bayesian GLM
m7_glm_poly_bb_cv <- k_fold_cv_glm(m7_glm_poly_bb, bony_bream_data, k = 10)
m7_bayesian_bb_cv <- k_fold_cv_bayesian(formula_bb, bony_bream_data, k = 10)

# Common Carp GLM and Bayesian GLM
m7_glm_poly_cc2_cv <- k_fold_cv_glm(m7_glm_poly_cc2, common_carp_data, k = 10)
m7_bayesian_cc_cv <- k_fold_cv_bayesian(formula_cc, common_carp_data, k = 10)

# Golden Perch GLM and Bayesian GLM
m7_glm_poly_gp_cv <- k_fold_cv_glm(m7_glm_poly_gp, golden_perch_data, k = 10)
m7_bayesian_gp_cv <- k_fold_cv_bayesian(formula_gp, golden_perch_data, k = 10)

# Aggregate results and compare
results <- list(
  BonyBream_GLM = colMeans(m7_glm_poly_bb_cv),
  BonyBream_Bayesian = colMeans(m7_bayesian_bb_cv),
  CommonCarp_GLM = colMeans(m7_glm_poly_cc2_cv),
  CommonCarp_Bayesian = colMeans(m7_bayesian_cc_cv),
  GoldenPerch_GLM = colMeans(m7_glm_poly_gp_cv),
  GoldenPerch_Bayesian = colMeans(m7_bayesian_gp_cv)
)

print(results)

