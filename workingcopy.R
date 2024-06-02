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

#library(LDdiag) 
#also removed from CRAN repository. devtools required to install it from the archives
#devtools::install_version("LDdiag","0.1")
library(LDdiag)



#########################################################################
## Importing and cleaning data from La Trobe U's Microsoft Power BI Solution


## Importing data downloaded from Microsoft Power BI
d1 <- read_excel("./Data/Bonybream.xlsx")
d2 <- read_excel("./Data/Carp.xlsx")
d3 <- read_excel("./Data/Goldenperch.xlsx")
strmflw <- read_excel("./Data/Streamflow.xlsx")
Gauges <- read_excel("./Data/Gauges.xlsx")


# Checking the general structure of the datasets
head(d1)
head(d2)
head(d3)
head(strmflw)
head(Gauges)

## Merging species datasets to form a combined dataset for all 3 species
comb_data <- rbind(d1, d2, d3)

## Incorporating environmental data from Power BI into the species dataset (comb_data)

# Checking the number of unique GaugeID entries
length(unique(strmflw$GaugeID))
# 10 entries. So 10 gauges for 11 sites
# There are 2 gauges for Site 7. 
# Excluding records with GaugeID '422209A' , the second gauge
Gauges <- Gauges %>%
  filter(!(GaugeID == "422209A" & SiteID == 7))


# Changing 'GaugeID' columns from both strmflw and Gauges datasets to the same type
# Because this column will be used to anchor the merging
strmflw$GaugeID <- as.character(strmflw$GaugeID)
Gauges$GaugeID <- as.character(Gauges$GaugeID)

# Merging the dataframes using the 'GaugeID' column
strmflw <- merge(strmflw, Gauges, by = "GaugeID", all.x = TRUE)

# Checking updated strmflw dataset to confirm the merge
head(strmflw)

# Making sure 'Year' and 'SiteID' columns are of the same type in both strmflw and comb_data datasets
comb_data$Year <- as.numeric(as.character(comb_data$Year))
comb_data$SiteID <- as.numeric(as.character(comb_data$SiteID))
strmflw$Year <- as.numeric(as.character(strmflw$Year))
strmflw$SiteID <- as.numeric(as.character(strmflw$SiteID))

# Merging the strmflw and comb_data datasets
full_data <- merge(comb_data, strmflw, by = c("Year", "SiteID"), all.x = TRUE)

# Checking the structure and the first few rows of the merged dataset
str(full_data)
head(full_data)

# Modifying FishID to remove the last three characters (year identifier)
# In case the variable is needed for comparison between years
full_data <- full_data %>%
  mutate(FishID_normalized = str_sub(FishID, 1, -4))


#########################################################################
## Incorporating environmental data from WMIP into full_data


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
  
  # Reading CSV files, skipping the first 3 rows and assuming the 4th row contains headers
  df <- read_csv(file_path, skip = 3, col_types = cols(
    `Date and time` = col_character(),
    Mean = col_double(),
    Min = col_double(),
    Max = col_double()
  ))
  
  df <- df %>%
    mutate(
      Date = dmy(sub("00:00:00 ", "", `Date and time`)),  # Remove the time part and convert to date
      Mean_level_m = Mean,
      Min_level_m = Min,
      Max_level_m = Max,
      GaugeID = gauge_id
    ) %>%
    dplyr::select(Date, Mean_level_m, Min_level_m, Max_level_m, GaugeID)
  
  data_frames[[gauge_id]] <- df
}

# Combine all data frames into one
comb_flow <- bind_rows(data_frames)


# Checking head/structure of combined data frame
head(comb_flow)
str(comb_flow)

# Filtering to retain only stream level data from 2019-2021
comb_flow <- comb_flow %>%
  filter(year(Date) >= 2019 & year(Date) <= 2021)


## Calculating the number of zero flow days and longest zero flow periods

# Summarising full_data dataset to make unique rows for each GaugeID
cease_to_flow <- full_data %>% 
  dplyr::select(GaugeID, Cease_to_flow_level) %>%
  distinct()

# Merging cease_to_flow information with `comb_flow`
comb_flow <- comb_flow %>%
  left_join(cease_to_flow, by = "GaugeID")

# Calculating the number of zero flow days per year for each gauge
zero_flow_days <- comb_flow %>%
  group_by(GaugeID, Year = year(Date)) %>%
  summarize(
    ZeroFlowDays = sum(Mean_level_m <= Cease_to_flow_level, na.rm = TRUE),
    .groups = 'drop'
  )


# Calculating longest no-flow period in each year

# Converting dataset to a data.table
setDT(comb_flow)

# Adding a logical column to identify no-flow days and extracting the year
comb_flow[, `:=`(No_flow = Mean_level_m <= Cease_to_flow_level, Year = year(Date))]

# Calculating the longest no-flow spell by gauge and year
longest_no_flow_spell <- comb_flow[, {
  # Filtering for no-flow days only
  no_flow_days <- .SD[No_flow == TRUE, .(Date)]
  
  # Calculating lengths of consecutive no-flow periods (in days)
  if (nrow(no_flow_days) > 0) {
    # Create a diff column to find breaks in no-flow periods
    no_flow_days[, diff_date := c(1, diff(Date))]
    
    # Each new spell starts when the diff is more than 1
    no_flow_days[, spell_id := cumsum(diff_date > 1)]
    
    # Calculating the length of each spell
    spell_lengths <- no_flow_days[, .(length = .N), by = spell_id][, max(length, na.rm = TRUE)]
  } else {
    spell_lengths = 0  # If no no-flow days, the longest spell is set to 0
  }
  
  # Return the longest spell for each group
  .(LongestZeroFlow = as.double(spell_lengths))
}, by = .(GaugeID, Year)]


# Check resultant dataset
head(longest_no_flow_spell)


# Merging ZeroFlowDays from zero_flow_days into full_data
full_data <- merge(full_data, zero_flow_days, by = c("GaugeID", "Year"), all.x = TRUE)

# Merge LongestZeroFlow from longest_no_flow_spell into full_data
full_data <- merge(full_data, longest_no_flow_spell, by = c("GaugeID", "Year"), all.x = TRUE)

# Checking the structure and head to confirm the merge
str(full_data)
head(full_data)



#########################################################################
## Incorporating environmental data from SILO into full_data


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

# Extracting data year from the YYYY.MM.DD column
comb_gauge <- comb_gauge %>%
  mutate(Year = year(ymd(YYYY.MM.DD)))

# Calculating annual averages for select variables, for each gauge
comb_gauge1 <- comb_gauge %>%
  group_by(GaugeID, Year) %>%
  summarize(
    AvgEvap_pan = mean(evap_pan, na.rm = TRUE),
    AvgDailyRain = mean(daily_rain, na.rm = TRUE),
    MaxTemp_degC_2021 = mean(max_temp, na.rm = TRUE),
    MinTemp_degC_2021 = mean(min_temp, na.rm = TRUE),
    MeanTemp = (mean(max_temp, na.rm = TRUE) + mean(min_temp, na.rm = TRUE)) / 2  # Calculating the mean annual temperature
  ) %>%
  ungroup()

# Adding SiteID column to comb_gauge1, to fascilitate merge with full_data
comb_gauge1 <- comb_gauge1 %>%
  left_join(Gauges %>% dplyr::select(GaugeID, SiteID), by = "GaugeID")

# Merging environmental data from SILO with full_data dataset
full_data <- full_data %>%
  left_join(dplyr::select(comb_gauge1, -GaugeID), by = c("Year", "SiteID"))

# Extracting average annual temperature per site, per year
annual_temps <- comb_gauge1 %>%
  filter(Year %in% c(2019, 2020)) %>%
  group_by(GaugeID, Year) %>%
  summarize(MeanTemp_degC = mean(MeanTemp, na.rm = TRUE), .groups = "drop")
# Last step to make sure there is only a single value per year per gauge

# Add Site ID columns to annual_temps dataset
annual_temps <- annual_temps %>%
  left_join(Gauges %>% dplyr::select(GaugeID, SiteID), by = "GaugeID")


annual_temps1 <- annual_temps %>%
  dplyr::select(GaugeID, SiteID, Year, MeanTemp_degC)

# Merge 2020 and 2019 temperature data

full_data <- full_data %>%
  left_join(dplyr::select(annual_temps, SiteID, MeanTemp_degC, Year) %>% filter(Year == 2020),
            by = "SiteID") %>%
  rename(MeanTemp_degC_2020 = MeanTemp_degC) %>%
  left_join(dplyr::select(annual_temps, SiteID, MeanTemp_degC, Year) %>% filter(Year == 2019),
            by = "SiteID") %>%
  rename(MeanTemp_degC_2019 = MeanTemp_degC)


# Checking resultant dataframe
head(full_data)
str(full_data)
# Still duplicates year column. Have tried multiple options.
# Only alternative left is to remove post merge

# Merge 2020 and 2019 zero flow days data
full_data <- full_data %>%
  left_join(zero_flow_days %>% filter(Year == 2019), by = "GaugeID", suffix = c("", "_2019")) %>%
  rename(ZeroFlowDays_2019 = ZeroFlowDays_2019) %>%
  left_join(zero_flow_days %>% filter(Year == 2020), by = "GaugeID", suffix = c("", "_2020")) %>%
  rename(ZeroFlowDays_2020 = ZeroFlowDays_2020)

# Merge 2020 and 2019 longest no flow spell data
full_data <- full_data %>%
  left_join(longest_no_flow_spell %>% filter(Year == 2019), by = "GaugeID", suffix = c("", "_2019")) %>%
  rename(LongestZeroFlow_2019 = LongestZeroFlow_2019) %>%
  left_join(longest_no_flow_spell %>% filter(Year == 2020), by = "GaugeID", suffix = c("", "_2020")) %>%
  rename(LongestZeroFlow_2020 = LongestZeroFlow_2020)


# Adding a Species column based on FishID_normalized prefixes
full_data <- full_data %>%
  mutate(Species = gsub("[^A-Za-z]", "", FishID_normalized))

# Checking the structure of the updated dataset
str(full_data)
head(full_data)


# Rearranging, renaming and selecting variables from full_data for clarity
full_data_final <- full_data %>%
  dplyr::select(
    Year = Year.x,
    SiteID,
    GaugeID,
    River,
    FishID_normalized,
    Age,
    Annual_Growth_Rate,
    Species,
    mean_StreamLvl_m,
    min_StreamLvl_m,
    max_StreamLvl_m,
    mean_StreamFlow_megalitre,
    min_StreamFlow_megalitre,
    max_StreamFlow_megalitre,
    MeanTemp_degC_2019,
    MeanTemp_degC_2020,
    MeanTemp_degC_2021 = MeanTemp,
    Catchment_area_sq_kms = `Catchment area_sq.kms`,
    Stream_Distance_from_station_to_mouth_km = `Stream Distance from station to mouth_km`,
    LongestZeroFlow_2019,
    LongestZeroFlow_2020,
    LongestZeroFlow_2021 = LongestZeroFlow,
    ZeroFlowDays_2019,
    ZeroFlowDays_2020,
    ZeroFlowDays_2021 = ZeroFlowDays,
    AvgEvap_pan,
    AvgDailyRain
  )

## Filtering data to only include fish in the correct age range (2yo in 2021)
fishData <- full_data_final %>%
  filter(Year == 2021, Age == 2)



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
  select(7, 9:ncol(.)) %>%
  cor(use = "complete.obs")

# Creating a mixed correlation plot

corrplot.mixed(
  corr_matrix,
  tl.cex = 0.01,  # Size of text labels
  number.cex = 0.7,  # Size of correlation coefficients
  upper = "ellipse",
  color = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
  tl.col = "white"
)

# Tabulating the corrmatrix for ease of reference
corr_data <- as.data.frame(corr_matrix) %>%
  rownames_to_column(var = "Variable1") %>%
  pivot_longer(cols = -Variable1, names_to = "Variable2", values_to = "Correlation") %>%
  filter(Variable1 != Variable2) 

## Scatterplot to look at spread of data and 'relationships'

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(GGally)

# Assigning colours to species
species_colors <- c("BB" = "lightblue", "GP" = "salmon", "CC" = "grey")


# Ensuring the dataset 'fishData' contains the 'Species' column used for coloring
plot_data <- fishData %>%
  select(7:ncol(.))

# Check names and adjust formula as needed based on actual names in 'plot_data'
names(plot_data)

plot_data <- plot_data %>%
  na.omit()
plot_data <- plot_data %>%
  filter(!is.na(Species))

# Creating the scatter plot matrix
pairs(~ Annual_Growth_Rate + mean_StreamLvl_m + min_StreamLvl_m + max_StreamLvl_m +
        mean_StreamFlow_megalitre + min_StreamFlow_megalitre + max_StreamFlow_megalitre +
        MeanTemp_degC_2019 + MeanTemp_degC_2020 + MeanTemp_degC_2021 +
        Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km +
        LongestZeroFlow_2019 + LongestZeroFlow_2020 + LongestZeroFlow_2021 +
        ZeroFlowDays_2019 + ZeroFlowDays_2020 + ZeroFlowDays_2021 + AvgEvap_pan + AvgDailyRain,
      data = plot_data,
      col = species_colors[plot_data$Species],
      pch = 16)

# Scatterplot matrix split into three sections for ease of visual examination
pairs(~ Annual_Growth_Rate + mean_StreamLvl_m + min_StreamLvl_m + max_StreamLvl_m +
        mean_StreamFlow_megalitre + min_StreamFlow_megalitre + max_StreamFlow_megalitre,
      data = plot_data,
      col = species_colors[plot_data$Species],
      pch = 16,
      main = "Scatter Plot Matrix for Stream Level and Flow Variables")

pairs(~ Annual_Growth_Rate + MeanTemp_degC_2019 + MeanTemp_degC_2020 + MeanTemp_degC_2021 +
        Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km +
        LongestZeroFlow_2019 + LongestZeroFlow_2020 + LongestZeroFlow_2021,
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

# Adding a legend
legend("topright", legend = c("Bony Bream", "Golden Perch", "Common Carp"), col = c("lightblue", "salmon", "grey"), pch = 16)

#########################################################################
## Exploring the structure of independent variables

## Creating histograms and box plots for continuous predictor variables

# List of continuous variables to plot
continuous_vars <- c("Annual_Growth_Rate", "mean_StreamLvl_m", "min_StreamLvl_m", "max_StreamLvl_m",
                     "mean_StreamFlow_megalitre", "min_StreamFlow_megalitre", "max_StreamFlow_megalitre",
                     "MeanTemp_degC_2019", "MeanTemp_degC_2020", "MeanTemp_degC_2021",
                     "Catchment_area_sq_kms", "Stream_Distance_from_station_to_mouth_km",
                     "LongestZeroFlow_2019", "LongestZeroFlow_2020", "LongestZeroFlow_2021",
                     "ZeroFlowDays_2019", "ZeroFlowDays_2020", "ZeroFlowDays_2021", "AvgEvap_pan", "AvgDailyRain")

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
    labs(title = paste("Boxplot of", var, "by Species"), x = "Species", y = var) +
    theme_minimal()
  print(p)
}

## Bar graph for Species
# Assigning colors to species
species_colors <- c("CC" = "grey", "GP" = "salmon", "BB" = "lightblue")

# Mapping for species names
species_names <- c("CC" = "Common Carp", "GP" = "Golden Perch", "BB" = "Bony Bream")

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

#########################################################################
## Model Fitting
#########################################################################

#########################################################################
## Model 1: Basic model
m1_fullDat <- lm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre + 
           Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km + 
           MeanTemp_degC_2021 + LongestZeroFlow_2021 + Species, data = full_data_final)

# This initial base model uses the least correlated stream variables and
# temperature/ zero flow variables from the study year.
summary(m1_fullDat)
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m1_fullDat))
# testing for outliers
testOutliers(simulateResiduals(m1_fullDat), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m1_fullDat))
# Extracting model equation
extract_eq(m1_fullDat)


## Basic model with only 2 yo fish in 2021: 
m1_fishDat <- lm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre + 
                   Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km + 
                   MeanTemp_degC_2021 + LongestZeroFlow_2021 + Species, data = fishData)
# This initial base model uses the full suite of predictor variables
summary(m1_fishDat)
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m1_fishDat))
# testing for outliers
testOutliers(simulateResiduals(m1_fishDat), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m1_fishDat))


#########################################################################
# Model 2: log transforming variables that are skewed 

# Check skewness of variables
# Removing rows with missing values for the skewness calculation
fishData_clean <- fishData %>%
  filter(!is.na(max_StreamLvl_m) & !is.na(max_StreamFlow_megalitre) & !is.na(LongestZeroFlow_2021))

# Calculating skewness
library(e1071)
skewness_data_clean <- sapply(fishData_clean[, c("Annual_Growth_Rate", "max_StreamLvl_m", "max_StreamFlow_megalitre", 
                                                 "Catchment_area_sq_kms", "Stream_Distance_from_station_to_mouth_km", 
                                                 "MeanTemp_degC_2021", "LongestZeroFlow_2021")], skewness)
print(skewness_data_clean)


# Ensure the fishData_clean is updated with selected variables
fishData <- fishData %>%
  mutate(log_mean_StreamFlow_megalitre = log1p(mean_StreamFlow_megalitre),
         log_Stream_Distance_from_station_to_mouth_km = log1p(Stream_Distance_from_station_to_mouth_km),
         log_LongestZeroFlow_2021 = log1p(LongestZeroFlow_2021))

# Fitting the model with selected and log-transformed variables
m2_transformed <- lm(Annual_Growth_Rate ~ mean_StreamLvl_m + log_mean_StreamFlow_megalitre + 
                       MeanTemp_degC_2021 + Catchment_area_sq_kms +
                       log_Stream_Distance_from_station_to_mouth_km + 
                       log_LongestZeroFlow_2021 + Species, 
                     data = fishData)
summary(m2_transformed)

# Diagnostics
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m2_transformed))
# testing for outliers
testOutliers(simulateResiduals(m2_transformed), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m2_transformed))

#No significant improvement from the un-transformed model


#########################################################################
# Model 3: Introducing interaction terms
# Because of limited success of m1/m2 in explaining variance (low R-squared)

m3_interaction <- lm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre + 
                       Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km + 
                       MeanTemp_degC_2021 * Species + LongestZeroFlow_2021 * Species + 
                       MeanTemp_degC_2021 * LongestZeroFlow_2021, data = fishData)

summary(m3_interaction)

# Diagnostics
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m3_interaction))
# testing for outliers
testOutliers(simulateResiduals(m3_interaction), type = "bootstrap")
# testing for over dispersion
testDispersion(simulateResiduals(m3_interaction))


#########################################################################
# Model 4: Generalised Linear Model with standard identity link

m4_glm <- glm(Annual_Growth_Rate ~ mean_StreamLvl_m + mean_StreamFlow_megalitre + 
                Catchment_area_sq_kms + Stream_Distance_from_station_to_mouth_km + 
                MeanTemp_degC_2021 + LongestZeroFlow_2021 + Species, 
              family = gaussian(), data = fishData)

summary(m4_glm)

#Diagnostics
# testing for Homoscedasticity
plot(simulateResiduals(m4_glm))
## Goodness of fit testing via pregibon
pregibon.glm(m4_glm)


# Testing simplified versions of the GLM
m4_glm_2021 <- glm(Annual_Growth_Rate ~ MeanTemp_degC_2021 + LongestZeroFlow_2021 + Species, 
              family = gaussian(), data = fishData)

summary(m4_glm_2021)

#Diagnostics
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m4_glm_2021))
## Goodness of fit testing via pregibon
pregibon.glm(m4_glm_2021)

m4_glm_2020 <- glm(Annual_Growth_Rate ~ MeanTemp_degC_2020 + LongestZeroFlow_2020 + Species, 
                   family = gaussian(), data = fishData)

summary(m4_glm_2020)

#Diagnostics
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m4_glm_2020))
## Goodness of fit testing via pregibon
pregibon.glm(m4_glm_2020)

m4_glm_2019 <- glm(Annual_Growth_Rate ~ MeanTemp_degC_2019 + LongestZeroFlow_2019 + Species, 
                   family = gaussian(), data = fishData)
summary(m4_glm_2019)

#Diagnostics
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m4_glm_2019))
## Goodness of fit testing via pregibon
pregibon.glm(m4_glm_2019)

# Simplified model with only Species
m4_glm_species <- glm(Annual_Growth_Rate ~ Species, family = gaussian(), data = fishData)
summary(m4_glm_species)

#Diagnostics
# testing for Homoscedasticity and Normality
plot(simulateResiduals(m4_glm_species))
## Goodness of fit testing via pregibon
pregibon.glm(m4_glm_species)

#########################################################################
# Model 5: PCA to reduce dimension

# Selecting the full suite of predictors
predictors <- fishData[, c("mean_StreamLvl_m", "min_StreamLvl_m", "max_StreamLvl_m", 
                           "mean_StreamFlow_megalitre", "min_StreamFlow_megalitre", 
                           "max_StreamFlow_megalitre", "MeanTemp_degC_2019", 
                           "MeanTemp_degC_2020", "MeanTemp_degC_2021", 
                           "Catchment_area_sq_kms", 
                           "Stream_Distance_from_station_to_mouth_km", 
                           "LongestZeroFlow_2019", "LongestZeroFlow_2020", 
                           "LongestZeroFlow_2021", "ZeroFlowDays_2019", 
                           "ZeroFlowDays_2020", "ZeroFlowDays_2021", 
                           "AvgEvap_pan", "AvgDailyRain")]

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
fishData_clean <- fishData[complete.cases(predictors_clean), ]

# Selecting the cleaned predictors
predictors_clean <- fishData_clean[, names(predictors_clean)]

# Applying PCA to the cleaned predictors
pca <- prcomp(predictors_clean, scale. = TRUE)

# Check the proportion of variance explained by each principal component
summary(pca)

# Use the first few principal components that explain most of the variance
fishData_pca <- data.frame(pca$x[, 1:3]) # Adjust the number of PCs as needed
fishData_pca$Annual_Growth_Rate <- fishData_clean$Annual_Growth_Rate
fishData_pca$Species <- fishData_clean$Species

# Fitting the model with principal components
m5_pca <- lm(Annual_Growth_Rate ~ PC1 + PC2 + PC3 + Species, data = fishData_pca)
summary(m5_pca)

# Diagnostics
plot(simulateResiduals(m5_pca))

pregibon.lm(m5_pca)


#########################################################################
# Model 6: GLM with PCA

# Fitting the model with principal components
m6_pca <- glm(Annual_Growth_Rate ~ PC1 + PC2 + PC3 + Species, data = fishData_pca)
summary(m6_pca)

# Diagnostics
plot(simulateResiduals(m6_pca))

pregibon.glm(m6_pca)



#########################################################################
# Model 7: Bayesian GLM

#Checking where g++ and make are located
Sys.which("g++")
Sys.which("make")

# Fitting model with PCs
m7_bayesian <- brm(
  Annual_Growth_Rate ~ PC1 + PC2 + PC3 + Species,
  data = fishData_pca,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),      # Priors for regression coefficients
    prior(student_t(3, 0, 10), class = "Intercept")  # Prior for the intercept
  ),
  iter = 4000, chains = 4, cores = 4,  # Adjust iterations and chains as needed
  control = list(adapt_delta = 0.95)  # Adjust control parameters if necessary
)

# Model summary
summary(m7_bayesian)

# Plot the posterior distributions
plot(m7_bayesian)

# Posterior predictive checks
pp_check(m7_bayesian)


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
# Extracting all model equations

extract_eq(m1_fullDat)
extract_eq(m1_fishDat)
extract_eq(m2_transformed)
extract_eq(m3_interaction)
extract_eq(m4_glm)
extract_eq(m4_glm_2019)
extract_eq(m4_glm_2020)
extract_eq(m4_glm_2021)
extract_eq(m4_glm_species)
extract_eq(m5_pca)
extract_eq(m6_pca)
extract_eq(m7_bayesian)


