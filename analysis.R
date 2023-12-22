#######################################################################
# ASV SPATIOTEMPORAL HETEROGENEITY -- MAIN TEXT
# Script written by Arsenault, Steele, Shingai, Cottingham et al. 2023
# Produces data analyses and figures for manuscript (main text)
# Relies on "setup.R" script, which loads libraries and data
#######################################################################

# Load libraries

library(tidyverse)
library(timetk)
library(lubridate)
library(zoo)
library(proj4)
library(rgdal)
library(rollply)
library(remotes)
library(mcr)
library(cowplot)
library(lemon)
library(gridExtra)
library(ggpubr)
library(ggmap)
library(forcats)
library(readr)
library(RColorBrewer)
library(ggnewscale)
library(Hmisc)
library(PerformanceAnalytics)
library(scales)
library(ggsn)
library(sf)
library(colorspace)

# Load data files

aub30Aug21 <- read.csv(file = "AUB/AUB_2021-08-30_a_asv_processed_v2022-03-23.csv", header=T)

chn06Aug21 <- read.csv(file = "CHN/CHN_2021-08-06_a_asv_processed_v2022-03-23.csv", header = T)[1:4005,] #modified to cut off after last loiter

chn28Sep21 <- read.csv(file = "CHN/CHN_2021-09-28_a_asv_processed_v2022-03-23.csv", header = T) %>% 
  select(-specificConductance_uscm) %>% #eliminate SpC data due to unsuccessful calibration
  mutate(specificConductance_uscm = NA)

chn05Oct21 <- read.csv(file = "CHN/CHN_2021-10-05_a_asv_processed_v2022-03-23.csv", header = T)[1:2128,] %>% #trim after last programmed waypoint
  select(-specificConductance_uscm) %>% #eliminate SpC data due to unsuccessful calibration
  mutate(specificConductance_uscm = NA)

chn12Oct21 <- read.csv(file = "CHN/CHN_2021-10-12_a_asv_processed_v2022-03-23.csv", header = T)

chn21Oct21 <- read.csv(file = "CHN/CHN_2021-10-21_a_asv_processed_v2022-03-23.csv", header = T)

sab20Aug21 <- read.csv(file = "SAB/SAB_2021-08-20_a_asv_processed_v2022-03-23.csv", header = T)

sab26Aug21 <- read.csv(file = "SAB/SAB_2021-08-26_a_asv_processed_v2022-03-23.csv", header = T)[697:15177,] #modified to trim out test path at beginning

sun11Jun21.HC <- read.csv(file = "SUN/SUN_2021-06-11_HC_asv_processed_v2022-03-23.csv", header = T)

sun11Jun21.NW <- read.csv(file = "SUN/SUN_2021-06-11_NW_asv_processed_v2022-03-23.csv", header = T)

sun15Jun21.HC <- read.csv(file = "SUN/SUN_2021-06-15_HC_asv_processed_v2022-03-23.csv", header = T)

sun15Jun21.NW <- read.csv(file = "SUN/SUN_2021-06-15_NW_asv_processed_v2022-03-23.csv", header = T)

sun21Jun21.HC <- read.csv(file = "SUN/SUN_2021-06-21_HC_asv_processed_v2022-03-23.csv", header = T)

sun21Jun21.NW <- read.csv(file = "SUN/SUN_2021-06-21_NW_asv_processed_v2022-03-23.csv", header = T)

sun01Jul21.HC <- read.csv(file = "SUN/SUN_2021-07-01_HC_asv_processed_v2022-03-23.csv", header = T)

sun22Jul21.HC <- read.csv(file = "SUN/SUN_2021-07-22_asv_processed_v2022-03-23.csv", header = T)[1:3222,] #trimmed to include planned path only

sun22Jul21.NW <- read.csv(file = "SUN/SUN_2021-07-22_asv_processed_v2022-03-23.csv", header = T)[3223:10040,] #trimmed to include planned path only

sun05Aug21.HC <- read.csv(file = "SUN/SUN_2021-08-05_HC_asv_processed_v2022-03-23.csv", header = T) 

sun05Aug21.NW <- read.csv(file = "SUN/SUN_2021-08-05_NW_asv_processed_v2022-03-23.csv", header = T)

sun10Aug21.HC <- read.csv(file = "SUN/SUN_2021-08-10_HC_asv_processed_v2022-03-23.csv", header = T)

sun10Aug21.NW <- read.csv(file = "SUN/SUN_2021-08-10_NW_asv_processed_v2022-03-23.csv", header = T)

sun27Aug21.HC <- read.csv(file = "SUN/SUN_2021-08-27_HC_asv_processed_v2022-03-23.csv", header = T)

sun27Aug21.NW <- read.csv(file = "SUN/SUN_2021-08-27_NW_asv_processed_v2022-03-23.csv", header = T) 

sun16Sep21.HC <- read.csv(file = "SUN/SUN_2021-09-16_HC_asv_processed_v2022-03-23.csv", header = T)

sun16Sep21.NW <- read.csv(file = "SUN/SUN_2021-09-16_NW_asv_processed_v2022-03-23.csv", header = T)

# Label data for analysis

travel_effects <- function(df, n_sec){ # requires a dataframe (df) and number of seconds to designate for movement, pre- and post-loiter
  
  # get a list of the loiter waypoints (wp) by HEADER SEQUENCE, waypoint_seq repeats sometimes, and we want each of these to be unique
  wp_loiter <- unique(df$header_seq[df$wp_command == 19 & !is.na(df$wp_command)])   
  
  # loop over way points
  # here 'i' is an iteration along 1 to the length of the value 'wp_loiter'; i therefore refers to the position 1:n of the list 'wp_loiter'
  for (i in 1:length(wp_loiter)){  
    
    # grab the time stamp for the loiter end using the wp_loiter list
    # here, wp_loiter[i] is grabbing the i'th position in the wp_loiter list.
    loiter_end <- df$timestamp_gps_sec[which(df$header_seq == wp_loiter[i])]   
    
    # calculate loiter start from the wp param (seconds of loiter) associated with the loiter
    loiter_start <- loiter_end - df$wp_param1[which(df$header_seq == wp_loiter[i])]
    
    # determine the start and end of the pre and post loiter periods 
    # uses the number of seconds away from the start or end of the loiter controlled by the user as sent into this function
    preloit_start <- loiter_start - n_sec
    postloit_end <- loiter_end + n_sec
    
    # make a data frame that only includes the info we're interested in for this waypoint
    df_test <- df %>% 
      # keep only those timestamps from the preloiter start through the post loiter end
      filter(timestamp_gps_sec >= preloit_start & timestamp_gps_sec < postloit_end)
    
    # calculate and label pre/post/during loiters
    df_test <- df_test %>% 
      mutate(loiter_tag = case_when(timestamp_gps_sec < loiter_start ~ 'pre-loiter',
                                    timestamp_gps_sec >= loiter_start & timestamp_gps_sec <= loiter_end ~ 'loiter',
                                    timestamp_gps_sec > loiter_end ~ 'post-loiter',
                                    TRUE ~ NA_character_))
    
    df_test$loiter_num = i
    df_test$loiterID = paste(df_test$loiter_tag, df_test$loiter_num, sep = '_')
    
    # if the first time through iteration, save test as df; otherwise append!
    if (i == 1) {
      loiter_df = df_test
    } else {
      loiter_df = full_join(loiter_df, df_test)
    }
  }
  return(loiter_df)
}

preset_loiter_time <- 60 # time in seconds pre- and post- loiter = travel time to compare to loiter time

# Run travel_effects function for China Lake - 5 deployments

chn06Aug21_labeled <- travel_effects(chn06Aug21, preset_loiter_time) %>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
chn28Sep21_labeled <- travel_effects(chn28Sep21, preset_loiter_time) %>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
chn05Oct21_labeled <- travel_effects(chn05Oct21, preset_loiter_time) %>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
chn12Oct21_labeled <- travel_effects(chn12Oct21, preset_loiter_time) %>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
chn21Oct21_labeled <- travel_effects(chn21Oct21, preset_loiter_time) %>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")

# Run travel_effects function for Sabattus Pond - 2 deployments

sab20Aug21_labeled <- travel_effects(sab20Aug21, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sab26Aug21_labeled <- travel_effects(sab26Aug21, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")

# Run travel_effects function for Auburn - 1 deployment

aub30Aug21_labeled <- travel_effects(aub30Aug21, preset_loiter_time) %>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")

# Run travel_effects function for Sunapee - 9 deployments

sun11Jun21.HC_labeled <- travel_effects(sun11Jun21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun11Jun21.NW_labeled <- travel_effects(sun11Jun21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun15Jun21.HC_labeled <- travel_effects(sun15Jun21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun15Jun21.NW_labeled <- travel_effects(sun15Jun21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun21Jun21.HC_labeled <- travel_effects(sun21Jun21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun21Jun21.NW_labeled <- travel_effects(sun21Jun21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun01Jul21.HC_labeled <- travel_effects(sun01Jul21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU") #only herrick cove
sun22Jul21.HC_labeled <- travel_effects(sun22Jul21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun22Jul21.NW_labeled <- travel_effects(sun22Jul21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun05Aug21.HC_labeled <- travel_effects(sun05Aug21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun05Aug21.NW_labeled <- travel_effects(sun05Aug21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun10Aug21.HC_labeled <- travel_effects(sun10Aug21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun10Aug21.NW_labeled <- travel_effects(sun10Aug21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun27Aug21.HC_labeled <- travel_effects(sun27Aug21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun27Aug21.NW_labeled <- travel_effects(sun27Aug21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun16Sep21.HC_labeled <- travel_effects(sun16Sep21.HC, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun16Sep21.NW_labeled <- travel_effects(sun16Sep21.NW, preset_loiter_time)%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")

# Avoid repeating site IDs (Sunapee data only) due to having two separate datasets (HC and NW) for each deployment date

sun11Jun21.HC_labeled <- rename(sun11Jun21.HC_labeled, loiternumb_old = loiter_num)
sun15Jun21.HC_labeled <- rename(sun15Jun21.HC_labeled, loiternumb_old = loiter_num)
sun21Jun21.HC_labeled <- rename(sun21Jun21.HC_labeled, loiternumb_old = loiter_num)
sun22Jul21.HC_labeled <- rename(sun22Jul21.HC_labeled, loiternumb_old = loiter_num)
sun05Aug21.HC_labeled <- rename(sun05Aug21.HC_labeled, loiternumb_old = loiter_num)
sun10Aug21.HC_labeled <- rename(sun10Aug21.HC_labeled, loiternumb_old = loiter_num)
sun27Aug21.HC_labeled <- rename(sun27Aug21.HC_labeled, loiternumb_old = loiter_num)
sun16Sep21.HC_labeled <- rename(sun16Sep21.HC_labeled, loiternumb_old = loiter_num)

# Replace loiter numbers 1,2 with 4,5 in the HC dataframes

sun11Jun21.HC_labeled_renamed <- mutate(sun11Jun21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun15Jun21.HC_labeled_renamed <- mutate(sun15Jun21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun21Jun21.HC_labeled_renamed <- mutate(sun21Jun21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun22Jul21.HC_labeled_renamed <- mutate(sun22Jul21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun05Aug21.HC_labeled_renamed <- mutate(sun05Aug21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun10Aug21.HC_labeled_renamed <- mutate(sun10Aug21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun27Aug21.HC_labeled_renamed <- mutate(sun27Aug21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")
sun16Sep21.HC_labeled_renamed <- mutate(sun16Sep21.HC_labeled, loiter_num = ifelse(loiternumb_old == 1, 4, 5))%>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")

# Combine columns from data frames to run analysis by parameter instead of by lake (combine all lakes)

df <- rbind(chn06Aug21_labeled, chn28Sep21_labeled, chn05Oct21_labeled, chn12Oct21_labeled, chn21Oct21_labeled, 
            sab20Aug21_labeled, sab26Aug21_labeled, aub30Aug21_labeled, 
            sun11Jun21.HC_labeled_renamed, sun11Jun21.NW_labeled, sun15Jun21.HC_labeled_renamed, sun15Jun21.NW_labeled, 
            sun21Jun21.HC_labeled_renamed, sun21Jun21.NW_labeled, sun01Jul21.HC_labeled, sun22Jul21.HC_labeled_renamed, sun22Jul21.NW_labeled,
            sun05Aug21.HC_labeled_renamed, sun05Aug21.NW_labeled, sun10Aug21.HC_labeled_renamed, sun10Aug21.NW_labeled, 
            sun27Aug21.HC_labeled_renamed, sun27Aug21.NW_labeled, sun16Sep21.HC_labeled_renamed, sun16Sep21.NW_labeled)

# Define "source" as an as.factor(date)

df$source <- as.factor(df$date)

# Filter the data to ask questions about differences between traveling and loitering

during_loiter_speed_threshold <- 0.5 # how fast the robot can be going in a loiter period
during_mvmt_speed_threshold <- 0.9 # how slow the robot can be going in a movement period

df_trim <- df %>%
  # keep only those rows rows which are during a movement period but gps speed is greater than the threshold
  # OR are during a loiter period & gps speed is below the too fast threshold
  filter(((loiter_tag=="pre-loiter" | loiter_tag=="post-loiter") & velocity_gps_mps > during_mvmt_speed_threshold) |
           (loiter_tag=="loiter" & velocity_gps_mps < during_mvmt_speed_threshold)) 

# Grouping pre/post to be just movement by site ID
df1 <-  df_trim %>%
  mutate(movement=case_when(str_detect(loiterID,"pre-loiter|post-loiter")~paste0("movement_",loiter_num),
                            TRUE~ paste0("loiter_",loiter_num)))

# pull out the site id as a separate column in a "data frame to functions"
df_to_fxns <- df1 %>% separate(movement,c("type","siteid"),"_") %>% mutate(siteid=parse_number(siteid)) 

# make sure everything is there (all lakes, all sites) - this creates a new column that ultimately won't be used for anything else other than to check
df_to_fxns$lake_site_date <- paste(df_to_fxns$lake, df_to_fxns$siteid, df_to_fxns$date)
df_to_fxns$lake_date <- paste(df_to_fxns$lake, df_to_fxns$date)

unique(df_to_fxns$lake_date)
unique(df_to_fxns$lake_site_date)

# Figure 2

plot_param <- function(df4, name, sensorError, label, startlim, endlim) {
  
  df4 <- df4 %>%
    # treat the source column as a date
    mutate(plotdate=as.Date.factor(source)) 
  
  ggplot(data=df4, aes(x=plotdate, y=median.diff, color = lake, fill = lake, shape = lake))+
    geom_hline(yintercept= 0, color = "grey")+
    geom_point(data=df4, aes(x=plotdate, y=median.diff), alpha = 0.6, size = 3)+
    geom_hline(yintercept = sensorError, linetype= "longdash")+
    geom_hline(yintercept= -(sensorError), linetype= "longdash")+
    scale_fill_manual(values = c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
    scale_color_manual(values = c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
    scale_shape_manual(values = c("AUB" = 21, "CHN" = 22, "SAB" = 23, "SUN" = 24))+    
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), axis.text=element_text(size=16, colour="black"), 
          legend.key = element_rect(fill = NA), plot.tag.position = c(0.8, 0.93),
          plot.tag = element_text(size = 16), axis.title = element_text(size=16))+    
    labs(x="", y="\u394 Loiter-Movement",str_replace(deparse(substitute(name)),"quo","")) +
    labs(tag = paste(label))+
    ylim(c(startlim, endlim))
}

# Function to calculate mean and median of loiter and movement periods
center.bydate <- function(df_in,name) {
  
  # calculate the mean travel/loiter at a site on a date, then get the difference between them
  # df_int = intermediate data frame
  df_int <- df_in %>%
    group_by(lake_date, lake, source, siteid, type) %>%
    # get the mean & then ungroup
    summarise(mean.x= mean(!!name)) %>% ungroup()
  
  # create a new data frame that makes separate columns for each type
  df_int_mean <- df_int %>%
    pivot_wider(id_cols=c("lake_date","lake", "source", "siteid"),names_from = type, values_from = mean.x) %>%
    # calculate the difference between loitering and moving at each site
    mutate(mean.loiter=loiter,
           mean.movement=movement,
           mean.diff=loiter-movement,
           mean.pctdiff=mean.diff/loiter) %>%
    select(lake_date, lake, source, siteid, mean.loiter, mean.movement, mean.diff, mean.pctdiff)
  
  # calculate the median travel/loiter at a site on a date, then get the difference between them
  # df_int = intermediate data frame
  df_int <- df_in %>%
    # group by sample date/source, site ID, and type of movement (loiter/travel)
    group_by(lake_date, lake, source, siteid, type) %>%
    # get the mean & then ungroup
    summarise(median.x= median(!!name)) %>% ungroup()
  
  # create a new data frame that makes separate columns for each type
  df_int_median <- df_int %>%
    pivot_wider(id_cols=c("lake_date","lake", "source","siteid"),names_from = type, values_from = median.x) %>%
    # calculate the difference between loitering and moving at each site
    mutate(median.loiter=loiter,
           median.movement=movement,
           median.diff=loiter-movement,
           median.pctdiff=median.diff/loiter) %>%
    select(lake_date, lake, source, siteid, median.loiter, median.movement, median.diff, median.pctdiff)
  
  # combine the means & medians
  df_int_both=full_join(df_int_mean, df_int_median, by=c("lake_date", "lake", "source","siteid"))
  
  # return the updated matrix
  center.bydate <- df_int_both
  
} # end function

# Define sensor accuracy ranges
dissolvedOxygen <- 0.1
turbidity <- 0.3
temperature <- 0.2
specificConductance <- 2.0
pH <- 0.2

chlorophyll_calc <- center.bydate(df_to_fxns, quo(chlorophyll_a_RFU))
chlorophyll_calc <- subset(chlorophyll_calc, (!is.na(chlorophyll_calc$mean.diff) & !is.na(chlorophyll_calc$median.diff)))

chlorophyll <- chlorophyll_calc %>% #chlorophyll gets its own code because sensorError = NA
  mutate(plotdate=as.Date.factor(source))

chlorophyll_plot <- ggplot(data=chlorophyll, aes(x=plotdate, y=median.diff, color = lake, fill = lake, shape = lake))+
  geom_point(data=chlorophyll, aes(x=plotdate, y=median.diff), alpha = 0.6, size = 3) +
  geom_hline(yintercept= 0, color = "grey")+
  scale_fill_manual(values = c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_color_manual(values = c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_shape_manual(values = c("AUB" = 21, "CHN" = 22, "SAB" = 23, "SUN" = 24))+
  theme(legend.position = c(0.5, 0.1), panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), axis.text=element_text(size=16, colour="black"), 
        legend.key = element_rect(fill = NA), plot.tag.position = c(0.8, 0.93),
        plot.tag = element_text(size = 16), axis.title = element_text(size=16),
        legend.direction = "horizontal", legend.title = element_blank(),
        legend.text = element_text(size=12))+   
  labs(x="", y="\u394 Loiter-Movement",str_replace(deparse(substitute(name)),"quo","")) +
  labs(tag = substitute(paste("Chl. ", italic("a "), "(RFU)")))+
  ylim(c(-0.5, 0.5))

dissolvedOxygen_calc <- center.bydate(df_to_fxns, quo(oxygenDissolved_mgl))
dissolvedOxygen_calc <- subset(dissolvedOxygen_calc, (!is.na(dissolvedOxygen_calc$mean.diff) & !is.na(dissolvedOxygen_calc$median.diff)))
dissolvedOxygen_plot <- plot_param(dissolvedOxygen_calc, quo(dissolvedOxygen_mgl), sensorError = dissolvedOxygen, "DO (mg/L)", -0.5, 0.5)

pH_calc <- center.bydate(df_to_fxns, quo(pH))
pH_calc <- subset(pH_calc, (!is.na(pH_calc$mean.diff) & !is.na(pH_calc$median.diff)))
pH_plot <- plot_param(pH_calc, quo(pH), sensorError = pH, "pH", -0.5, 0.5)

specificConductance_calc <- center.bydate(df_to_fxns, quo(specificConductance_uscm))
specificConductance_calc <- subset(specificConductance_calc, (!is.na(specificConductance_calc$mean.diff) & !is.na(specificConductance_calc$median.diff)))
specificConductance_plot <- plot_param(specificConductance_calc, quo(specificConductance_uscm), sensorError = specificConductance, paste0("Sp.C. (", "\u00B5", "S/cm)"), -2, 2)

temperature_calc <- center.bydate(df_to_fxns, quo(temperatureWater_degC))
temperature_calc <- subset(temperature_calc, (!is.na(temperature_calc$mean.diff) & !is.na(temperature_calc$median.diff)))
temperature_plot <- plot_param(temperature_calc, quo(temperatureWater_degC), sensorError = temperature, "Temp. (ºC)", -0.5, 0.5)

turbidity_calc <- center.bydate(df_to_fxns, quo(turbidity_NTU))
turbidity_calc <- subset(turbidity_calc, (!is.na(turbidity_calc$mean.diff) & !is.na(turbidity_calc$median.diff)))

turbidity_calc_plot <- turbidity_calc %>%
  mutate(plotdate=as.Date.factor(source))

turbidity_plot <- ggplot(data=turbidity_calc_plot, aes(x=plotdate, y=median.diff, color = lake, fill = lake, shape = lake))+
  geom_point(data=turbidity_calc_plot, aes(x=plotdate, y=median.diff), alpha = 0.6, size = 3) +
  geom_hline(yintercept= 0, color = "grey")+
  geom_hline(yintercept = turbidity, linetype= "longdash")+
  geom_hline(yintercept= -(turbidity), linetype= "longdash")+
  scale_fill_manual(values = c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_color_manual(values = c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_shape_manual(values = c("AUB" = 21, "CHN" = 22, "SAB" = 23, "SUN" = 24))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), axis.text=element_text(size=16, colour="black"), 
        legend.key = element_rect(fill = NA), plot.tag.position = c(0.8, 0.93),
        plot.tag = element_text(size = 16), axis.title = element_text(size=16),
        legend.direction = "horizontal", legend.title = element_blank(),
        legend.text = element_text(size=12))+   
  labs(x="", y="\u394 Loiter-Movement",str_replace(deparse(substitute(name)),"quo","")) +
  scale_x_date(date_breaks = "1 month", 
               labels=date_format("%b"),
               limits = as.Date(c('2021-06-11','2021-10-21'))) +
  labs(tag = substitute(paste("Turb. (NTU)")))+
  ylim(c(-2.1, 2.1))

# Calculate number outside sensor range of error
dissolvedOxygen_calc_traveleffects <- dissolvedOxygen_calc %>% 
  mutate(pass = case_when(median.diff < -dissolvedOxygen | median.diff > dissolvedOxygen ~ 1,
            TRUE ~ 0))
dissolvedOxygen_calc_traveleffects_fails <- sum(dissolvedOxygen_calc_traveleffects$pass)

pH_calc_traveleffects <- pH_calc %>% 
  mutate(pass = case_when(median.diff < -pH | median.diff > pH ~ 1,
                          TRUE ~ 0))
pH_calc_traveleffects_fails <- fail <- sum(pH_calc_traveleffects$pass)

specificConductance_calc_traveleffects <- specificConductance_calc %>% 
  mutate(pass = case_when(median.diff < -specificConductance | median.diff > specificConductance ~ 1,
                          TRUE ~ 0))
specificConductance_calc_traveleffects_fails <- sum(specificConductance_calc_traveleffects$pass)

temperature_calc_traveleffects <- temperature_calc %>% 
  mutate(pass = case_when(median.diff < -temperature | median.diff > temperature ~ 1,
                          TRUE ~ 0))
temperature_calc_traveleffects_fails <- sum(temperature_calc_traveleffects$pass)

turbidity_calc_traveleffects <- turbidity_calc %>% 
  mutate(pass = case_when(median.diff < -turbidity | median.diff > turbidity ~ 1,
                          TRUE ~ 0))
turbidity_calc_traveleffects_fails <- sum(turbidity_calc_traveleffects$pass)

sum(c(dissolvedOxygen_calc_traveleffects_fails, pH_calc_traveleffects_fails,
    specificConductance_calc_traveleffects_fails, temperature_calc_traveleffects_fails,
    turbidity_calc_traveleffects_fails)) / sum(c(nrow(dissolvedOxygen_calc), 
                                                nrow(pH_calc),
                                                nrow(specificConductance_calc),
                                                nrow(temperature_calc),
                                                nrow(turbidity_calc)))*100

fig2 <- ggarrange(temperature_plot, dissolvedOxygen_plot, pH_plot, specificConductance_plot, turbidity_plot, chlorophyll_plot,
          nrow=2, ncol=3, align = "v", labels = c("a","b","c", "d", "e", "f"))

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/Fig2.png", fig2, height = 10, width = 13)

# Prepare files for next set of plots

sun11Jun21.HC <- sun11Jun21.HC %>% 
  mutate(turbidity_NTU = NA) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg) 

sun11Jun21.NW <- sun11Jun21.NW %>% 
  mutate(turbidity_NTU = NA) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun11Jun21 <- rbind(sun11Jun21.HC, sun11Jun21.NW) %>% 
  mutate(turbidity_NTU = NA) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun15Jun21 <- rbind(sun15Jun21.HC, sun15Jun21.NW) %>% 
  mutate(turbidity_NTU = NA) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun21Jun21 <- rbind(sun21Jun21.HC, sun21Jun21.NW) %>% 
  mutate(turbidity_NTU = NA) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun01Jul21.HC <- sun01Jul21.HC %>% 
  mutate(turbidity_NTU = NA) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun22Jul21 <- rbind(sun22Jul21.HC, sun22Jul21.NW) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun05Aug21 <- rbind(sun05Aug21.HC, sun05Aug21.NW) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun10Aug21 <- rbind(sun10Aug21.HC, sun10Aug21.NW) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun27Aug21.HC <- sun27Aug21.HC %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun27Aug21.NW <- sun27Aug21.NW %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun27Aug21 <- rbind(sun27Aug21.HC, sun27Aug21.NW) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

sun16Sep21 <- rbind(sun16Sep21.HC, sun16Sep21.NW) %>% 
  select(date, lake, chlorophyll_a_RFU, oxygenDissolved_mgl, pH, 
         specificConductance_uscm, temperatureWater_degC, turbidity_NTU, 
         timestamp_sonde_sec, velocity_gps_mps, latitude_gps_deg, longitude_gps_deg)

# Join dataframes by lake for plotting
auburn_all <- aub30Aug21
china_all <- rbind(chn06Aug21, chn28Sep21, chn05Oct21, chn12Oct21, chn21Oct21)
sabattus_all <- rbind(sab20Aug21, sab26Aug21)
sunapee_all <- rbind(sun11Jun21, sun15Jun21, sun21Jun21, sun01Jul21.HC, 
                     sun22Jul21, sun05Aug21, sun10Aug21, sun27Aug21, sun16Sep21)

# Moving window for visualizations

asv_window <- function(df, parameter, window_size_sec){ #parameter should not be entered within quotation marks
  
  parameter <- enquo(parameter)
  
  df2 <- df %>% 
    mutate(date = as.POSIXct(timestamp_sonde_sec, origin = "1970-01-01"))
  
  window <- rollply(df2, ~ date, wdw.size = window_size_sec,
                    dplyr::summarize, mean = mean(!!parameter, na.rm = TRUE),)
  
  return(window)
}

# Set boundaries for maps

auburn <- c(left = -70.2504, bottom = 44.1450, right = -70.2345, top = 44.1600)
china <- c(left = -69.6071, bottom = 44.4430, right = -69.6000, top = 44.4470)
sabattus <- c(left = -70.1180, bottom = 44.1156, right = -70.0780, top = 44.1760)
sunapee <- c(left = -72.0562, bottom = 43.4073, right = -72.0300, top = 43.4200)

# Figure 3

sabattus_large <- c(left = -70.1300, bottom = 44.1156, right = -70.0700, top = 44.1760)

sab_temp <- get_stadiamap(sabattus_large, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$temperatureWater_degC, na.rm = T)-min(sab26Aug21$temperatureWater_degC, na.rm = T))/2)+min(sab26Aug21$temperatureWater_degC, na.rm = T),
                        name="\u00B0C") +
  annotate("point", x = -70.09111, y = 44.17221, colour = "blue", size = 3) + # Hooper Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle("Temp. (\u00B0C)") +
  xlab("Longitude") +
  ylab("Latitude")

aub_temp <- get_stadiamap(auburn, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
             data = aub30Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$temperatureWater_degC, na.rm = T)-min(aub30Aug21$temperatureWater_degC, na.rm = T))/2)+min(aub30Aug21$temperatureWater_degC, na.rm = T),
                        name="\u00B0C") +
  annotate("point", x = -70.24266, y = 44.15925, colour = "blue", size = 3) + # Townsend Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5), plot.tag = element_text(color = "white"))+
  ggtitle("Temp. (\u00B0C)") +
  xlab("Longitude") +
  ylab("Latitude")

sab_spc <- get_stadiamap(sabattus_large, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$specificConductance_uscm, na.rm = T)-min(sab26Aug21$specificConductance_uscm, na.rm = T))/2)+min(sab26Aug21$specificConductance_uscm, na.rm = T),
                        name=paste0("Sp.C. (", "\u00B5", "S/cm)")) +
  annotate("point", x = -70.09111, y = 44.17221, colour = "blue", size = 3) + # Hooper Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle(paste0("Sp.C. (", "\u00B5", "S/cm)")) +
  xlab("Longitude") +
  ylab("Latitude")

aub_spc <- get_stadiamap(auburn, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = aub30Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$specificConductance_uscm, na.rm = T)-min(aub30Aug21$specificConductance_uscm, na.rm = T))/2)+min(aub30Aug21$specificConductance_uscm, na.rm = T),
                        name=paste0("Sp.C. (", "\u00B5", "S/cm)")) +
  annotate("point", x = -70.24266, y = 44.15925, colour = "blue", size = 3) + # Townsend Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5), plot.tag = element_text(color = "white"))+
  ggtitle(paste0("Sp.C. (", "\u00B5", "S/cm)")) +
  xlab("Longitude") +
  ylab("Latitude")

sab_turb <- get_stadiamap(sabattus_large, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$turbidity_NTU, na.rm = T)-min(sab26Aug21$turbidity_NTU, na.rm = T))/2)+min(sab26Aug21$turbidity_NTU, na.rm = T),
                        name="NTU") +
  annotate("point", x = -70.09111, y = 44.17221, colour = "blue", size = 3) + # Hooper Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle("Turb. (NTU)") +
  xlab("Longitude") +
  ylab("Latitude")

aub_turb <- get_stadiamap(auburn, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
             data = aub30Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$turbidity_NTU, na.rm = T)-min(aub30Aug21$turbidity_NTU, na.rm = T))/2)+min(aub30Aug21$turbidity_NTU, na.rm = T),
                        name="NTU") +
  annotate("point", x = -70.24266, y = 44.15925, colour = "blue", size = 3) + # Townsend Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5), plot.tag = element_text(color = "white"))+
  ggtitle("Turb. (NTU)") +
  xlab("Longitude") +
  ylab("Latitude")

sab_chl <- get_stadiamap(sabattus_large, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$chlorophyll_a_RFU, na.rm = T)-min(sab26Aug21$chlorophyll_a_RFU, na.rm = T))/2)+min(sab26Aug21$chlorophyll_a_RFU, na.rm = T),
                        name="RFU") +
  annotate("point", x = -70.09111, y = 44.17221, colour = "blue", size = 3) + # Hooper Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle(substitute(paste("Chl. ", italic("a "), "(RFU)"))) +
  xlab("Longitude") +
  ylab("Latitude")

aub_chl <- get_stadiamap(auburn, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
             data = aub30Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$chlorophyll_a_RFU, na.rm = T)-min(aub30Aug21$chlorophyll_a_RFU, na.rm = T))/2)+min(aub30Aug21$chlorophyll_a_RFU, na.rm = T),
                        name="RFU") +
  annotate("point", x = -70.24266, y = 44.15925, colour = "blue", size = 3) + # Townsend Brook
  theme(legend.position = "bottom", text = element_text(size=14, color = "black"),
        legend.key.width = unit(6, "mm"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.38, 0.45),
        axis.text.x = element_text(angle = 45, vjust = 0.5), plot.tag = element_text(color = "white"))+
  ggtitle(substitute(paste("Chl. ", italic("a "), "(RFU)"))) +
  xlab("Longitude") +
  ylab("Latitude")

fig3 <- ggarrange(sab_temp, sab_spc, sab_turb, sab_chl, 
                  aub_temp, aub_spc, aub_turb, aub_chl,
                  ncol=4, nrow = 2, align = "h",
                  labels = c("a","b","c", "d", "e", "f", "g", "h")) +
  theme(plot.margin = margin(.1, .1, .1, .1, "cm"))

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/Fig3.png", fig3,
       height = 9, width = 12)

# Figure 4

aub_do <- get_stadiamap(auburn, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
             data = aub30Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$oxygenDissolved_mgl, na.rm = T)-min(aub30Aug21$oxygenDissolved_mgl, na.rm = T))/2)+min(aub30Aug21$oxygenDissolved_mgl, na.rm = T),
                        name="mg/L") +
  annotate("point", x = -70.24266, y = 44.15925, colour = "blue", size = 3) + # Townsend Brook
  theme(legend.position = "none", text = element_text(size=12, color = "black"),
        legend.key.width = unit(6, "mm"), 
        plot.title = element_text(size = 12, color = "black"), axis.text=element_text(colour="black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.90, 0.90)) +
  geom_segment(aes(x =-70.23900, y = 44.15262, xend = -70.24000, yend = 44.15264),
               arrow = arrow(length = unit(0.1,"cm"), angle = 30), color = "white", size = 0.5) +
  scale_x_continuous(labels = label_number(accuracy = 0.01), expand = c(0,0), 
                     breaks = seq(min(aub30Aug21$longitude_gps_deg), max(aub30Aug21$longitude_gps_deg), by = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = c(0,0), 
                     breaks = seq(min(aub30Aug21$latitude_gps_deg), max(aub30Aug21$latitude_gps_deg), by = 0.01)) +
  xlab("Longitude") +
  ylab("Latitude")

aub_do_time <- ggplot() + 
  geom_point(data = aub30Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, 
                                    color = oxygenDissolved_mgl), size = 2)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$oxygenDissolved_mgl, na.rm = T)-min(aub30Aug21$oxygenDissolved_mgl, na.rm = T))/2)+min(aub30Aug21$oxygenDissolved_mgl, na.rm = T),
                        name="(mg/L)") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(aub30Aug21, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), 
            color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("DO (mg/L)")

aub_pH <- get_stadiamap(auburn, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
             data = aub30Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$pH, na.rm = T)-min(aub30Aug21$pH, na.rm = T))/2)+min(aub30Aug21$pH, na.rm = T),
                        name="") +
  annotate("point", x = -70.24266, y = 44.15925, colour = "blue", size = 3) + # Townsend Brook
  theme(legend.position = "none", text = element_text(size=12, color = "black"),
        legend.key.width = unit(6, "mm"), 
        plot.title = element_text(size = 12, color = "black"), axis.text=element_text(colour="black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.90, 0.90))+
  geom_segment(aes(x =-70.23900, y = 44.15262, xend = -70.24000, yend = 44.15264),
               arrow = arrow(length = unit(0.1,"cm"), angle = 30), color = "white", size = 0.5) +
  scale_x_continuous(labels = label_number(accuracy = 0.01), expand = c(0,0), 
                     breaks = seq(min(aub30Aug21$longitude_gps_deg), max(aub30Aug21$longitude_gps_deg), by = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = c(0,0), 
                     breaks = seq(min(aub30Aug21$latitude_gps_deg), max(aub30Aug21$latitude_gps_deg), by = 0.01)) +
  xlab("Longitude") +
  ylab("Latitude")

aub_pH_time <- ggplot() + 
  geom_point(data = aub30Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 2)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$pH, na.rm = T)-min(aub30Aug21$pH, na.rm = T))/2)+min(aub30Aug21$pH, na.rm = T),
                        name="pH") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(aub30Aug21, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("pH")

aub_spc <- get_stadiamap(auburn, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = aub30Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$specificConductance_uscm, na.rm = T)-min(aub30Aug21$specificConductance_uscm, na.rm = T))/2)+min(aub30Aug21$specificConductance_uscm, na.rm = T),
                        name="uS/cm") +
  annotate("point", x = -70.24266, y = 44.15925, colour = "blue", size = 3) + # Townsend Brook
  theme(legend.position = "none", text = element_text(size=12, color = "black"),
        legend.key.width = unit(6, "mm"),
        plot.title = element_text(size = 12, color = "black"), axis.text=element_text(colour="black"),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.90, 0.90))+
  geom_segment(aes(x =-70.23900, y = 44.15262, xend = -70.24000, yend = 44.15264),
               arrow = arrow(length = unit(0.1,"cm"), angle = 30), color = "white", size = 0.5) +
  scale_x_continuous(labels = label_number(accuracy = 0.01), expand = c(0,0), 
                     breaks = seq(min(aub30Aug21$longitude_gps_deg), max(aub30Aug21$longitude_gps_deg), by = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = c(0,0), 
                     breaks = seq(min(aub30Aug21$latitude_gps_deg), max(aub30Aug21$latitude_gps_deg), by = 0.01)) +
  xlab("Longitude") +
  ylab("Latitude")

aub_spc_time <- ggplot() + 
  geom_point(data = aub30Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), 
                                    y = specificConductance_uscm, color = specificConductance_uscm), size = 2)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(aub30Aug21$specificConductance_uscm, na.rm = T)-min(aub30Aug21$specificConductance_uscm, na.rm = T))/2)+min(aub30Aug21$specificConductance_uscm, na.rm = T),
                        name="Spec. Cond.") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(aub30Aug21, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab(paste0("Sp.C. (", "\u00B5", "S/cm)"))

fig4 <- ggarrange(aub_do_time, aub_do, aub_pH_time, aub_pH, aub_spc_time, aub_spc,
                  ncol=2, nrow = 3, align = "h", labels = c("a", "", "b", "", "c")) #Figure 4

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/Fig4.png", fig4,
       height = 9, width = 6)

# Auburn hot spot parameter correlations (to support Figure 4)

aub30Aug21_hotspot <-  aub30Aug21 %>% 
  filter(timestamp_gps_sec > 1630340529 - 60 & timestamp_gps_sec < 1630340529 + 60) # 60 seconds on either end of the peak

aub30Aug21.anoms <- aub30Aug21_hotspot %>% 
  select(temperatureWater_degC, pH, specificConductance_uscm, chlorophyll_a_RFU, 
         oxygenDissolved_mgl)

chart.Correlation(aub30Aug21.anoms, histogram = TRUE, pch = 19, method = "spearman")

# Figure 5

Jun11 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun11Jun21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_blank(),
        axis.text.y = element_text(color = c("black", "black", "black", "black", "black", "white"), hjust = 0, size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "June 11, 2021") +
  labs(tag=expression(bold("a")), hjust = -0.1)

Jun15 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun15Jun21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(color = "white", size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "June 15, 2021") +
  labs(tag=expression(bold("b")), hjust = -0.1)

Jun21 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun21Jun21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(color = "white", size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "June 21, 2021") +
  labs(tag=expression(bold("c")), hjust = -0.1)

Jul01 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun01Jul21.HC) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(color = c("black", "black", "black", "black", "black", "white"), hjust = 0, size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "July 1, 2021") +
  labs(tag=expression(bold("d")), hjust = -0.1)

Jul22 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun22Jul21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(color = "white", size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "July 22, 2021") +
  labs(tag=expression(bold("e")), hjust = -0.1)

Aug05 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun05Aug21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(color = "white", size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "August 5, 2021") +
  labs(tag=expression(bold("f")), hjust = -0.1)

Aug10 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun10Aug21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 90, vjust = 0, color = "black", size = 12),
        axis.ticks = element_blank(), axis.text.y = element_text(color = c("black", "black", "black", "black", "black", "white"), hjust = 0, size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "August 10, 2021") +
  labs(tag=expression(bold("g")), hjust = -0.1)

Aug27 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun27Aug21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 90, vjust = 0, color = "black", size = 12),
        axis.ticks = element_blank(), axis.text.y = element_text(color = "white", size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "August 27, 2021") +
  labs(tag=expression(bold("h")), hjust = -0.1)

Sep16 <- get_stadiamap(sunapee, zoom = 16, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sun16Sep21) +
  annotate("point", x = -72.03080, y = 43.40967, colour = "blue", size = 5) + # Herrick Cove South
  annotate("point", x = -72.03589, y = 43.41281, colour = "blue", size = 5) + # Herrick Cove North
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sunapee_all$specificConductance_uscm, na.rm = T)-min(sunapee_all$specificConductance_uscm, na.rm = T))/2)+min(sunapee_all$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)"), limits = c(103.8537, 112.3589)) +
  theme(legend.position = "bottom", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(1, 1, 1, 1), plot.tag.position = c(0.1, 0.9),
        legend.title = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 90, vjust = 0, color = "black", size = 12),
        axis.ticks = element_blank(), axis.text.y = element_text(color = "white", size = 12),
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 14, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  ggtitle(label = "September 16, 2021") +
  labs(tag=expression(bold("i")), hjust = -0.1)

fig5 <- ggarrange(Jun11, Jun15, Jun21, Jul01, Jul22, Aug05, Aug10, Aug27, Sep16,
                  ncol=3, nrow = 3, align = "h", common.legend = T, legend="bottom") #Figure 5
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/Fig5.png", fig5,
       height = 10, width = 12)

#######################################################################
# SUPPORTING INFORMATION
#######################################################################

# Calculate total deployment time for each path for Supp Table 2
aub30Aug21.time <- difftime(as.POSIXct(aub30Aug21$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(aub30Aug21$timestamp_gps_sec[nrow(aub30Aug21)], origin = "1970-01-01"), units="hours")
difftime(as.POSIXct(sab20Aug21$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sab20Aug21$timestamp_gps_sec[nrow(sab20Aug21)], origin = "1970-01-01"), units="hours")
difftime(as.POSIXct(sab26Aug21$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sab26Aug21$timestamp_gps_sec[nrow(sab26Aug21)], origin = "1970-01-01"), units="hours")
chn06Aug21.time <- difftime(as.POSIXct(chn06Aug21$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(chn06Aug21$timestamp_gps_sec[nrow(chn06Aug21)], origin = "1970-01-01"), units="hours")
difftime(as.POSIXct(chn05Oct21$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(chn05Oct21$timestamp_gps_sec[nrow(chn05Oct21)], origin = "1970-01-01"), units="hours")
difftime(as.POSIXct(chn12Oct21$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(chn12Oct21$timestamp_gps_sec[nrow(chn12Oct21)], origin = "1970-01-01"), units="hours")
difftime(as.POSIXct(chn21Oct21$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(chn21Oct21$timestamp_gps_sec[nrow(chn21Oct21)], origin = "1970-01-01"), units="hours")

abs(difftime(as.POSIXct(sun11Jun21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun11Jun21.HC$timestamp_gps_sec[nrow(sun11Jun21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun11Jun21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun11Jun21.NW$timestamp_gps_sec[nrow(sun11Jun21.NW)], origin = "1970-01-01"), units="hours"))
abs(difftime(as.POSIXct(sun15Jun21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun15Jun21.HC$timestamp_gps_sec[nrow(sun15Jun21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun15Jun21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun15Jun21.NW$timestamp_gps_sec[nrow(sun15Jun21.NW)], origin = "1970-01-01"), units="hours"))
abs(difftime(as.POSIXct(sun21Jun21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun21Jun21.HC$timestamp_gps_sec[nrow(sun21Jun21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun21Jun21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun21Jun21.NW$timestamp_gps_sec[nrow(sun21Jun21.NW)], origin = "1970-01-01"), units="hours"))
difftime(as.POSIXct(sun01Jul21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun01Jul21.HC$timestamp_gps_sec[nrow(sun01Jul21.HC)], origin = "1970-01-01"), units="hours")
abs(difftime(as.POSIXct(sun22Jul21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun22Jul21.HC$timestamp_gps_sec[nrow(sun22Jul21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun22Jul21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun22Jul21.NW$timestamp_gps_sec[nrow(sun22Jul21.NW)], origin = "1970-01-01"), units="hours"))
abs(difftime(as.POSIXct(sun05Aug21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun05Aug21.HC$timestamp_gps_sec[nrow(sun05Aug21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun05Aug21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun05Aug21.NW$timestamp_gps_sec[nrow(sun05Aug21.NW)], origin = "1970-01-01"), units="hours"))
abs(difftime(as.POSIXct(sun10Aug21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun10Aug21.HC$timestamp_gps_sec[nrow(sun10Aug21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun10Aug21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun10Aug21.NW$timestamp_gps_sec[nrow(sun10Aug21.NW)], origin = "1970-01-01"), units="hours"))
abs(difftime(as.POSIXct(sun27Aug21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun27Aug21.HC$timestamp_gps_sec[nrow(sun27Aug21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun27Aug21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun27Aug21.NW$timestamp_gps_sec[nrow(sun27Aug21.NW)], origin = "1970-01-01"), units="hours"))
abs(difftime(as.POSIXct(sun16Sep21.HC$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun16Sep21.HC$timestamp_gps_sec[nrow(sun16Sep21.HC)], origin = "1970-01-01"), units="hours")+
      difftime(as.POSIXct(sun16Sep21.NW$timestamp_gps_sec[1], origin = "1970-01-01"), as.POSIXct(sun16Sep21.NW$timestamp_gps_sec[nrow(sun16Sep21.NW)], origin = "1970-01-01"), units="hours"))

# Supp Fig 2

china_suppfig2 <- c(left = -69.6070, bottom = 44.4420, right = -69.6010, top = 44.4470)

chn_map <- get_stadiamap(china_suppfig2, zoom = 18, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
             data = chn05Oct21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(chn05Oct21$pH, na.rm = T)-min(chn05Oct21$pH, na.rm = T))/2)+min(chn05Oct21$pH, na.rm = T),
                        name="pH") +
  theme(text = element_text(size = 16, color = "black"),
        legend.key.width = unit(6, "mm"), 
        legend.title = element_text(size=16), legend.position = "bottom",
        legend.text = element_text(angle = 0, hjust = 0, vjust = 0, size = 16), axis.text=element_text(colour="black", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.25, 0.92),
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  xlab("Longitude") +
  ylab("Latitude")

chn05Oct21_labeled_plot <- travel_effects(chn05Oct21, preset_loiter_time) %>% select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU", "longitude_gps_deg", "latitude_gps_deg")

chn05Oct21_labeled_plot$loiter_tag <- factor(chn05Oct21_labeled_plot$loiter_tag, levels=c("pre-loiter", "loiter", "post-loiter"))

chn_loitermap <- get_stadiamap(china_suppfig2, zoom = 18, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = loiter_tag), size=1, shape = 10,
             data = chn05Oct21_labeled_plot) +
  scale_color_manual(values = c("blue", "yellow", "red"), name = "Tag", labels = c("Pre-Loiter", "Loiter", "Post-Loiter")) +
  theme(text = element_text(size = 16, color = "black"), legend.position = "bottom",
        legend.key.width = unit(6, "mm"), legend.key = element_rect(fill = NA),
        legend.title = element_text(size= 16), 
        legend.text = element_text(size = 16), axis.text=element_text(colour="black", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.25, 0.92),
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  guides(color = guide_legend(override.aes = list(size = 5, shape = 16))) +
  xlab("Longitude") +
  ylab("Latitude")

chn_time <- ggplot() + 
  geom_point(data = chn05Oct21, aes(x = as.POSIXct(timestamp_gps_sec, origin="1970-01-01", tz="EST"), y = pH, 
                                    color=velocity_gps_mps), size = 3)+theme_bw()+xlab("Timestamp") + 
  theme_bw()+
  xlab("Timestamp") + 
  theme(legend.position = "bottom", axis.title=element_text(size = 16, color = "black"), text = element_text(size = 16, color = "black"),
        legend.key.width = unit(10, "mm"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_text(size=16), axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(angle = 0, hjust = 0, vjust = 0, size = 16), panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 2),plot.tag.position = c(0.08, 0.9))+
  scale_x_datetime(breaks = scales::date_breaks("5 mins"), date_labels = "%H:%M")+
  scale_color_gradient(low = "grey80", high = "black", breaks=c(0, 0.5, 1),
                       limits=c(0,1.5), guide = guide_colorbar(title = "m/s",
                                                               title.vjust = 0.5)) +
  xlab("Time (EDT)")+
  ylab("pH")

suppfig2 <- ggarrange(ggarrange(chn_loitermap, chn_map, ncol = 2, labels = c("a", "b")),
                  chn_time, align = "v",
                  nrow = 2, heights = c(3, 2),
                  labels = c("", "c")) #Supp Fig 2

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig2.png", suppfig2,
       height = 10, width = 12)

# Supp Fig 3

sab26Aug21_labeled <- travel_effects(sab26Aug21, preset_loiter_time) %>% 
select("lake", "date", "timestamp_gps_sec", "loiter_tag", "loiterID", "loiter_num", "velocity_gps_mps", "temperatureWater_degC", "pH", "chlorophyll_a_RFU", "specificConductance_uscm", "oxygenDissolved_mgl", "turbidity_NTU")

sab26Aug21.labeled.chl <- sab26Aug21_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number"))%>%
  mutate(date = "2021-08-26", path = "perimeter", lake = "Sabattus")

# add loiter 6.5a and 6.5b, which were not technically loiters but suffice for this purpose, just stopped in place for 2 minutes while dragging ASV with our boat

sab26Aug21_loiter6.5a <- sab26Aug21 %>% # 44.16268, -70.09155
  dplyr::filter(timestamp_gps_sec > 1629995759 & timestamp_gps_sec < 1629995882) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU))

sab26Aug21_loiter6.5b <- sab26Aug21 %>% # 44.1679, -70.08805
  dplyr::filter(timestamp_gps_sec > 1629996779 & timestamp_gps_sec < 1629996902) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU))

aub30Aug21.labeled.chl <- aub30Aug21_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number"))%>%
  mutate(date = "2021-08-30", path = "Townsend", lake = "Auburn")

chn05Oct21.labeled.chl <- chn05Oct21_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>% 
  mutate(date = "2021-10-05", path = "capstone", lake = "China")

chn12Oct21.labeled.chl <- chn12Oct21_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>% 
  mutate(date = "2021-10-12", path = "capstone", lake = "China")

chn21Oct21.labeled.chl <- chn21Oct21_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>% 
  mutate(date = "2021-10-21", path = "capstone", lake = "China")

sun11Jun21.HC.labeled.chl <- sun11Jun21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>% 
  mutate(date = "2021-06-11", path = "HC", lake = "Sunapee")

sun11Jun21.NW.labeled.chl <- sun11Jun21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-06-11", path = "NW", lake = "Sunapee")

sun15Jun21.HC.labeled.chl <- sun15Jun21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-06-15", path = "HC", lake = "Sunapee")

sun15Jun21.NW.labeled.chl <- sun15Jun21.NW_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-06-15", path = "NW", lake = "Sunapee")

sun21Jun21.HC.labeled.chl <- sun21Jun21.HC_labeled %>%
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-06-21", path = "HC", lake = "Sunapee")

sun21Jun21.NW.labeled.chl <- sun21Jun21.NW_labeled %>%   
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-06-21", path = "NW", lake = "Sunapee")

sun01Jul21.HC.labeled.chl <- sun01Jul21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-07-01", path = "HC", lake = "Sunapee")

sun22Jul21.HC.labeled.chl <- sun22Jul21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-07-22", path = "HC", lake = "Sunapee")

sun22Jul21.NW.labeled.chl <- sun22Jul21.NW_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-07-22", path = "NW",lake = "Sunapee")

sun05Aug21.HC.labeled.chl <- sun05Aug21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-08-05", path = "HC",lake = "Sunapee")

sun05Aug21.NW.labeled.chl <- sun05Aug21.NW_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-08-05", path = "NW", lake = "Sunapee")

sun10Aug21.HC.labeled.chl <- sun10Aug21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-08-10", path = "HC", lake = "Sunapee")

sun10Aug21.NW.labeled.chl <- sun10Aug21.NW_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-08-10", path = "NW", lake = "Sunapee")

sun27Aug21.HC.labeled.chl <- sun27Aug21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-08-27", path = "HC", lake = "Sunapee")

sun27Aug21.NW.labeled.chl <- sun27Aug21.NW_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-08-27", path = "NW", lake = "Sunapee")

sun16Sep21.HC.labeled.chl <- sun16Sep21.HC_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-09-16", path = "HC", lake = "Sunapee")

sun16Sep21.NW.labeled.chl <- sun16Sep21.NW_labeled %>% 
  dplyr::group_by(loiterID) %>% 
  dplyr::summarize(median_sondechl_120 = median(chlorophyll_a_RFU)) %>% 
  separate(loiterID, c("flag", "number"), "_") %>% 
  filter(flag == "loiter") %>% 
  unite(loiterID, c("flag", "number")) %>%
  mutate(date = "2021-09-16", path = "NW", lake = "Sunapee")

alldates.sondechl <- rbind(sun11Jun21.HC.labeled.chl, sun11Jun21.NW.labeled.chl, sun15Jun21.HC.labeled.chl, 
                           sun15Jun21.NW.labeled.chl, sun21Jun21.HC.labeled.chl, sun21Jun21.NW.labeled.chl,
                           sun01Jul21.HC.labeled.chl, sun22Jul21.HC.labeled.chl, sun22Jul21.NW.labeled.chl, 
                           sun05Aug21.HC.labeled.chl, sun05Aug21.NW.labeled.chl, sun10Aug21.HC.labeled.chl, 
                           sun10Aug21.NW.labeled.chl, sun27Aug21.HC.labeled.chl, sun27Aug21.NW.labeled.chl, 
                           sun16Sep21.HC.labeled.chl, sun16Sep21.NW.labeled.chl, aub30Aug21.labeled.chl, 
                           sab26Aug21.labeled.chl, chn05Oct21.labeled.chl,
                           chn12Oct21.labeled.chl, chn21Oct21.labeled.chl) %>% 
  mutate(daymonth = format(as.Date(date), "%m/%d/%Y")) %>% 
  mutate(row_ID = paste(lake, date, path, loiterID, sep = '_'))

# Read in lab extracted chl-a data
labchl <- read.csv("~/Dropbox/Bates/Manuscripts/ASVLimno_MS/data/chlorophyll/ASV_insitu_AUB_SAB_SUN.csv", header = T)

labchl_sab26Aug21_loiter6.5a <- labchl %>% 
  dplyr::filter(loiterID == "loiter_6.5a") %>% 
  select(conc_ugL_lake)

labchl_sab26Aug21_loiter6.5b <- labchl %>% 
  dplyr::filter(loiterID == "loiter_6.5b") %>% 
  select(conc_ugL_lake)

labchl$date <- as.Date(labchl$date)

labchl <- labchl %>% 
  filter(field_rep == "A") %>% 
  mutate(row_ID = paste(lake, date, path, loiterID, sep = '_'))

# Join sonde and lab chl-a dataframes
chl.regression <- merge(alldates.sondechl, labchl, by = "row_ID")
chl.regression.aub <- chl.regression[chl.regression$lake.x == "Auburn",]
chl.regression.sab <- chl.regression[chl.regression$lake.x == "Sabattus",]
chl.regression.chn <- chl.regression[chl.regression$lake.x == "China",]
chl.regression.sun <- chl.regression[chl.regression$lake.x == "Sunapee",]

# Linear regressions for sonde vs lab measured chl a
aub.lm <- lm(median_sondechl_120~conc_ugL_lake, data = chl.regression.aub)
sab.lm <- lm(median_sondechl_120~conc_ugL_lake, data = chl.regression.sab)
chn.lm <- lm(median_sondechl_120~conc_ugL_lake, data = chl.regression.chn)
sun.lm <- lm(median_sondechl_120~conc_ugL_lake, data = chl.regression.sun)

plot.chl.all <- ggplot() +
  geom_abline(slope = aub.lm$coefficients[2], intercept = aub.lm$coefficients[1], color = "#FF6DB6", size = 0.75, alpha = 0.5) +
  geom_abline(slope = sab.lm$coefficients[2], intercept = sab.lm$coefficients[1], color = "#22CF22", size = 0.75, alpha = 0.5) +
  geom_abline(slope = sun.lm$coefficients[2], intercept = sun.lm$coefficients[1], color = "#006DDB", size = 0.75, alpha = 0.5) +
  geom_point(data = chl.regression.aub, aes(x = conc_ugL_lake, y = median_sondechl_120, color=lake.x), size=4, alpha = 0.5)+
  geom_point(data = chl.regression.sab, aes(x = conc_ugL_lake, y = median_sondechl_120, color=lake.x), size=4, alpha = 0.5)+
  geom_point(data = chl.regression.sun, aes(x = conc_ugL_lake, y = median_sondechl_120, color=lake.x), size=4, alpha = 0.5) + #use #228833 to add china
  annotate("point", x = as.numeric(labchl_sab26Aug21_loiter6.5a), y = as.numeric(sab26Aug21_loiter6.5a), color = "#22CF22", size = 4, alpha = 0.5) +
  annotate("point", x = as.numeric(labchl_sab26Aug21_loiter6.5b), y = as.numeric(sab26Aug21_loiter6.5b), color = "#22CF22", size = 4, alpha = 0.5) +
  scale_color_manual(values = c("#FF6DB6", "#22CF22", "#006DDB"))+
  theme(legend.position = "bottom", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black"),
        legend.key = element_rect(fill = NA), text=element_text(size=14, color = "black"), legend.title = element_blank())+ 
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 2))+
  geom_rect(aes(xmin = 0, xmax = 3, ymin = 0, ymax = 0.5), fill = NA, color="black")+
  xlab(expression("Lab-Analyzed Chlorophyll"~italic("a")~"(\u00B5g/L)"))+
  ylab(expression("Sonde-Measured Chlorophyll"~italic("a")~"(RFU)"))

plot.chl.aub.sun <- ggplot() +  
  geom_abline(slope = aub.lm$coefficients[2], intercept = aub.lm$coefficients[1], color = "#FF6DB6", size = 0.75, alpha = 0.4) +
  geom_abline(slope = sun.lm$coefficients[2], intercept = sun.lm$coefficients[1], color = "#006DDB", size = 0.75, alpha = 0.4) +
  geom_point(data = chl.regression.aub, aes(x = conc_ugL_lake, y = median_sondechl_120, color=lake.x), size=4, alpha = 0.7)+
  geom_point(data = chl.regression.sun, aes(x = conc_ugL_lake, y = median_sondechl_120, color=lake.x), size=4, alpha = 0.7) +
  scale_color_manual(values = c("#FF6DB6", "#006DDB"))+
  theme(legend.position = "bottom", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black"),
        legend.key = element_rect(fill = NA), text=element_text(size=14, color = "black"), legend.title = element_blank())+ 
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 0.5))+
  xlab("")+
  ylab("")

suppfig3 <- grid_arrange_shared_legend(plot.chl.all, plot.chl.aub.sun, ncol=2, widths = c(1.4, 0.6), position = "bottom")

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig3.png", height = 8, width = 10, suppfig3)

# Supp Fig 4

df.list <- list(aub30Aug21, sab20Aug21, sab26Aug21, chn06Aug21, chn28Sep21, chn05Oct21, chn12Oct21, chn21Oct21,
                sun11Jun21, sun15Jun21, sun21Jun21, sun01Jul21.HC, sun22Jul21, sun05Aug21, sun10Aug21, sun27Aug21,
                sun16Sep21)

names(df.list[[10]])
nms <- names(df.list[[1]])

namelist <- list()

for (ii in 1:length(df.list)){
  nms <- colnames(df.list[[ii]])
  namelist[[paste0("df",ii)]] <- nms
}

namelist
comnames <- namelist %>% reduce(intersect)
comnames

df <- df.list %>% reduce(full_join,by=comnames) %>% 
  dplyr::select(date,lake,pH,chlorophyll_a_RFU,turbidity_NTU,specificConductance_uscm,
                temperatureWater_degC,oxygenDissolved_mgl)

levels(as.factor(df$date))
df$date <- as.Date(df$date,format="%Y-%m-%d")  #assign "Date" format

chl <- ggplot(data = df, aes(date, chlorophyll_a_RFU)) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 16, angle = 315, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 18,color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_boxplot(aes(color = lake, group = date), width = 8, show.legend = FALSE) +
  scale_color_manual(values=c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_x_date(breaks = "1 month", date_labels = "%b")+
  scale_y_continuous(lim = c(0, 3), breaks = seq(0, 3, 0.5))+
  labs(x = "", y = expression(paste("Chl. ", italic(" a")," (RFU)")))

ph <- ggplot(data = df, aes(date, pH)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, angle = 315, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 18,color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_boxplot(aes(color = lake, group = date), width = 8, show.legend = FALSE) +
  scale_color_manual(values=c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_x_date(breaks = "1 month", date_labels = "%b")+
  labs(x = "", y = expression(paste("pH")))

spc <- ggplot(data = df, aes(date, specificConductance_uscm)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 16, angle = 315, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 18,color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_boxplot(aes(color = lake, group = date), width = 8, show.legend = FALSE) +  
  scale_color_manual(values=c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_x_date(breaks = "1 month", date_labels = "%b")+
  labs(x = "", y = paste0("Sp. C. (", "\u00B5", "S/cm)"))

do <- ggplot(data = df, aes(date, oxygenDissolved_mgl)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 16, angle = 315, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 18,color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_boxplot(aes(color = lake, group = date), width = 8, show.legend = FALSE) +    
  scale_color_manual(values=c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_x_date(breaks = "1 month", date_labels = "%b")+
  labs(x = "", y = expression(paste("DO (mg/L)")))

temp <- ggplot(data = df, aes(date, temperatureWater_degC, group=date)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 16, angle = 315, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 18,color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_boxplot(aes(color = lake, group = date), width = 8, show.legend = FALSE) + 
  scale_color_manual(values=c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_x_date(breaks = "1 month", date_labels = "%b")+
  labs(x = "", y = expression(paste("Temp. (ºC)")), color="Lake")

turb <- ggplot(data = df, aes(date, turbidity_NTU)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 16, angle = 315, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 18,color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_boxplot(aes(color = lake, group = date), width = 8, show.legend = TRUE) + 
  scale_color_manual(values=c("AUB" = "#FF6DB6", "CHN" = "#B66DFF", "SAB" = "#22CF22", "SUN" = "#006DDB"))+
  scale_x_date(breaks = "1 month", date_labels = "%b",
               limits = as.Date(c('2021-06-11','2021-10-21')))+
  labs(x="",y=expression(paste("Turb. (NTU)")))

suppfig4 <- ggarrange(temp, do, ph, spc, turb, chl, 
                                  nrow=3, ncol=2, common.legend = TRUE, legend= "bottom",
                      labels = c("a", "b", "c", "d", "e", "f"))

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig4.png", height = 12, width = 10, suppfig4)


# Supp Figs 5-21

# Supp Fig 5 (SABATTUS)

a <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
             data = sab20Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$pH, na.rm = T)-min(sabattus_all$pH, na.rm = T))/2)+min(sabattus_all$pH, na.rm = T),
                        name="pH") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2), plot.tag.position = c(0.9, 0.95))

aa <- ggplot() + 
  geom_point(data = sab20Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$pH,na.rm = T)-min(sabattus_all$pH, na.rm = T))/2)+min(sabattus_all$pH, na.rm = T),
                        name="pH") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab20Aug21, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("pH")

b <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
             data = sab20Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$temperatureWater_degC, na.rm = T)-min(sabattus_all$temperatureWater_degC, na.rm = T))/2)+min(sabattus_all$temperatureWater_degC, na.rm = T),
                        name="Temperature (ºC)") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

bb <- ggplot() + 
  geom_point(data = sab20Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$temperatureWater_degC, na.rm = T)-min(sabattus_all$temperatureWater_degC, na.rm = T))/2)+min(sabattus_all$temperatureWater_degC, na.rm = T),
                        name="Temperature (ºC)") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab20Aug21, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("Temperature (ºC)")

c <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_mscm), size=1, shape = 10,
             data = sab20Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab20Aug21$specificConductance_mscm, na.rm = T)-min(sab20Aug21$specificConductance_mscm, na.rm = T))/2)+min(sab20Aug21$specificConductance_mscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

cc <- ggplot() + 
  geom_point(data = sab20Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_mscm, color = specificConductance_mscm), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab20Aug21$specificConductance_mscm, na.rm = T)-min(sab20Aug21$specificConductance_mscm, na.rm = T))/2)+min(sab20Aug21$specificConductance_mscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab20Aug21, specificConductance_mscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab(paste0("Specific conductance (", "\u00B5", "S/cm)"))

d <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_ugl), size=1, shape = 10,
             data = sab20Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab20Aug21$chlorophyll_a_ugl, na.rm = T)-min(sab20Aug21$chlorophyll_a_ugl, na.rm = T))/2)+min(sab20Aug21$chlorophyll_a_ugl, na.rm = T),
                        name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

dd <- ggplot() + 
  geom_point(data = sab20Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_ugl, color = chlorophyll_a_ugl), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab20Aug21$chlorophyll_a_ugl, na.rm = T)-min(sab20Aug21$chlorophyll_a_ugl, na.rm = T))/2)+min(sab20Aug21$chlorophyll_a_ugl, na.rm = T),
                        name=expression(paste("Chlorophyll ", italic("a"), " (", "\u00B5", "g/L)"))) +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab20Aug21, chlorophyll_a_ugl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab(expression(paste("Chlorophyll ", italic("a"), " (", "\u00B5", "g/L)")))

e <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
             data = sab20Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$oxygenDissolved_mgl, na.rm = T)-min(sabattus_all$oxygenDissolved_mgl, na.rm = T))/2)+min(sabattus_all$oxygenDissolved_mgl, na.rm = T),
                        name="DO (mg/L)") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

ee <- ggplot() + 
  geom_point(data = sab20Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$oxygenDissolved_mgl, na.rm = T)-min(sabattus_all$oxygenDissolved_mgl,na.rm = T))/2)+min(sabattus_all$oxygenDissolved_mgl, na.rm = T),
                        name="Dissolved oxygen (mg/L)") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab20Aug21, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("Dissolved oxygen (mg/L)")

f <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
             data = sab20Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$turbidity_NTU, na.rm = T)-min(sabattus_all$turbidity_NTU, na.rm = T))/2)+min(sabattus_all$turbidity_NTU, na.rm = T),
                        name="Turb (NTU)") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

ff <- ggplot() + 
  geom_point(data = sab20Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$turbidity_NTU, na.rm = T)-min(sabattus_all$turbidity_NTU, na.rm = T))/2)+min(sabattus_all$turbidity_NTU, na.rm = T),
                        name="Turbidity (NTU)") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab20Aug21, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("Turbidity (NTU)")

suppfig5 <- ggarrange(aa, bb, cc, dd, ee, ff,
                      a, b, c, d, e, f,
                      ncol=6, nrow = 2,
                      labels = c("a", "b", "c", "d", "e", "f"))
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig5.png", suppfig5,
       height = 8.5, width = 13)

# Supp Fig 6 (SABATTUS)

a <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$pH, na.rm = T)-min(sabattus_all$pH, na.rm = T))/2)+min(sabattus_all$pH, na.rm = T),
                        name="pH") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

aa <- ggplot() + 
  geom_point(data = sab26Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$pH,na.rm = T)-min(sabattus_all$pH, na.rm = T))/2)+min(sabattus_all$pH, na.rm = T),
                        name="pH") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab26Aug21, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("pH")

b <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$temperatureWater_degC, na.rm = T)-min(sabattus_all$temperatureWater_degC, na.rm = T))/2)+min(sabattus_all$temperatureWater_degC, na.rm = T),
                        name="Temperature (ºC)") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

bb <- ggplot() + 
  geom_point(data = sab26Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$temperatureWater_degC, na.rm = T)-min(sabattus_all$temperatureWater_degC, na.rm = T))/2)+min(sabattus_all$temperatureWater_degC, na.rm = T),
                        name="Temperature (ºC)") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab26Aug21, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("Temperature (ºC)")

c <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$specificConductance_uscm, na.rm = T)-min(sab26Aug21$specificConductance_uscm, na.rm = T))/2)+min(sab26Aug21$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

cc <- ggplot() + 
  geom_point(data = sab26Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$specificConductance_uscm, na.rm = T)-min(sab26Aug21$specificConductance_uscm, na.rm = T))/2)+min(sab26Aug21$specificConductance_uscm, na.rm = T),
                        name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab26Aug21, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab(paste0("Specific conductance (", "\u00B5", "S/cm)"))

d <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$chlorophyll_a_RFU, na.rm = T)-min(sab26Aug21$chlorophyll_a_RFU, na.rm = T))/2)+min(sab26Aug21$chlorophyll_a_RFU, na.rm = T),
                        name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

dd <- ggplot() + 
  geom_point(data = sab26Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sab26Aug21$chlorophyll_a_RFU, na.rm = T)-min(sab26Aug21$chlorophyll_a_RFU, na.rm = T))/2)+min(sab26Aug21$chlorophyll_a_RFU, na.rm = T),
                        name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab26Aug21, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab(expression(paste("Chlorophyll ", italic("a"), " (RFU)")))

e <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$oxygenDissolved_mgl, na.rm = T)-min(sabattus_all$oxygenDissolved_mgl, na.rm = T))/2)+min(sabattus_all$oxygenDissolved_mgl, na.rm = T),
                        name="DO (mg/L)") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

ee <- ggplot() + 
  geom_point(data = sab26Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$oxygenDissolved_mgl, na.rm = T)-min(sabattus_all$oxygenDissolved_mgl,na.rm = T))/2)+min(sabattus_all$oxygenDissolved_mgl, na.rm = T),
                        name="Dissolved oxygen (mg/L)") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab26Aug21, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("Dissolved oxygen (mg/L)")

f <- get_stadiamap(sabattus, zoom = 14, maptype = "stamen_terrain") %>% ggmap()+
  geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
             data = sab26Aug21) +
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$turbidity_NTU, na.rm = T)-min(sabattus_all$turbidity_NTU, na.rm = T))/2)+min(sabattus_all$turbidity_NTU, na.rm = T),
                        name="Turb (NTU)") +
  theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
        legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, color = "black"),
        legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(colour="black"),
        panel.border = element_rect(color = "black", fill = NA, size = 2))

ff <- ggplot() + 
  geom_point(data = sab26Aug21, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
  scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                        midpoint = ((max(sabattus_all$turbidity_NTU, na.rm = T)-min(sabattus_all$turbidity_NTU, na.rm = T))/2)+min(sabattus_all$turbidity_NTU, na.rm = T),
                        name="Turbidity (NTU)") +
  theme_bw()+
  xlab("Timestamp") + 
  geom_line(data = asv_window(sab26Aug21, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(), text=element_text(size=12, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 2))+
  scale_x_datetime(breaks = scales::date_breaks("1 hour"), date_labels = "%H:%M")+
  xlab("Time (EDT)")+
  ylab("Turbidity (NTU)")

suppfig6 <- ggarrange(aa, bb, cc, dd, ee, ff,
                      a, b, c, d, e, f,
                      ncol=6, nrow = 2,
                      labels = c("a", "b", "c", "d", "e", "f"))

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig6.png", suppfig6,
       height = 8.5, width = 13)

# Supp Fig 7 (AUBURN)

asv_plotter_auburn <- function(map_bounds, deployment, lake_all, zoom){ 
  #where "deployment" is the dataframe of that particular deployment date
  #and "lake_all" is a combined dataframe that contains all deployments for a single lake
  #and "map_bounds" is a set boundary of the spatial extent to be mapped
  #and "zoom" looks best at 12, 14, or 16, depending on spatial scale of deployment
  
  a <- get_stadiamap(bbox = map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  aa <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("45 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  b <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  bb <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("45 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  c <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  cc <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("45 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Specific conductance (", "\u00B5", "S/cm)"))
  
  d <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  dd <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("45 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chlorophyll ", italic("a"), " (RFU)")))
  
  e <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ee <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="Dissolved oxygen (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("45 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Dissolved oxygen (mg/L)")
  
  f <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turb (NTU)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ff <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turbidity (NTU)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("45 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Turbidity (NTU)")
  
  plot <- ggarrange(aa, bb, cc, dd, ee, ff,
                    a, b, c, d, e, f,
                    ncol=6, nrow = 2,
                    labels = c("a", "b", "c", "d", "e", "f"))
  
  return(plot)
}

suppfig7 <- asv_plotter_auburn(auburn, aub30Aug21, auburn_all, 16)

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig7.png", suppfig7,
       height = 7, width = 13)

# Supp Figs 8-12 (CHINA)

asv_plotter_china_allParam <- function(map_bounds, deployment, lake_all, zoom){ 
  #where "deployment" is the dataframe of that particular deployment date
  #and "lake_all" is a combined dataframe that contains all deployments for a single lake
  #and "map_bounds" is a set boundary of the spatial extent to be mapped
  #and "zoom" looks best at 12, 14, or 16, depending on spatial scale of deployment
  
  a <- get_stadiamap(bbox = map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  aa <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  b <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  bb <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  c <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  cc <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Specific conductance (", "\u00B5", "S/cm)"))
  
  d <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  dd <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chlorophyll ", italic("a"), " (RFU)")))
  
  e <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ee <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="Dissolved oxygen (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Dissolved oxygen (mg/L)")
  
  f <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turb (NTU)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ff <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turbidity (NTU)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Turbidity (NTU)")
  
  plot <- ggarrange(aa, bb, cc, dd, ee, ff,
                    a, b, c, d, e, f,
                    ncol=6, nrow = 2,
                    labels = c("a", "b", "c", "d", "e", "f"))
  
  return(plot)
}

suppfig8 <- asv_plotter_china_allParam(china, chn06Aug21, china_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig8.png", suppfig8,
       height = 7, width = 14)

asv_plotter_china_noSpC <- function(map_bounds, deployment, lake_all, zoom){ 
  #where "deployment" is the dataframe of that particular deployment date
  #and "lake_all" is a combined dataframe that contains all deployments for a single lake
  #and "map_bounds" is a set boundary of the spatial extent to be mapped
  #and "zoom" looks best at 12, 14, or 16, depending on spatial scale of deployment
  
  a <- get_stadiamap(bbox = map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  aa <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  b <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  bb <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  d <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  dd <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chlorophyll ", italic("a"), " (RFU)")))
  
  e <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ee <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="Dissolved oxygen (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Dissolved oxygen (mg/L)")
  
  f <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turb (NTU)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y.left = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ff <- ggplot() + 
    geom_point(data = deployment, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turbidity (NTU)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(deployment, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          plot.tag.position = c(0.95, 0.95), panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("25 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Turbidity (NTU)")
  
  plot <- ggarrange(aa, bb, dd, ee, ff,
                    a, b, d, e, f,
                    ncol=5, nrow = 2,
                    labels = c("a", "b", "c", "d", "e"))
  
  return(plot)
}

suppfig9 <- asv_plotter_china_noSpC(china, chn28Sep21, china_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig9.png", suppfig9,
       height = 7, width = 14)

suppfig10 <- asv_plotter_china_noSpC(china, chn05Oct21, china_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig10.png", suppfig10,
       height = 7, width = 14)

suppfig11 <- asv_plotter_china_allParam(china, chn12Oct21, china_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig11.png", suppfig11,
       height = 7, width = 14)

suppfig12 <- asv_plotter_china_allParam(china, chn21Oct21, china_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig12.png", suppfig12,
       height = 7, width = 14)

# Supp Figs 13-21 (SUNAPEE)

asv_plotter_sunapee_preTurb <- function(map_bounds, deployment, path1, path2, lake_all, zoom){ 
  #early dates before inclusion of turbidity sensor as part of data stream
  #where deployment is the dataframe of that particular deployment date
  #and lake_all is a combined dataframe that contains all deployments for a single lake
  #and map_bounds is a set boundary of the spatial extent to be mapped
  # and zoom is either 12 for SAB or 14 for AUB, CHN
  
  a <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1), axis.text.y = element_text(color = "black"),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  aa <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  aaa <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  b <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"),  axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6), axis.text.y = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  bb <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  bbb <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  c <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1), axis.text.y = element_text(color = "black"),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  cc <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Sp. C. (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Sp. C. (", "\u00B5", "S/cm)"))
  
  ccc <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Sp. conductance (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Sp. C. (", "\u00B5", "S/cm)"))
  
  d <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6), axis.text.y = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  dd <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chl. ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chl. ", italic("a"), " (RFU)")))
  
  ddd <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chl. ", italic("a"), " (RFU)")))
  
  e <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y = element_text(color = "black"), axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ee <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("DO (mg/L)")
  
  eee <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("DO (mg/L)")
  
  plot <- ggarrange(aa, bb, cc, dd, ee,
                    aaa, bbb, ccc, ddd, eee,
                    a, b, c, d, e,
                    ncol=5, nrow = 3,
                    labels = c("a", "b", "c", "d", "e"))
  
  return(plot)
}

suppfig13 <- asv_plotter_sunapee_preTurb(sunapee, sun11Jun21, sun11Jun21.HC, sun11Jun21.NW, sunapee_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig13.png", suppfig13,
       height = 7, width = 12)

suppfig14 <- asv_plotter_sunapee_preTurb(sunapee, sun15Jun21, sun15Jun21.HC, sun15Jun21.NW, sunapee_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig14.png", suppfig14,
       height = 7, width = 12)

suppfig15 <- asv_plotter_sunapee_preTurb(sunapee, sun21Jun21, sun21Jun21.HC, sun21Jun21.NW, sunapee_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig15.png", suppfig15,
       height = 7, width = 12)

asv_plotter_sunapee_preTurb_oneSite <- function(map_bounds, deployment, path1, lake_all, zoom){ 
  #dates when deployments only happened at one of the two Sunapee sites, preTurb
  #where deployment is the dataframe of that particular deployment date
  #and lake_all is a combined dataframe that contains all deployments for a single lake
  #and map_bounds is a set boundary of the spatial extent to be mapped
  # and zoom is either 12 for SAB or 14 for AUB, CHN
  
  a <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1), axis.text.y = element_text(color = "black"),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  aa <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  b <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"),  axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6), axis.text.y = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  bb <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  c <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1), axis.text.y = element_text(color = "black"),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  cc <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Sp. C. (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Sp. C. (", "\u00B5", "S/cm)"))
  
  d <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6), axis.text.y = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  dd <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chlorophyll ", italic("a"), " (RFU)")))
  
  e <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, colour = "black"), axis.text.y = element_text(color = "black"), axis.ticks = element_blank(),
          legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ee <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("DO (mg/L)")
  
  plot <- ggarrange(aa, bb, cc, dd, ee,
                    a, b, c, d, e,
                    ncol=5, nrow = 2,
                    labels = c("a", "b", "c", "d", "e"))
  
  return(plot)
}

suppfig16 <- asv_plotter_sunapee_preTurb_oneSite(sunapee, sun01Jul21.HC, sun01Jul21.HC, sunapee_all, 16)
ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig16.png", suppfig16,
       height = 5, width = 12)

asv_plotter_sunapee_postTurb_oneSite <- function(map_bounds, deployment, path1, lake_all, zoom){ #where deployment is the dataframe of that particular deployment date
  #and lake_all is a combined dataframe that contains all deployments for a single lake
  #and map_bounds is a set boundary of the spatial extent to be mapped
  # and zoom is either 12 for SAB or 14 for AUB, CHN
  
  a <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  aa <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  b <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  bb <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  c <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  cc <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Sp. C. (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Sp. C. (", "\u00B5", "S/cm)"))
  
  d <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  dd <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chl. ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chl. ", italic("a"), " (RFU)")))
  
  e <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ee <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("DO (mg/L)")
  
  f <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turb (NTU)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ff <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turbidity (NTU)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Turbidity (NTU)")
  
  plot <- ggarrange(aa, bb, cc, dd, ee, ff,
                    a, b, c, d, e, f,
                    ncol=6, nrow = 2,
                    labels = c("a", "b", "c", "d", "e", "f"))
  
  return(plot)
}

suppfig17 <- asv_plotter_sunapee_postTurb_oneSite(sunapee, sun22Jul21.HC, sun22Jul21.HC, sunapee_all, 16)

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig17.png", suppfig17,
       height = 5, width = 12)

asv_plotter_sunapee_postTurb <- function(map_bounds, deployment, path1, path2, lake_all, zoom){ #where deployment is the dataframe of that particular deployment date
  #and lake_all is a combined dataframe that contains all deployments for a single lake
  #and map_bounds is a set boundary of the spatial extent to be mapped
  # and zoom is either 12 for SAB or 14 for AUB, CHN
  
  a <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = pH), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  aa <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  aaa <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = pH, color = pH), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$pH, na.rm = T)-min(lake_all$pH, na.rm = T))/2)+min(lake_all$pH, na.rm = T),
                          name="pH") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, pH, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("pH")
  
  b <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = temperatureWater_degC), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  bb <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  bbb <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = temperatureWater_degC, color = temperatureWater_degC), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$temperatureWater_degC, na.rm = T)-min(lake_all$temperatureWater_degC, na.rm = T))/2)+min(lake_all$temperatureWater_degC, na.rm = T),
                          name="Temperature (ºC)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, temperatureWater_degC, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Temperature (ºC)")
  
  c <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = specificConductance_uscm), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  cc <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Specific conductance (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Sp. C. (", "\u00B5", "S/cm)"))
  
  ccc <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = specificConductance_uscm, color = specificConductance_uscm), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$specificConductance_uscm, na.rm = T)-min(lake_all$specificConductance_uscm, na.rm = T))/2)+min(lake_all$specificConductance_uscm, na.rm = T),
                          name=paste0("Sp. C. (", "\u00B5", "S/cm)")) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, specificConductance_uscm, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(paste0("Sp. C. (", "\u00B5", "S/cm)"))
  
  d <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = chlorophyll_a_RFU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chlorophyll ", italic("a"), " (RFU)"))) +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  dd <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chl. ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chl. ", italic("a"), " (RFU)")))
  
  ddd <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = chlorophyll_a_RFU, color = chlorophyll_a_RFU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$chlorophyll_a_RFU, na.rm = T)-min(lake_all$chlorophyll_a_RFU, na.rm = T))/2)+min(lake_all$chlorophyll_a_RFU, na.rm = T),
                          name=expression(paste("Chl. ", italic("a"), " (RFU)"))) +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, chlorophyll_a_RFU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab(expression(paste("Chl. ", italic("a"), " (RFU)")))
  
  e <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = oxygenDissolved_mgl), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ee <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("DO (mg/L)")
  
  eee <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = oxygenDissolved_mgl, color = oxygenDissolved_mgl), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$oxygenDissolved_mgl, na.rm = T)-min(lake_all$oxygenDissolved_mgl, na.rm = T))/2)+min(lake_all$oxygenDissolved_mgl, na.rm = T),
                          name="DO (mg/L)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, oxygenDissolved_mgl, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("DO (mg/L)")
  
  f <- get_stadiamap(map_bounds, zoom = zoom, maptype = "stamen_terrain") %>% ggmap()+
    geom_point(aes(x = longitude_gps_deg, y = latitude_gps_deg, color = turbidity_NTU), size=1, shape = 10,
               data = deployment) +
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turb (NTU)") +
    theme(legend.position = "none", axis.title=element_blank(), text = element_text(size=12, color = "black"),
          legend.key.width = unit(3, "mm"), plot.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          legend.text = element_text(angle = 45, hjust = 0, vjust = 0, size = 6), axis.text=element_text(colour="black"),
          panel.border = element_rect(color = "black", fill = NA, size = 2))
  
  ff <- ggplot() + 
    geom_point(data = path1, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turbidity (NTU)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path1, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Turbidity (NTU)")
  
  fff <- ggplot() + 
    geom_point(data = path2, aes(x = as.POSIXct(timestamp_sonde_sec, origin="1970-01-01", tz = "EST"), y = turbidity_NTU, color = turbidity_NTU), size = 3)+
    scale_color_gradient2(low = 'yellow', mid = 'purple', high = 'blue', 
                          midpoint = ((max(lake_all$turbidity_NTU, na.rm = T)-min(lake_all$turbidity_NTU, na.rm = T))/2)+min(lake_all$turbidity_NTU, na.rm = T),
                          name="Turbidity (NTU)") +
    theme_bw()+
    xlab("Timestamp") + 
    geom_line(data = asv_window(path2, turbidity_NTU, 120), aes(x = as.POSIXct(date, origin="1970-01-01", tz="EST"), y = mean), color = 'black', size=1)+
    theme(legend.position = "none", panel.background=element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(), text=element_text(size=12, color = "black"),
          axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 2))+
    scale_x_datetime(breaks = scales::date_breaks("40 mins"), date_labels = "%H:%M")+
    xlab("Time (EDT)")+
    ylab("Turbidity (NTU)")
  
  plot <- ggarrange(aa, bb, cc, dd, ee, ff,
                    aaa, bbb, ccc, ddd, eee, fff,
                    a, b, c, d, e, f,
                    ncol=6, nrow = 3,
                    labels = c("a", "b", "c", "d", "e", "f"))
  
  return(plot)
}

suppfig18 <- asv_plotter_sunapee_postTurb(sunapee, sun05Aug21, sun05Aug21.HC, sun05Aug21.NW, sunapee_all, 16)

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig18.png", suppfig18,
       height = 8, width = 14)

suppfig19 <- asv_plotter_sunapee_postTurb(sunapee, sun10Aug21, sun10Aug21.HC, sun10Aug21.NW, sunapee_all, 16)

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig19.png", suppfig19,
       height = 8, width = 14)

suppfig20 <- asv_plotter_sunapee_postTurb(sunapee, sun27Aug21, sun27Aug21.HC, sun27Aug21.NW, sunapee_all, 16)

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig20.png", suppfig20,
       height = 8, width = 14)

suppfig21 <- asv_plotter_sunapee_postTurb(sunapee, sun16Sep21, sun16Sep21.HC, sun16Sep21.NW, sunapee_all, 16)

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig21.png", suppfig21,
       height = 8, width = 14)

# Supp Fig 22

precip <- read.csv(file = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/data/precipitation/PrecipAll_LebanonMunicipalAirport.csv", header = T)

precip$DATE <- as.Date(precip$DATE, format = "%m/%d/%y")

precip_2021 <- precip %>% 
  mutate(year = year(DATE)) %>% 
  filter(year == 2021) %>% 
  filter(DATE > "2021-04-01" & DATE < "2021-10-01")

deployment_dates <- data.frame(date = as.Date(c("2021-06-11", "2021-06-15", "2021-06-21", 
                                                "2021-07-01", "2021-07-22", "2021-07-22", 
                                                "2021-08-10", "2021-08-27", "2021-09-16")))

suppfig22 <- ggplot() +
  geom_point(data = precip_2021, aes(x = DATE, y = PRCP)) +
  geom_path(data = precip_2021, aes(x = DATE, y = PRCP)) +
  scale_x_date(date_labels = "%b", date_breaks = "months") +
  geom_segment(data = deployment_dates,
               aes(x = date, xend = date, y = 0, yend = 29),
               linetype = "dashed", color = "navy") +
  annotate("text", x = deployment_dates$date, y = 32, label = deployment_dates$date, 
           angle = 90, size = 5) +
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black"),
        legend.key = element_rect(fill = NA), text=element_text(size=18, color = "black"), legend.title = element_blank())+
  xlab("") +
  ylab("Precipitation (mm)")

ggsave("~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig22.png", suppfig22, height = 8, width = 10)

# Supp Fig 23

historical_spc <- read.csv("~/Dropbox/Bates/Manuscripts/ASVLimno_MS/data/LSPALMP_1986-2020_v2021-03-29.csv", header = T)

suppfig23 <- historical_spc %>% 
  filter(station == 505 | station == 510 | station == 830 | station == 835) %>% 
  filter(parameter == "cond_uScm") %>% 
  ggplot()+
  geom_boxplot(aes(x=as.factor(station), y=value, fill=as.factor(station)))+
  scale_fill_manual(values = c("#77AADD", "#44BB99", "#EE8866", "#FFAABB"))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black"),
        legend.key = element_rect(fill = NA), text=element_text(size=18, color = "black"), legend.title = element_blank())+
  xlab("")+
  ylab("Conductivity" ~ "(\u00B5S/cm)")

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig23.png", height = 8, width = 10, suppfig23)

# Supp Fig 24

cond.sun <- read.csv("~/Dropbox/Bates/Manuscripts/ASVLimno_MS/data/conductivity/master_database_10Nov2009_v24Oct2023.csv", header = T)

na.cl.plot <- ggplot()+
  geom_point(data = cond.sun, aes(x=Na_mgl, Cl_mgl, color=cond_uScm), size=4)+
  scale_colour_gradient(low = "blue", high = "red") +
  theme(legend.position = c(0.9, 0.3), panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black"),
        legend.key = element_rect(fill = NA), text=element_text(size=16, color = "black"), legend.title = element_blank())+
  xlab("Na (mg/L)")+
  ylab("Cl (mg/L)")

cond.cl.plot <- ggplot()+
  geom_point(data = cond.sun, aes(x=cond_uScm, Cl_mgl), shape=15, size=4)+
  theme(legend.position = "bottom", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black"),
        legend.key = element_rect(fill = NA), text=element_text(size=16, color = "black"), legend.title = element_blank())+
  xlab("Conductivity" ~ "(\u00B5S/cm)")+
  ylab("Cl (mg/L)")

cond.na.plot <- ggplot()+
  geom_point(data = cond.sun, aes(x=cond_uScm, Na_mgl), shape=17, size=4)+
  theme(legend.position = "bottom", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black"),
        legend.key = element_rect(fill = NA), text=element_text(size=16, color = "black"), legend.title = element_blank())+
  xlab("Conductivity" ~ "(\u00B5S/cm)")+
  ylab("Na (mg/L)")

suppfig24 <- ggarrange(na.cl.plot, arrangeGrob(cond.cl.plot, cond.na.plot, ncol=2), ncol=1, nrow = 2)

ggsave(filename = "~/Dropbox/Bates/Manuscripts/ASVLimno_MS/documents/SuppFig24.png", height = 8, width = 10, suppfig24)

cor.test(x = cond.sun$Na_mgl, y = cond.sun$Cl_mgl)
cor.test(x = cond.sun$cond_uScm, y = cond.sun$Cl_mgl)
cor.test(x = cond.sun$cond_uScm, y = cond.sun$Na_mgl)

#######################################################################