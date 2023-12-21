#######################################################################
# ASV SPATIOTEMPORAL HETEROGENEITY -- SETUP
# Script written by Arsenault, Steele, Shingai, Cottingham et al. 2023
# Reads in data and loads libraries required for subsequent scripts
#######################################################################

##### LOAD LIBRARIES

library(tidyverse)
library(timetk)
library(lubridate)
library(zoo)
library(proj4)
library(rgdal)
library(rollply)
library(remotes)
library(deming)
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

##### LOAD DATA FILES

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

#######################################################################