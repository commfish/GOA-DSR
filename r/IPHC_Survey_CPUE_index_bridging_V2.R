################################################################################
# Aaron Lambert, ADFG
# June 23, 2026
# 
# NOTE:
# This script is copied over from the DSR Yelloweye repo and was created by Caitlin Stern
# from a Phil Joy script. I have gone through it and it is working as it should. 
# 
# Purpose:
# This code is used to process and organize the IPHC survey data for use
# as a secondary index of yelloweye rockfish abundance in the SEO. This specific 
# version is to create IPHC indices of specific flavors for the bridging analysis:
# 
# # The seven bridging models for IPHC CPUE
#
#  mod.0: "22.2"
#         nominal CPUE in numbers per hook
#
#  mod.1: "22.2 + hook adj."
#         nominal CPUE in numbers per hook, 
#         adjusted for hook saturation and competition
#         
#  mod.2: "22.2 + lon"
#         nominal CPUE in numbers per hook, 
#         adjusted for hook saturation and competition, 
#         with an error in the definition of the Eastern Yakutat (EYKT)
#          management district corrected (140 degrees longitude is the correct boundary, 
#          not 137 degrees as was used in the 2022 assessment)
#
#  mod.3: "22.2 + kg/hook"
#         nominal CPUE in kg per hook, 
#         adjusted for hook saturation and competition, 
#         with the correct EYKT boundary 
#  
#  mod.4  "22.2+ correct stations"
#         nominal CPUE in kg per hook, 
#         adjusted for hook saturation and competition, 
#         with the correct EYKT boundary &
#         the correct lat for EYKT and omitting SEI stations
#         
#  mod.5: "22.2 + std. CPUE"
#         CPUE in kg per hook standardized using a GAM model with the
#         Tweedie distribution in order to accommodate zero inflation, 
#         adjusted for hook saturation and competition, 
#         with the correct EYKT boundary
#   
#  mod.6: "22.2 + std  CPUE + correct stations"
#         CPUE in kg per hook standardized using a GAM model with the
#         Tweedie distribution in order to accommodate zero inflation, 
#         adjusted for hook saturation and competition, 
#         with the correct EYKT boundary & 
#         the correct lat for EYKT and omitting SEI stations.
#
# Notes from the past drivers:
# Methods include the original bootstrap methods used in the 2022 assessment
# and the new Tweedie estimator developed after the 2023 CIE review
# Phil Joy
# Oct. 2023


# Updated 9/24/24 by RKE and fixed some errors and issues 
# you will need to fix file location because I did not push all of my files to github
# and i set my folder structure up differently due to confusion with current organization 

################################################################################
# Set the year of the assessment
YEAR <- 2026

library(dplyr)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(readr)
library(vroom)
library(scales)
library(ggpubr)
library(mgcv)
library(MuMIn)
library(GGally)
library(mgcViz)
library(wesanderson)
library(RColorBrewer)
library(tweedie)
library(here)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# Deprecated function
# source(file = here("r_helper/Port_bio_function.R"))

#*******************************************************************************
# Import data ------------------------------------------------------------------
# RKE truncated this code and cleaned up so you only need to import the full data files from IPHC once and
# avoid join issues 

# Halibut catch by station. Needed to join with Yelloweye catch to get all the 
# 0 catch stations. Additionally, this has potential covariates for the GAM model
# such as pH, temp, etc...
BUT2C3A <- read_csv(here("data/raw/Set and Pacific halibut data 4June2026.csv")) %>% 
  dplyr::rename_all(~make.names(.))

# YE catch 
YE2C3A <- read_csv(here("data/raw/Non-Pacific halibut data 4June2026.csv")) %>% 
  dplyr::rename_all(~make.names(.))


# RKE corrected EYKT bounds to 137-140 Longitude as that is the management area
# Phil had 137-139 so it was exluding 1/3 of the area. Not sure if that was intentional or an error 
# or if this change impacts the overall results at all 
Survey_prep <- BUT2C3A %>% 
  full_join(YE2C3A, by=c("Stlkey","Station")) %>% 
  filter(IPHC.Stat.Area <= 200, 
         -MidLon.fished <= 140, 
         Eff == "Y") %>%  
  mutate(SEdist = ifelse(IPHC.Stat.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                         ifelse(IPHC.Stat.Area %in% c(142:144,152,153),"SSEI",
                                ifelse(-MidLon.fished >= 137 & -MidLon.fished <= 140, "EYKT", 
                                       ifelse(-MidLon.fished <= 137 & MidLat.fished >= 57.5,"NSEO",
                                              ifelse(-MidLon.fished <= 137 & MidLat.fished <57.5 & MidLat.fished >= 56,"CSEO",
                                                     ifelse(-MidLon.fished <= 137 & MidLat.fished < 56,"SSEO",NA)))))), 
         In.Out = ifelse(SEdist %in% c("NSEO","CSEO","SSEO","EYKT"),"SEO",
                         ifelse(SEdist %in% c("NSEI","SSEI"),"SEI",NA))) %>% 
  rename(YE.obs = Number.Observed,
         AvgDepth.fm = AvgDepth..fm.,
         Year = Year.x) %>%
  mutate(depth_bin = cut(AvgDepth.fm, breaks = seq(0,400,50),
                         labels = paste (seq(50,400,50))),
         depth_m = AvgDepth.fm * 1.8288,
         depth_bin_m = cut(depth_m, breaks = seq(0,800,50),
                           labels = paste (seq(50,800,50)))) %>% 
  select(-c(Row.number.x,Row.number.y))

# Using Phils old Longitude for the 2024 (presented in 2026) bridging analysis (mod.2.data)
Survey_prep_incorrectLON <- BUT2C3A %>% 
  full_join(YE2C3A, by=c("Stlkey","Station")) %>% 
  filter(IPHC.Stat.Area <= 200, 
         -MidLon.fished <= 140, 
         Eff == "Y") %>%  
  mutate(SEdist = ifelse(IPHC.Stat.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                         ifelse(IPHC.Stat.Area %in% c(142:144,152,153),"SSEI",
                                ifelse(-MidLon.fished >= 137 & -MidLon.fished <= 139, "EYKT", 
                                       ifelse(-MidLon.fished <= 137 & MidLat.fished >= 57.5,"NSEO",
                                              ifelse(-MidLon.fished <= 137 & MidLat.fished <57.5 & MidLat.fished >= 56,"CSEO",
                                                     ifelse(-MidLon.fished <= 137 & MidLat.fished < 56,"SSEO",NA)))))), 
         In.Out = ifelse(SEdist %in% c("NSEO","CSEO","SSEO","EYKT"),"SEO",
                         ifelse(SEdist %in% c("NSEI","SSEI"),"SEI",NA))) %>% 
  rename(YE.obs = Number.Observed,
         AvgDepth.fm = AvgDepth..fm.,
         Year = Year.x) %>%
  mutate(depth_bin = cut(AvgDepth.fm, breaks = seq(0,400,50),
                         labels = paste (seq(50,400,50))),
         depth_m = AvgDepth.fm * 1.8288,
         depth_bin_m = cut(depth_m, breaks = seq(0,800,50),
                           labels = paste (seq(50,800,50)))) %>% 
  select(-c(Row.number.x,Row.number.y))

# This is the current corrected (2026) calculations that fix the EYKT lat and omit SEI stations that were in the 2024 and earlier assessments
Survey_prep_2026correct <- BUT2C3A %>% 
  full_join(YE2C3A, by=c("Stlkey","Station")) %>% 
  filter(IPHC.Stat.Area <= 200, 
         -MidLon.fished <= 140, 
         !Station %in% c(3120, 3072, 3067,  3071, 3061, 3060, 3059, 3065, 3216, 3213, 3212), # AWL and RKE added this, but if this is ever used for SEI, then it will need to be reassessed...
         Eff == "Y") %>%  
  mutate(SEdist = ifelse(IPHC.Stat.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                         ifelse(IPHC.Stat.Area %in% c(143:144,152,153),"SSEI", # AWL and RKE changed to 143 from 142
                                ifelse(-MidLon.fished >= 137 & -MidLon.fished <= 140 & MidLat.fished <= 59.5, "EYKT", 
                                       ifelse(-MidLon.fished <= 137 & MidLat.fished >= 57.5,"NSEO",
                                              ifelse(-MidLon.fished <= 137 & -MidLon.fished >=134.6 & MidLat.fished <57.5 & MidLat.fished >= 56,"CSEO",
                                                     ifelse(-MidLon.fished <= 137 & -MidLon.fished>= 132.69 & MidLat.fished < 56,"SSEO",NA)))))), 
         In.Out = ifelse(SEdist %in% c("NSEO","CSEO","SSEO","EYKT"),"SEO",
                         ifelse(SEdist %in% c("NSEI","SSEI"),"SEI",NA))) %>% 
  rename(YE.obs = Number.Observed,
         AvgDepth.fm = AvgDepth..fm.,
         Year = Year.x) %>%
  mutate(depth_bin = cut(AvgDepth.fm, breaks = seq(0,400,50),
                         labels = paste (seq(50,400,50))),
         depth_m = AvgDepth.fm * 1.8288,
         depth_bin_m = cut(depth_m, breaks = seq(0,800,50),
                           labels = paste (seq(50,800,50)))) %>% 
  select(-c(Row.number.x,Row.number.y))

# Bring in and join hook adjustment factor to control for hook saturation when calculating CPUE
# Note from PJ: The data was already posted: https://www.iphc.int/data/fiss-survey-raw-survey-data/
hadj <- read_csv(here("data/raw/iphc-2025-fiss-hadj.csv")) %>%
  # hadj <- read_csv(here("data/raw/iphc-2024-fiss-hadj.csv")) %>% 
  dplyr::rename_all(~(make.names(.))) %>% 
  filter(IPHC.Reg.Area %in% c("2C", "3A"), 
         IPHC.Stat.Area <= 200)  %>%
  select(Stlkey = stlkey, Year, Station, AdjHooks.Observed = Hooks.Observed, h.adj)

str(hadj)

# Join the survey prep and hadj factor
Survey <- left_join(Survey_prep,hadj,by=c("Year","Stlkey","Station")) %>% 
  mutate(YE.obs = ifelse(is.na(YE.obs),0,YE.obs),
         # RKE used AdjHooks.Observed from hadj file to fill in missing 
         # hook observed fields instead of static "140" from Phil's code which is unclear 
         # where that came from given a different number of skates are sampled 
         # he said 20 hooks observed x 7 skates sampled but some sets are a total of 6 skates set/hauled
         # so those would be 120 hooks instead of 140. Regardless, if you use the Hooks Observed
         # from the hadj file, it's IPHC data for that Stlkey/Station! So I assume it's correct. 
         HooksObserved = ifelse(is.na(HooksObserved), AdjHooks.Observed, HooksObserved), 
         HooksRetrieved = as.numeric(HooksRetrieved), 
         HooksObserved = as.numeric(HooksObserved), 
         YE.exp = HooksRetrieved/HooksObserved)

# Join the survey prep and hadj factor
Survey.wronglon <- left_join(Survey_prep_incorrectLON,hadj,by=c("Year","Stlkey","Station")) %>% 
  mutate(YE.obs = ifelse(is.na(YE.obs),0,YE.obs),
         # RKE used AdjHooks.Observed from hadj file to fill in missing 
         # hook observed fields instead of static "140" from Phil's code which is unclear 
         # where that came from given a different number of skates are sampled 
         # he said 20 hooks observed x 7 skates sampled but some sets are a total of 6 skates set/hauled
         # so those would be 120 hooks instead of 140. Regardless, if you use the Hooks Observed
         # from the hadj file, it's IPHC data for that Stlkey/Station! So I assume it's correct. 
         HooksObserved = ifelse(is.na(HooksObserved), AdjHooks.Observed, HooksObserved), 
         HooksRetrieved = as.numeric(HooksRetrieved), 
         HooksObserved = as.numeric(HooksObserved), 
         YE.exp = HooksRetrieved/HooksObserved)

Survey_2026_correct <- left_join(Survey_prep_2026correct ,hadj,by=c("Year","Stlkey","Station")) %>% 
  mutate(YE.obs = ifelse(is.na(YE.obs),0,YE.obs),
         # RKE used AdjHooks.Observed from hadj file to fill in missing 
         # hook observed fields instead of static "140" from Phil's code which is unclear 
         # where that came from given a different number of skates are sampled 
         # he said 20 hooks observed x 7 skates sampled but some sets are a total of 6 skates set/hauled
         # so those would be 120 hooks instead of 140. Regardless, if you use the Hooks Observed
         # from the hadj file, it's IPHC data for that Stlkey/Station! So I assume it's correct. 
         HooksObserved = ifelse(is.na(HooksObserved), AdjHooks.Observed, HooksObserved), 
         HooksRetrieved = as.numeric(HooksRetrieved), 
         HooksObserved = as.numeric(HooksObserved), 
         YE.exp = HooksRetrieved/HooksObserved)

# # Sanity check. Should return 1 (indicating that each obs has 1 row)
# {y<-sample(unique(Survey$Year),1)
#   d<-sample(unique(Survey$depth_bin[Survey$Year==y]),1)
#   s<-sample(unique(Survey$Station[Survey$Year==y & Survey$depth_bin==d]),1)
#   Dat<-Survey[Survey$Year == y & Survey$depth_bin == d & Survey$Station == s,]
#   unique(Dat); nrow(Dat)}
# 
# str(Survey); head(Survey,10)
# 
# # str(HA.Harv)
# mean(Survey$AvgDepth.fm)
# 
# # write.csv(Survey,"./Data_processing/Data/IPHC_survey_1998-2024_cas.csv")
# 
# # Exploratory data analysis ----------------------------------------------------
# # Look at some depth stuff and identify stations that have seen yelloweye...
# hist(Survey$BeginDepth..fm.)
# hist(Survey$EndDepth..fm.)
# hist(Survey$AvgDepth.fm)
# 
# plot(Survey$YE.obs ~ Survey$AvgDepth.fm)
# plot(Survey$O32.Pacific.halibut.count ~ Survey$AvgDepth.fm)
# plot(Survey$U32.Pacific.halibut.count ~ Survey$AvgDepth.fm)
# 
# # This is a table of ye catch by depth strata
# # Some data for PWS YE assessment: 
# with(Survey, table(YE.obs,depth_bin_m))
# ye_by_depth_m <- Survey %>% 
#   group_by(depth_bin_m) %>% 
#   summarize(ye_obs = sum(YE.obs), n_stations = n_distinct(Station)) %>% 
#   mutate(YE_per_station = ye_obs / n_stations)
# 
# #write.csv(ye_by_depth_m,"Data_processing/Data/SE_YE_FISS_depth_forPWS.csv")
# 
# Summary of obs YE data by station
Station.sum <- Survey %>% group_by(Station) %>%
  dplyr::summarise(Lat = mean(MidLat.fished, na.rm=T),
                   Lon = mean(MidLon.fished, na.rm=T),
                   Depth = mean(AvgDepth.fm, na.rm=T),
                   obs = n(),
                   # mean.hal.O32.count = mean(O32.Pacific.halibut.count, na.rm=T),
                   # mean.hal.O32.wt = mean(O32.Pacific.halibut.weight, na.rm=T),
                   # mean.hal.U32.count = mean(U32.Pacific.halibut.count, na.rm=T),
                   # mean.hal.U32.wt = mean(U32.Pacific.halibut.weight, na.rm=T),
                   mean.skates = mean(No..skates.set, na.rm=T),
                   mean.hooks.ob = mean(HooksObserved, na.rm=T),
                   mean.hooks.ret = mean(HooksRetrieved, na.rm=T),
                   mean.YE = mean(YE.obs, na.rm=T),
                   var.YE = var(YE.obs, na.rm=T),
                   mean.YEcpue = mean(YE.obs/HooksObserved, na.rm=T),
                   var.YEcpue = var(YE.obs/HooksObserved, na.rm=T),
                   noYE.count = sum(YE.obs == 0),
                   prop.0 = noYE.count/obs
  )

# After building Survey_2026_correct but BEFORE filtering to Survey_non0_2026,
# rebuild Station.sum from the corrected survey
Station.sum_2026 <- Survey_2026_correct %>%
  group_by(Station) %>%
  dplyr::summarise(
    Lat = mean(MidLat.fished, na.rm = TRUE),
    Lon = mean(MidLon.fished, na.rm = TRUE),
    Depth = mean(AvgDepth.fm, na.rm = TRUE),
    obs = n(),
    mean.YE = mean(YE.obs, na.rm = TRUE),
    mean.YEcpue = mean(YE.obs / HooksObserved, na.rm = TRUE),
    noYE.count  = sum(YE.obs == 0),
    prop.0 = noYE.count / obs
  )

YE.stations_2026 <- Station.sum_2026$Station[Station.sum_2026$mean.YEcpue != 0]
# 
# nrow(Station.sum)
# 
# # Histogram of mean YE catch at stations
# hist(Station.sum$mean.YE, breaks=50) 
# 
# # Proportion of stations that encountered YE at least once
# 1-(nrow(Station.sum[Station.sum$mean.YE == 0,])/nrow(Station.sum))
# 
# # Mean station YE catch by depth with 0 catch stations
# ggplot(Station.sum,aes(Depth, mean.YE)) +
#   geom_point() + geom_smooth(linewidth = 1.5, se = FALSE, linetype=1, alpha=0.75) +
#   ylab("Mean number of yelloweye encountered at FISS stations") +
#   xlab("Depth (fathoms)") + theme_bw() +
#   scale_x_continuous(breaks=seq(0,350,25)) 
# 
# # ggsave(paste0("output/FISS_YE_x_depth_cas_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in") 
# 
# # Mean station YE catch by depth, excluding 0 catch stations
# ggplot(Station.sum %>% filter(mean.YE > 0),
#        aes(Depth, mean.YE)) +
#   geom_point() + geom_smooth(size = 1.5, se = FALSE, linetype=1, alpha=0.75) +
#   ylab("Mean number of yelloweye encountered at FISS stations") +
#   xlab("Depth (fathoms)") + theme_bw() +
#   scale_x_continuous(breaks=seq(0,350,25))
# 
# # ggsave(paste0("output/FISS_YE_x_depth_non0_cas_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in") 
# 
# ye_caught_prop<-percent(1-(nrow(Station.sum[Station.sum$mean.YE == 0,])/nrow(Station.sum)))
# 
# # Plot showing dist of mean YE catch at stations with catch prop printed
# all_stations_hist <- ggplot(Station.sum) +
#   geom_histogram(aes(mean.YE), binwidth=0.25, color="grey", fill="darkgrey") +
#   annotate("text", y=65, # 
#            x=12.5, color="darkorange",
#            label=c(paste0("Yelloweye caught at least once at ",ye_caught_prop," of stations"))) +
#   annotate("text", y=75, 
#            x=12.5, color="black",fontface =2, 
#            label=c(paste0("All SEO FISS stations"))) +
#   ylab("Number of FISS stations") +
#   xlab("") + 
#   theme_bw()
# 
# all_stations_hist
# 
# # Same plot as above with station means = 0 removed
# ggplot(Station.sum %>% filter(mean.YE > 0)) +
#   geom_histogram(aes(mean.YE), binwidth=0.25, color="grey", fill="darkgrey") +
#   annotate("text", y=15, 
#            x=12.5, color="black",fontface =2,
#            label=c(paste0("SEO FISS stations that encountered yelloweye at least once"))) +
#   ylab("Number of FISS stations") +
#   xlab("mean number of yelloweye encountered per year") + theme_bw() -> non0_stations_hist
# 
# ggarrange(all_stations_hist, non0_stations_hist,
#           nrow=2)
# 
# ggsave(paste0("outputs/FISS_YE_histograms_cas_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")
# 
# hist(Station.sum$mean.YEcpue, breaks=25); nrow(Station.sum[Station.sum$mean.YEcpue == 0,])/nrow(Station.sum)
# hist(Station.sum$noYE.count, breaks=25)
# hist(Station.sum$prop.0, breaks=25)
# hist(Station.sum$prop.0[Station.sum$prop.0<1], breaks=25)
# max(Station.sum$prop.0[Station.sum$prop.0<1])
# 
# 
# blanks.stations<-Station.sum$Station[Station.sum$mean.YEcpue == 0]
# 
# # This is used below to subset stations. All stations is used from 2024 assessment
# # and onward, as a CIE review recommended including 0 catch observations when standarziing CPUE
# All.stations <- Station.sum$Station
# 
# # This is used to subset stations with YE catch. This was how the 
# # CPUE was calculated pre 2024. 
# YE.stations<-Station.sum$Station[Station.sum$mean.YEcpue != 0]
# 
# # Look at proportions of stations with low or not catch over years
# YE.stations_10p<-Station.sum$Station[Station.sum$prop.0 < 0.9]
# YE.stations_20p<-Station.sum$Station[Station.sum$prop.0 < 0.8]
# YE.stations_25p<-Station.sum$Station[Station.sum$prop.0 < 0.75]
# YE.stations_40p<-Station.sum$Station[Station.sum$prop.0 < 0.6]
# 
# length(All.stations);length(YE.stations);length(YE.stations_10p);length(YE.stations_20p);length(YE.stations_25p);length(YE.stations_40p)
# 
# plot_stations <- st_as_sf(Station.sum,coords = c("Lon","Lat"))
# 
# plot_stations <- plot_stations %>% st_set_crs(4326)
# 
# # Get land polygons
# land <- ne_countries(scale = "medium", returnclass = "sf")
# states <- ne_states(country = c("united states of america", "canada"), returnclass = "sf")
# 
# # Define bounding box for SE Alaska
# se_ak_bbox <- c(xmin = -140, xmax = -129, ymin = 54, ymax = 60)
# 
# # Function to add land to station plots
# plot_with_land <- function(data, color_var, size_var, palette, plot_title) {
#   ggplot() +
#     geom_sf(data = land, fill = "grey80", color = "grey60", linewidth = 0.3) +
#     geom_sf(data = states, fill = "grey80", color = "grey60", linewidth = 0.2) +
#     geom_sf(data = data %>% filter(Station %in% use_stations),
#             aes(color = .data[[color_var]], size = .data[[size_var]]),
#             alpha = 0.7) +
#     scale_color_viridis_c(option = palette, begin = 0, name = color_var) +
#     coord_sf(xlim = c(se_ak_bbox["xmin"], se_ak_bbox["xmax"]),
#              ylim = c(se_ak_bbox["ymin"], se_ak_bbox["ymax"])) +
#     annotation_scale(location = "bl") +
#     labs(title = plot_title) +
#     theme_bw() +
#     theme(panel.grid.major = element_line(color = "grey90"),
#           legend.position = "right")
# }
# 
# use_stations<-All.stations
# 
# # Mean YE catch
# plot_with_land(plot_stations, "mean.YE", "mean.YE", "magma", "Mean yelloweye catch per station")
# 
# # Mean YE CPUE
# plot_with_land(plot_stations, "mean.YEcpue", "mean.YEcpue", "plasma", "Mean yelloweye CPUE per station")
# 
# # Variance in YE CPUE
# plot_with_land(plot_stations, "var.YEcpue", "var.YEcpue", "viridis", "Variance in yelloweye CPUE per station")
# 
# # Proportion zero catches
# plot_with_land(plot_stations, "prop.0", "prop.0", "plasma", "Proportion of zero catches per station")
# 
# # # use all stations
# # #use_stations<-YE.stations
# # use_stations<-All.stations
# # 
# # ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
# #   geom_sf(aes(color = mean.YE, size = mean.YE)) +
# #   scale_color_viridis_c(option = "magma",begin = 0)
# # 
# # ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
# #   geom_sf(aes(color = mean.YEcpue, size = mean.YEcpue)) +
# #   scale_color_viridis_c(option = "plasma",begin = 0)
# # 
# # ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
# #   geom_sf(aes(color = var.YEcpue, size = var.YEcpue)) +
# #   scale_color_viridis_c(option = "viridis",begin = 0)
# # 
# # ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
# #   geom_sf(aes(color = prop.0, size = prop.0)) +
# #   scale_color_viridis_c(option = "plasma",begin = 0)

#*********************************************************************************
# Import port sample weight to get WCPUE ---------------------------------------
# Need Yelloweye weights to get wcpue estimates
#********************************************************************************* 

# Load the port samples that were downloaded from oceansAK; 
# Port sampling biological data ------------------------------------------------
Port1 <- read.csv("data/port sampling data/SEO_YE_port_sampling_bio_data_1980-1989.csv")
Port2 <- read.csv("data/port sampling data/SEO_YE_port_sampling_bio_data_1990-1999.csv")
Port3 <- read.csv("data/port sampling data/SEO_YE_port_sampling_bio_data_2000-2009.csv")
Port4 <- read.csv("data/port sampling data/SEO_YE_port_sampling_bio_data_2010-2019.csv")
Port5 <- read.csv(paste0("data/port sampling data/SEO_YE_port_sampling_bio_data_2020-", 2025, ".csv"))

Port <- rbind(Port1, Port2, Port3, Port4, Port5)
Port$Year <- as.integer(Port[,1])

# Recode EYAK to EYKT for consistency
Port$Groundfish.Management.Area.Code[Port$Groundfish.Management.Area.Code == "EYAK"] <- "EYKT"

# Assign groundfish management units using stat areas where GFMA code is missing
statareas <- read.csv("data/port sampling data/g_stat_area.csv")

SSc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEO"])
CSc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "CSEO"])
NSc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEO"])
EYc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "EYKT"])

Port <- Port %>% 
  mutate(GFMU = ifelse(Groundfish.Management.Area.Code == "",
                       ifelse(Groundfish.Stat.Area %in% SSc |
                                substr(Groundfish.Stat.Area.Group, 1, 4) %in% "SSEO", "SSEO",
                              ifelse(Groundfish.Stat.Area %in% CSc |
                                       substr(Groundfish.Stat.Area.Group, 1, 4) %in% "CSEO", "CSEO",
                                     ifelse(Groundfish.Stat.Area %in% NSc |
                                              substr(Groundfish.Stat.Area.Group, 1, 4) %in% "NSEO", "NSEO",
                                            ifelse(Groundfish.Stat.Area %in% EYc |
                                                     substr(Groundfish.Stat.Area.Group, 1, 4) %in% "EYKT", "EYKT", NA)))),
                       Groundfish.Management.Area.Code),
         Sex = case_when(Sex.Code == 1 ~ "Male",
                         Sex.Code == 2 ~ "Female"),
         Sex = as.factor(Sex))

# Filter to random samples only and remove missing weights
# Rhea edited so that samples are Random. Want to avoid the
#  select samples for maturity/vonB that were collected 
Port <- Port %>% 
  filter(Sample.Type == "Random") %>%
  filter(!is.na(Weight.Kilograms)) %>% 
  filter(Groundfish.Management.Area.Code != "")

str(Port)
unique(Port$Groundfish.Management.Area.Code)  #EYAK = EYKT

# Port<-subset(Port, !is.na(Weight.Kilograms))

# Get mean weights by year and management area
uYEkg <- Port %>% 
  group_by(Year,Groundfish.Management.Area.Code) %>%
  summarize(uYEkg = mean(Weight.Kilograms, na.rm=TRUE),
            vYEkg = var(Weight.Kilograms, na.rm=TRUE),
            NyeKG = length(Weight.Kilograms))


uYEkg<-as.data.frame(uYEkg); uYEkg

# Initialize columns for loop below
Survey$mean.YE.kg <- NA
Survey$var.YE.kg <- NA
Survey$N.YE.kg <-NA

Survey.wronglon$mean.YE.kg <- NA
Survey.wronglon$var.YE.kg <- NA
Survey.wronglon$N.YE.kg <-NA

Survey_2026_correct$mean.YE.kg <- NA
Survey_2026_correct$var.YE.kg <- NA
Survey_2026_correct$N.YE.kg <-NA

# Get best data for average weights and add columns to survey data
# Rhea updated this so that we get 200 samples per year per management area 
# looking at recent years (2022-2024) we have very few samples. For example,
# from EYKT there are 3 boats sampled in 2022, 3 boats in 2023, and 2 boats in 2024. I don't think 
# these are representative. For our port sampling goals, we aim for 550 per year per management area, 
# with 50-75 samples per boat to spread distribution of samples across time and area. 
Years <- unique(Survey$Year)
GFMA <- unique(Survey$SEdist)

# Loop for the correct lon
for (i in Years){    #i<-Years[1]
  for (j in GFMA){   #j<-GFMA[1]
    
    P <- Port[Port$Groundfish.Management.Area.Code == j,]
    
    #reach back from year i and see if you can get 150 weights...rhea bumped this to 200
    Nw <- nrow(P[P$Year == i,])
    
    k <- 0
    while (Nw < 200 & i-k >= min(Port$Year) ){  #go back as far as needed to get 200 weights... 
      k <- k+1
      Nw <- nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
    }
    m <- 0
    while (Nw < 200 & i+m <= max(Port$Year)) {    #go forward if you failed to find enough weights... 
      m <- m+1                                                        
      Nw <- nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
    }
    #get average weights of YE...
    Sample <- P$Weight.Kilograms[P$Year >= i-k & P$Year <= i+m]
    
    #PHIL'S CODE FAILS HERE DUE TO WEIRD SUBSETTING CODE - FIXED IN FOLLOWING 3 LINES
    # Survey[,"mean.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
    # Survey[,"var.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-var(Sample)
    # Survey[,"N.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-length(Sample)
    Survey$mean.YE.kg[Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
    Survey$var.YE.kg[Survey$Year == i & Survey$SEdist == j]<-var(Sample)
    Survey$N.YE.kg[Survey$Year == i & Survey$SEdist == j]<-length(Sample)
    
    Survey.wronglon$mean.YE.kg[Survey.wronglon$Year == i & Survey.wronglon$SEdist == j]<-mean(Sample)
    Survey.wronglon$var.YE.kg[Survey.wronglon$Year == i & Survey.wronglon$SEdist == j]<-var(Sample)
    Survey.wronglon$N.YE.kg[Survey.wronglon$Year == i & Survey.wronglon$SEdist == j]<-length(Sample)
    
    Survey_2026_correct$mean.YE.kg[Survey_2026_correct$Year == i & 
                                     Survey_2026_correct$SEdist == j] <- mean(Sample)
    Survey_2026_correct$var.YE.kg[Survey_2026_correct$Year == i & 
                                    Survey_2026_correct$SEdist == j] <- var(Sample)
    Survey_2026_correct$N.YE.kg[Survey_2026_correct$Year == i & 
                                  Survey_2026_correct$SEdist == j] <- length(Sample)
  }
}


Survey$mean.YE.kg
Survey.wronglon$mean.YE.kg

Survey %>% group_by(SEdist,Year) %>%
  dplyr::summarise(mean_wt = mean(mean.YE.kg)) %>%
  ggplot() + 
  geom_point(aes(Year,mean_wt,col=SEdist)) +
  geom_line(aes(Year,mean_wt,col=SEdist)) 


#*********************************************************************************
# Boot strap CPUE estimates
#*********************************************************************************
# FUNCTION for generating cpue estimates from IPHC surveys for different management areas

YEHA.fxn<-function(Survey=Survey, Area="SEdist",Deep=250, Shallow=0,  nboot=1000){
  col.ref<-which(colnames(Survey)==Area)
  
  IPHC.cpue<-data.frame()
  Subs<-unique(Survey$SEdist) # FIXED FROM PHIL'S CODE - DIDN'T SPECIFY COLUMN
  Years<-unique(Survey$Year)
  
  j<-1
  for (y in Years) {  #y<-Years[26]
    for (s in Subs){  #s<-Subs[3]
      Dat<-Survey[Survey$Year == y & 
                    Survey[,col.ref] == s &
                    Survey$AvgDepth.fm > Shallow & 
                    Survey$AvgDepth.fm < Deep &
                    Survey$Eff == "Y",]
      YE_w <-unique(Dat$mean.YE.kg)
      if (nrow(Dat)>0){
        Stations<-unique(Dat$Station)
        CPUEi<-vector()
        WCPUEi<-vector()
        i<-1
        for (st in Stations){    #st<-Stations[1]   length(Stations)  st<-1
          Stat.Dat<-Dat[Dat$Station == st,] #; Stat.Dat
          #debug
          #if (nrow(Stat.Dat) > 1){aaa} else {}
          CPUE<-mean(Stat.Dat$YE.obs/Stat.Dat$HooksObserved)
          #hook adjustment factor
          CPUE<-CPUE*unique(Stat.Dat$h.adj)
          WCPUE<-CPUE*YE_w
          
          if (CPUE == 0){
            C<-0
          } else {
            C<-CPUE*Stat.Dat$HooksRetrieved
          }
          
          CPUEi[i]<-CPUE
          WCPUEi[i]<-WCPUE
          i<-i+1
        }
        
        Out<-data.frame()
        #CPUE.out<-data.frame()
        for (ii in 1:nboot){ #nboot<-1000  ii<-1
          Resamp3<-sample(CPUEi,length(CPUEi),replace=T)
          Out[ii,"CPUE"]<-mean(Resamp3, na.rm=T)
          Out[ii,"CPUE.var"]<-var(Resamp3, na.rm=T)
          Resamp4<-sample(WCPUEi,length(WCPUEi),replace=T)
          Out[ii,"WCPUE"]<-mean(Resamp4, na.rm=T)
          Out[ii,"WCPUE.var"]<-var(Resamp4, na.rm=T)
        }
        
        #hist(Out$KgHa, breaks = 25)
        IPHC.cpue[j,"Year"]<-y
        IPHC.cpue[j,"mngmt.divisions"]<-Area
        IPHC.cpue[j,"mngmt.area"]<-s
        IPHC.cpue[j,"deep.bound"]<-Deep
        IPHC.cpue[j,"shallow.bound"]<-Shallow
        IPHC.cpue[j,"no.stations"]<-length(Stations)
        
        IPHC.cpue[j,"CPUE.mean"]<-mean(CPUEi)
        IPHC.cpue[j,"CPUE.bootmean"]<-unname(quantile(Out$CPUE,c(0.5)))  #mean WCPUE from Tribuzio - AWL: actually the median...
        IPHC.cpue[j,"CPUE.lo95ci"]<-unname(quantile(Out$CPUE,c(0.025)))
        IPHC.cpue[j,"CPUE.hi95ci"]<-unname(quantile(Out$CPUE,c(0.975)))
        IPHC.cpue[j,"CPUE.var"]<-var(Out$CPUE)
        IPHC.cpue[j,"CPUE.cv"]<-sd(Out$CPUE)/mean(Out$CPUE)
        #    IPHC.wcpue[j,"WCPUE32.cv"]<-sqrt(var(WCPUEi.32))/mean(WCPUEi.32)
        IPHC.cpue[j,"WCPUE.mean"]<-mean(WCPUEi)
        IPHC.cpue[j,"WCPUE.bootmean"]<-unname(quantile(Out$WCPUE,c(0.5)))  #mean WCPUE from Tribuzio - AWL: actually the median...
        IPHC.cpue[j,"WCPUE.lo95ci"]<-unname(quantile(Out$WCPUE,c(0.025)))
        IPHC.cpue[j,"WCPUE.hi95ci"]<-unname(quantile(Out$WCPUE,c(0.975)))
        IPHC.cpue[j,"WCPUE.var"]<-var(Out$WCPUE)
        IPHC.cpue[j,"WCPUE.cv"]<-sd(Out$WCPUE)/mean(Out$WCPUE)
        
        j<-j+1
      } else {}
    }
  }
  
  
  return(IPHC.cpue)
}

# Hook adjustment factor removed
YEHA.fxn.noHA<-function(Survey=Survey, Area="SEdist",Deep=250, Shallow=0,  nboot=1000){
  col.ref<-which(colnames(Survey)==Area)
  
  IPHC.cpue<-data.frame()
  Subs<-unique(Survey$SEdist) # FIXED FROM PHIL'S CODE - DIDN'T SPECIFY COLUMN
  Years<-unique(Survey$Year)
  
  j<-1
  for (y in Years) {  #y<-Years[26]
    for (s in Subs){  #s<-Subs[3]
      Dat<-Survey[Survey$Year == y & 
                    Survey[,col.ref] == s &
                    Survey$AvgDepth.fm > Shallow & 
                    Survey$AvgDepth.fm < Deep &
                    Survey$Eff == "Y",]
      YE_w <-unique(Dat$mean.YE.kg)
      if (nrow(Dat)>0){
        Stations<-unique(Dat$Station)
        CPUEi<-vector()
        WCPUEi<-vector()
        i<-1
        for (st in Stations){    #st<-Stations[1]   length(Stations)  st<-1
          Stat.Dat<-Dat[Dat$Station == st,] #; Stat.Dat
          #debug
          #if (nrow(Stat.Dat) > 1){aaa} else {}
          CPUE<-mean(Stat.Dat$YE.obs/Stat.Dat$HooksObserved)
          
          WCPUE<-CPUE*YE_w
          
          if (CPUE == 0){
            C<-0
          } else {
            C<-CPUE*Stat.Dat$HooksRetrieved
          }
          
          CPUEi[i]<-CPUE
          WCPUEi[i]<-WCPUE
          i<-i+1
        }
        
        Out<-data.frame()
        #CPUE.out<-data.frame()
        for (ii in 1:nboot){ #nboot<-1000  ii<-1
          Resamp3<-sample(CPUEi,length(CPUEi),replace=T)
          Out[ii,"CPUE"]<-mean(Resamp3, na.rm=T)
          Out[ii,"CPUE.var"]<-var(Resamp3, na.rm=T)
          Resamp4<-sample(WCPUEi,length(WCPUEi),replace=T)
          Out[ii,"WCPUE"]<-mean(Resamp4, na.rm=T)
          Out[ii,"WCPUE.var"]<-var(Resamp4, na.rm=T)
        }
        
        #hist(Out$KgHa, breaks = 25)
        IPHC.cpue[j,"Year"]<-y
        IPHC.cpue[j,"mngmt.divisions"]<-Area
        IPHC.cpue[j,"mngmt.area"]<-s
        IPHC.cpue[j,"deep.bound"]<-Deep
        IPHC.cpue[j,"shallow.bound"]<-Shallow
        IPHC.cpue[j,"no.stations"]<-length(Stations)
        
        IPHC.cpue[j,"CPUE.mean"]<-mean(CPUEi)
        IPHC.cpue[j,"CPUE.bootmean"]<-unname(quantile(Out$CPUE,c(0.5)))  #mean WCPUE from Tribuzio - AWL: actually the median...
        IPHC.cpue[j,"CPUE.lo95ci"]<-unname(quantile(Out$CPUE,c(0.025)))
        IPHC.cpue[j,"CPUE.hi95ci"]<-unname(quantile(Out$CPUE,c(0.975)))
        IPHC.cpue[j,"CPUE.var"]<-var(Out$CPUE)
        IPHC.cpue[j,"CPUE.cv"]<-sd(Out$CPUE)/mean(Out$CPUE)
        #    IPHC.wcpue[j,"WCPUE32.cv"]<-sqrt(var(WCPUEi.32))/mean(WCPUEi.32)
        IPHC.cpue[j,"WCPUE.mean"]<-mean(WCPUEi)
        IPHC.cpue[j,"WCPUE.bootmean"]<-unname(quantile(Out$WCPUE,c(0.5)))  #mean WCPUE from Tribuzio - AWL: actually the median...
        IPHC.cpue[j,"WCPUE.lo95ci"]<-unname(quantile(Out$WCPUE,c(0.025)))
        IPHC.cpue[j,"WCPUE.hi95ci"]<-unname(quantile(Out$WCPUE,c(0.975)))
        IPHC.cpue[j,"WCPUE.var"]<-var(Out$WCPUE)
        IPHC.cpue[j,"WCPUE.cv"]<-sd(Out$WCPUE)/mean(Out$WCPUE)
        
        j<-j+1
      } else {}
    }
  }
  
  
  return(IPHC.cpue)
}


colnames(Survey)

# # This is used below to subset stations. All stations is used from 2024 assessment
# # and onward, as a CIE review recommended including 0 catch observations when standarziing CPUE
All.stations <- Station.sum$Station
# 
# # This is used to subset stations with YE catch. This was how the 
# # CPUE was calculated pre 2024. 
YE.stations<-Station.sum$Station[Station.sum$mean.YEcpue != 0]

#Survey_40p<-Survey %>% filter(Station %in% c(YE.stations_40p))
Survey_all <- Survey %>% filter(Station %in% c(All.stations) &
                                  SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

Survey_non0 <- Survey %>% filter(Station %in% c(YE.stations) &
                                   SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

Survey_non0_wronglon <- Survey.wronglon %>% filter(Station %in% c(YE.stations) &
                                                     SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

Survey_non0_2026 <- Survey_2026_correct %>% filter(Station %in% c(YE.stations_2026) &
                                                     SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

Survey_all_2026 <- Survey_2026_correct %>% filter(Station %in% c(All.stations) &
                                  SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))
# Exclude year × stratum combinations with insufficient sampling effort
# This was motivated by EYKT having only 5 stations sampled in 2025, well below the 
# historical average of 47 stations. All 5 were zero catch, leading the model to predict/ fit to 
# 0 CPUE. 
# 
# Threshold: < 25% of historical mean stations per stratum (mean = 1998–2023)
# Affects: 2024 EYKT (0 stations), 2025 NSEO (0 stations), 
#          2025 EYKT (5 stations, ~11% of historical mean ~47)
#          
# Applied to bootstrap, GAM, and sdmTMB indices

# historical_mean <- Survey_non0 %>%
#   filter(Year <= 2023) %>%
#   group_by(SEdist, Year) %>%
#   summarise(n = n_distinct(Station), .groups = "drop") %>%
#   group_by(SEdist) %>%
#   summarise(mean_stations = mean(n), .groups = "drop")
# 
# historical_mean
# 
# insufficient_combos <- Survey_non0 %>%
#   group_by(SEdist, Year) %>%
#   summarise(n = n_distinct(Station), .groups = "drop") %>%
#   left_join(historical_mean, by = "SEdist") %>%
#   filter(n < 0.25 * mean_stations) %>%
#   select(SEdist, Year)
# 
# insufficient_combos
# 
# # Apply filter to all survey datasets
# Survey_non0 <- Survey_non0 %>%
#   anti_join(insufficient_combos, by = c("SEdist", "Year"))
# 
# Survey_non0_wronglon <- Survey_non0_wronglon %>%
#   anti_join(insufficient_combos, by = c("SEdist", "Year"))
# 
# Survey_non0_2026 <- Survey_non0_2026 %>%
#   anti_join(insufficient_combos, by = c("SEdist", "Year"))
# 
# Survey_all <- Survey_all %>%
#   anti_join(insufficient_combos, by = c("SEdist", "Year"))
# 
# Survey_all_2026 <- Survey_all_2026 %>%
#   anti_join(insufficient_combos, by = c("SEdist", "Year"))


# For model = mod.0
IPHC.index.mod0 <- YEHA.fxn.noHA(Survey=Survey_non0, # If you want 0 station use survey_non0
                                 Area="SEdist",
                                 Deep=250, 
                                 Shallow=0,
                                 nboot=1000)

# For model = mod.1 and mod.3
IPHC.index.mod1_3 <- YEHA.fxn(Survey=Survey_non0, # If you want 0 station use survey_non0
                              Area="SEdist",
                              Deep=250, 
                              Shallow=0,
                              nboot=1000)

# For model = mod.2
IPHC.index.mod2 <- YEHA.fxn(Survey=Survey_non0_wronglon, # If you want 0 station use survey_non0
                            Area="SEdist",
                            Deep=250, 
                            Shallow=0,
                            nboot=1000)

# For model = mod.4
IPHC.index.mod4 <- YEHA.fxn(Survey=Survey_non0_2026, # If you want 0 station use survey_non0
                            Area="SEdist",
                            Deep=250, 
                            Shallow=0,
                            nboot=1000)


IPHC.index.mod0$name <- "22.2"
IPHC.index.mod1_3$name <- "22.2_hhokadj_kghook"
IPHC.index.mod2$name <- "22.2_lon"
IPHC.index.mod4$name <- "22.2 + correct stations"

IPHC.index <- rbind(IPHC.index.mod0,IPHC.index.mod1_3,IPHC.index.mod2,IPHC.index.mod4)

str(IPHC.index)

# Save for briding analysis
# write.csv(IPHC.index, file = here("outputs/bridging analysis/bridge.cpue.boot 16July2026.csv"))


# Compare boot WCPUE and CPUE indices by area
ggplot(IPHC.index, aes(x=Year)) +
  geom_point(aes(y=CPUE.mean),size=2) +
  geom_point(aes(y=WCPUE.mean),size=2,col="blue") +
  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci),
                col="black", alpha=0.5) +
  geom_errorbar(aes(ymin = WCPUE.lo95ci, ymax= WCPUE.hi95ci),
                col="blue", alpha=0.5) +
  facet_grid(name~mngmt.area) +
  xlab("\nYear") +
  ylab("Yelloweye CPUE / WCPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank()) +
  labs(title ="Yelloweye cpue in IPHC FISS",
       subtitle = "All stations")

# ggsave(paste0("outputs/YE_IPHC_cpue_non0_boot_16July", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")

#get average number of stations in each area 
IPHC.index %>% group_by(mngmt.area,name) %>%
  summarize(mean_no_stations = mean(no.stations),
            min_no_stations = min(no.stations),
            max_no_stations = max(no.stations))


#*********************************************************************************
# Model based CPUE estimate using Tweedie model
#*********************************************************************************
# # We will still restrict stations to those that have encountered YE at least once
# # in the time series and not consider stations below 250 fathoms
# Survey_non0<-Survey %>% filter(Station %in% c(YE.stations) &
#                                  SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))
# 
# # For 2024: we will use all stations < 250 fathoms
# Survey_all <- Survey %>% filter(Station %in% c(All.stations) &
#                                   SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))
# 
# #IPHC_tweed <- Survey_non0 %>% filter(AvgDepth.fm <= 250) %>% 
# IPHC_tweed <- Survey_all %>% filter(AvgDepth.fm <= 250) %>% 
#   group_by(Year,Station) %>%
#   mutate(CPUE = h.adj*YE.obs/HooksObserved,
#          WCPUE = CPUE*mean.YE.kg) %>%
#   select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,        
#          Purpose.Code, Date, Eff, 
#          Lat = BeginLat, Lon = BeginLon, Depth = AvgDepth.fm, 
#          Soak = Soak.time..min.., Temp.C, Pres = Max.Pressure..db., 
#          pH, Sal = Salinity.PSU, O2_1 = Oxygen_ml,
#          O2_2 = Oxygen_umol, Oxygen_sat,
#          CPUE, WCPUE) %>% 
#   mutate(Year = as.factor(Year),
#          SEdist = as.factor(SEdist),
#          Soak = as.numeric(Soak),
#          Temp.C = as.numeric(Temp.C),
#          Sal = as.numeric(Sal),
#          O2_1 = as.numeric(O2_1),
#          O2_2 = as.numeric(O2_2),
#          Oxygen_sat = as.numeric(Oxygen_sat)) %>%
#   unique() %>% data.frame()
# 
# head(data.frame(IPHC_tweed),30)
# colnames(IPHC_tweed)
# 
# str(IPHC_tweed)
# 
# IPHC_tweed %>% filter(Year == 2021 & Station == 3211)
# 
# nrow(IPHC_tweed %>% filter(is.na(Soak)))
# nrow(IPHC_tweed %>% filter(is.na(Depth)))
# nrow(IPHC_tweed %>% filter(is.na(Temp.C)))
# nrow(IPHC_tweed %>% filter(is.na(Pres)))
# nrow(IPHC_tweed %>% filter(is.na(pH)))
# nrow(IPHC_tweed %>% filter(is.na(Sal)))
# nrow(IPHC_tweed %>% filter(is.na(O2_1)))
# nrow(IPHC_tweed %>% filter(is.na(Lat)))
# nrow(IPHC_tweed %>% filter(is.na(Lon)))
# 
# ################################################################################
# ## Visually see how wcpue is related to variables:
# 
# # Depth
# ggplot(IPHC_tweed, aes(Depth, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # Soak
# ggplot(IPHC_tweed, aes(Soak, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # Temp
# hist(IPHC_tweed$Temp.C)
# ggplot(IPHC_tweed, aes(Temp.C, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # Pressure
# ggplot(IPHC_tweed, aes(Pres, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # pH
# ggplot(IPHC_tweed, aes(pH, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # Salinity
# ggplot(IPHC_tweed, aes(Sal, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # Oxygen 1
# ggplot(IPHC_tweed, aes(O2_1, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # Oxygen 2
# ggplot(IPHC_tweed, aes(O2_2, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # Oxygen saturation
# ggplot(IPHC_tweed, aes(Oxygen_sat, WCPUE)) + geom_point(shape = 20) + 
#   geom_smooth(size = 2, se = FALSE)
# 
# # LatLon
# ggplot(IPHC_tweed, aes(Lat, WCPUE)) +
#   geom_smooth(method = 'loess', span = 1, se = FALSE)
# 
# ggplot(IPHC_tweed, aes(Lon, WCPUE)) +
#   geom_smooth(method = 'loess', span = 1, se = FALSE)
# 
# ################################################################################
# ## Check for colinearity
# 
# IPHC_tweed %>% 
#   select(Depth, Soak, Lat, Lon) %>% 
#   GGally::ggpairs() 
# 
# #some strong correlations here...
# # almost everything correlated with depth! 
# # Temp, Pres pH and Salinity all correlated with each other... :( 
# 
# # probably best to just consider latlon, depth, and soak  colnames(IPHC_tweed)
# IPHC_tweed <- IPHC_tweed %>% 
#   select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,        
#          Purpose.Code, Date, Eff, 
#          Lat, Lon, Depth, 
#          Soak, 
#          CPUE, WCPUE)
# 
# 
# # Look at distribution of CPUE data
# ggplot(IPHC_tweed, aes(WCPUE)) + geom_density(alpha = 0.4, fill = 4)
# ggplot(IPHC_tweed, aes(log(WCPUE+0.01))) + geom_density(alpha = 0.4, fill = 4)
# ggplot(IPHC_tweed, aes(log(WCPUE+0.1*mean(IPHC_tweed$WCPUE,na.rm=T)))) + geom_density(alpha = 0.4, fill = 4)
# 
# #need complete data sets for running models: 
# nrow(IPHC_tweed)
# nrow(IPHC_tweed[complete.cases(IPHC_tweed),])
# 
# #only use complete cases... 
# fulldat<-IPHC_tweed[complete.cases(IPHC_tweed),]
# 
# m0 <- gam(WCPUE ~ Year * SEdist, data=fulldat, gamma=1.4, family=tw(), method = "REML")
# m.depth <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4), data=fulldat, gamma=1.4, family=tw(), method = "REML")
# 
# m.soak <- gam(WCPUE ~ Year * SEdist + s(Soak, k=4), data=fulldat, gamma=1.4, family=tw(), method = "REML")
# m.ll <- gam(WCPUE ~ Year * SEdist + te(Lon, Lat), data=fulldat, gamma=1.4, family=tw(), method = "REML")
# m.depth_soak <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + s(Soak, k=4), 
#                     data=fulldat, gamma=1.4, family=tw(), method = "REML")
# m.depth_ll <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + te(Lon, Lat), 
#                   data=fulldat, gamma=1.4, family=tw(), method = "REML")
# m.soak_ll <- gam(WCPUE ~ Year * SEdist + s(Soak, k=4) + te(Lon, Lat) ,
#                  data=fulldat, gamma=1.4, family=tw(), method = "REML")
# m.global <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + s(Soak, k=4) + te(Lon, Lat) , 
#                 data=fulldat, gamma=1.4, family=tw(), method = "REML")
# 
# model.list<-list(m0,#m0.trip,
#                  m.depth,m.soak,m.ll,
#                  m.depth_soak,m.depth_ll,m.soak_ll,m.global
# )
# 
# names(model.list)<-c("m0",#"trip",
#                      "depth","soak","latlon","depth+soak","depth+latlon",
#                      "soak+latlon","global")
# modsum0<-data.frame(); j<-1
# for (i in model.list) {
#   #mod<-i
#   modsum0[j,"model"]<-names(model.list[j])
#   modsum0[j,"aic"]<-AIC(i)
#   modsum0[j,"bic"]<-BIC(i)
#   modsum0[j,"qaic"]<-QAIC(i, chat= 1/(1-summary(i)$r.sq))
#   modsum0[j,"dev"]<-summary(i)$dev.expl
#   modsum0[j,"rsq"]<-summary(i)$r.sq
#   modsum0[j,"dev_exp"]<-summary(i)$dev.expl-summary(m0)$dev.expl
#   j<-j+1
# }
# 
# modsum0 %>% arrange(aic)
# modsum0 %>% arrange(qaic)
# modsum0 %>% arrange(bic)
# modsum0 %>% arrange(-dev)  
# modsum0 %>% arrange(-rsq) 
# 
# plot(m.global, page = 1, shade = TRUE, resid = TRUE, all = TRUE)
# summary(m.global)
# 
# # No residual patterns, but may be some outliers
# plot(fitted(m.global), resid(m.global))
# abline(h = 0, col = "red", lty = 2)
# 
# plot(m.global)
# 
# plot(m.depth_soak)
# plot(fitted(m.depth_soak), resid(m.depth_soak))
# 
# # Check for outliers
# which(fitted(m.global) < -1.5)   
# 
# gam.check(m.global)
# 
# mod_std_tweed <- m.global
# 
# # Cant have interaction between 
# 
# #cpue_dat<-cpue_nom %>% select(-c(Gear,Drifts,Stat, Depth))
# cpue_dat<-IPHC_tweed[complete.cases(IPHC_tweed),]
# cpue_dat <- as.data.frame(IPHC_tweed)
# 
# # Predictions ----
# 
# # Data set of average variables to predict CPUE from:
# std_dat_tweed <- expand.grid(Year = as.factor(unique(cpue_dat$Year)),
#                              SEdist = as.factor(unique(cpue_dat$SEdist)), #table(cpue_dat$Gear)
#                              Depth = mean(cpue_dat$Depth),
#                              Soak = mean(cpue_dat$Soak, na.rm=T), 
#                              #Lat = mean(cpue_dat$Lat),
#                              #Lon = mean(cpue_dat$Lon),
#                              dum = 1,
#                              dumstat = 1) %>%
#   mutate(Lat = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "EYKT"]),
#                          SEdist == "NSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "NSEO"]),
#                          SEdist == "CSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "CSEO"]),
#                          SEdist == "SSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "SSEO"])),
#          Lon = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "EYKT"]),
#                          SEdist == "NSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "NSEO"]),
#                          SEdist == "CSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "CSEO"]),
#                          SEdist == "SSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "SSEO"])))
# 
# # Predict CPUE
# pred_cpue <- predict(mod_std_tweed, std_dat_tweed, type = "link", se = TRUE)
# 
# #checking my code with Jane's... checks out :)
# # preds <- predict.gam(mod_std_tweed, type="response", std_dat_tweed, se = TRUE)
# 
# pred_cpue
# 
# #Put the standardized CPUE and SE into the data frame and convert to
# #backtransformed (bt) CPUE
# alpha <- 0.05  # for a 95% confidence interval on bycatch and discard estimates
# z <- qnorm(1 - alpha / 2)  # Z value for 95% CI
# 
# std_dat_tweed %>% 
#   mutate(fit = pred_cpue$fit,
#          se = pred_cpue$se.fit,
#          lower_link = fit - z * se,
#          upper_link = fit + z * se,
#          lower = mod_std_tweed$family$linkinv(lower_link),
#          upper = mod_std_tweed$family$linkinv(upper_link),
#          #upper = fit + (2 * se),
#          #lower = fit - (2 * se),
#          bt_cpue = exp(fit),
#          bt_upper = upper, #exp(upper),
#          bt_lower = lower, #exp(lower),
#          bt_se = (bt_upper - bt_cpue) / 2  ,
#          bt_cv = bt_se/bt_cpue
#   ) -> std_dat_tweed
# 
# 
# std_dat_tweed %>% 
#   select(Year, cpue = bt_cpue, upper = bt_upper, lower = bt_lower, SEdist) %>% 
#   mutate(CPUE = "Tweedie index") %>% data.frame() -> plot_dat
# 
# plot_dat %>% data.frame() %>%
#   ggplot(aes(group = 1)) + #geom_ribbon(aes(Year, ymin = lower, ymax = upper)) +
#   geom_ribbon(aes(Year, ymin = lower, ymax = upper), fill = wes_palette("Darjeeling1")[c(5)],
#               colour = NA, alpha = 0.5) +
#   geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2, show.legend = F) +
#   geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
#   facet_wrap(~ SEdist, scales = "free") +
#   scale_colour_manual(values =  wes_palette("Darjeeling1")[c(5,4)],
#                       aesthetics = c("colour","fill"), name = "IPHC CPUE") +
#   labs(x = "Year", y = "FISS CPUE (round kg / hook)\n") +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + #c(0.85, 0.9)) +
#   expand_limits(y = 0)
# 
# ggsave(paste0("output/IPHC_cpue_tweedie_cas_",YEAR,".png"), dpi=300, height=6, width=6, units="in")
# 
# # nominal CPUE?
# nom.CPUE.mean <- fulldat %>%
#   group_by(Year, SEdist) %>%
#   summarise(cpue = mean(WCPUE), no.stations = length(unique(Station))) %>%
#   mutate(CPUE = "nominal", upper = NA, lower = NA, cv = NA)
# 
# # Nominal CPUE for comparison ----
# 
# # IPHC_tweed %>%
# #   group_by(Year,SEdist) %>%
# #   do(data.frame(rbind(Hmisc::smean.cl.boot(.$WCPUE)))) %>%
# #   mutate(calc = "set lvl kg/hook",
# #          fsh_cpue = Mean,
# #          upper = Upper,
# #          lower = Lower,
# #          cv = (upper-lower)/1.96) -> fsh_sum_boot

# 2026 GAM for SAFE Bridging analysis--------------------------------------------
# The 2024 and 2025 IPHC surveys was very limited and didnt survey EYKT in 
# 2024 and NSEO in 2025. This GAM is fit using the 2024 selected model as a 
# bridging exercise for the 2026 SAFE. The GAM will not be used for future 
# assessments
# 
# 
IPHC_tweed <- Survey_all %>% filter(AvgDepth.fm <= 250) %>%
  group_by(Year,Station) %>%
  mutate(CPUE = h.adj*YE.obs/HooksObserved,
         WCPUE = CPUE*mean.YE.kg) %>%
  select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,
         Purpose.Code, Date, Eff,
         Lat = BeginLat, Lon = BeginLon, Depth = AvgDepth.fm,
         Soak = Soak.time..min.., Temp.C, Pres = Max.Pressure..db.,
         pH, Sal = Salinity.PSU, O2_1 = Oxygen_ml,
         O2_2 = Oxygen_umol, Oxygen_sat,
         CPUE, WCPUE) %>%
  mutate(Year = as.factor(Year),
         SEdist = as.factor(SEdist),
         Soak = as.numeric(Soak),
         Temp.C = as.numeric(Temp.C),
         Sal = as.numeric(Sal),
         O2_1 = as.numeric(O2_1),
         O2_2 = as.numeric(O2_2),
         Oxygen_sat = as.numeric(Oxygen_sat)) %>%
  unique() %>% data.frame()

IPHC_tweed <- IPHC_tweed %>%
  select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,
         Purpose.Code, Date, Eff,
         Lat, Lon, Depth,
         Soak,
         CPUE, WCPUE)


IPHC_tweed_2026 <- Survey_all_2026 %>% filter(AvgDepth.fm <= 250) %>%
  group_by(Year,Station) %>%
  mutate(CPUE = h.adj*YE.obs/HooksObserved,
         WCPUE = CPUE*mean.YE.kg) %>%
  select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,
         Purpose.Code, Date, Eff,
         Lat = BeginLat, Lon = BeginLon, Depth = AvgDepth.fm,
         Soak = Soak.time..min.., Temp.C, Pres = Max.Pressure..db.,
         pH, Sal = Salinity.PSU, O2_1 = Oxygen_ml,
         O2_2 = Oxygen_umol, Oxygen_sat,
         CPUE, WCPUE) %>%
  mutate(Year = as.factor(Year),
         SEdist = as.factor(SEdist),
         Soak = as.numeric(Soak),
         Temp.C = as.numeric(Temp.C),
         Sal = as.numeric(Sal),
         O2_1 = as.numeric(O2_1),
         O2_2 = as.numeric(O2_2),
         Oxygen_sat = as.numeric(Oxygen_sat)) %>%
  unique() %>% data.frame()

IPHC_tweed_2026 <- IPHC_tweed_2026 %>%
  select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,
         Purpose.Code, Date, Eff,
         Lat, Lon, Depth,
         Soak,
         CPUE, WCPUE)
# save for sdmTMB analysis
# write.csv(IPHC_tweed_2026, file = paste0(here(),"/outputs/IPHC_cpue_notstandardized_fixed_no25p_",YEAR,".csv"))


fulldat <- IPHC_tweed[complete.cases(IPHC_tweed),]
fulldat_2026 <- IPHC_tweed_2026[complete.cases(IPHC_tweed_2026),]

# Add combined year-stratum factor for strata specific indices
fulldat_bridge <- fulldat %>%
  mutate(
    Year_f_SEdist = droplevels(interaction(Year, SEdist, sep = "."))
  )

fulldat_bridge_2026 <- fulldat_2026 %>%
  mutate(
    Year_f_SEdist = droplevels(interaction(Year, SEdist, sep = "."))
  )

# Fit best model (global) with combined factor
m.global_bridge <- gam(
  WCPUE ~ Year_f_SEdist + s(Depth, k=4) + s(Soak, k=4) + te(Lon, Lat),
  data   = fulldat_bridge,
  gamma  = 1.4,
  family = tw(),
  method = "REML"
)

m.global_bridge_2026 <- gam(
  WCPUE ~ Year_f_SEdist + s(Depth, k=4) + s(Soak, k=4) + te(Lon, Lat),
  data   = fulldat_bridge_2026,
  gamma  = 1.4,
  family = tw(),
  method = "REML"
)

#
cpue_dat<-IPHC_tweed[complete.cases(IPHC_tweed),]
cpue_dat_2026<-IPHC_tweed_2026[complete.cases(IPHC_tweed_2026),]

# cpue_dat <- as.data.frame(IPHC_tweed)

# Prediction grid that excludes missing year/strata combos
std_dat_bridge <- expand.grid(
  Year   = as.factor(unique(fulldat_bridge$Year)),
  SEdist = as.factor(unique(fulldat_bridge$SEdist)),
  Depth  = mean(cpue_dat$Depth),
  Soak   = mean(cpue_dat$Soak, na.rm = TRUE),
  dum    = 1,
  dumstat = 1
) %>%
  mutate(
    Lat = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "EYKT"]),
                    SEdist == "NSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "NSEO"]),
                    SEdist == "CSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "CSEO"]),
                    SEdist == "SSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "SSEO"])),
    Lon = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "EYKT"]),
                    SEdist == "NSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "NSEO"]),
                    SEdist == "CSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "CSEO"]),
                    SEdist == "SSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "SSEO"])),
    Year_f_SEdist = factor(paste(Year, SEdist, sep = "."),
                           levels = levels(fulldat_bridge$Year_f_SEdist))
  ) %>%
  filter(!is.na(Year_f_SEdist))  # removes 2024 EYKT and 2025 NSEO

std_dat_bridge_2026 <- expand.grid(
  Year   = as.factor(unique(fulldat_bridge_2026$Year)),
  SEdist = as.factor(unique(fulldat_bridge_2026$SEdist)),
  Depth  = mean(cpue_dat_2026$Depth),
  Soak   = mean(cpue_dat_2026$Soak, na.rm = TRUE),
  dum    = 1,
  dumstat = 1
) %>%
  mutate(
    Lat = case_when(SEdist == "EYKT" ~ mean(cpue_dat_2026$Lat[cpue_dat_2026$SEdist == "EYKT"]),
                    SEdist == "NSEO" ~ mean(cpue_dat_2026$Lat[cpue_dat_2026$SEdist == "NSEO"]),
                    SEdist == "CSEO" ~ mean(cpue_dat_2026$Lat[cpue_dat_2026$SEdist == "CSEO"]),
                    SEdist == "SSEO" ~ mean(cpue_dat_2026$Lat[cpue_dat_2026$SEdist == "SSEO"])),
    Lon = case_when(SEdist == "EYKT" ~ mean(cpue_dat_2026$Lon[cpue_dat_2026$SEdist == "EYKT"]),
                    SEdist == "NSEO" ~ mean(cpue_dat_2026$Lon[cpue_dat_2026$SEdist == "NSEO"]),
                    SEdist == "CSEO" ~ mean(cpue_dat_2026$Lon[cpue_dat_2026$SEdist == "CSEO"]),
                    SEdist == "SSEO" ~ mean(cpue_dat_2026$Lon[cpue_dat_2026$SEdist == "SSEO"])),
    Year_f_SEdist = factor(paste(Year, SEdist, sep = "."),
                           levels = levels(fulldat_bridge_2026$Year_f_SEdist))
  ) %>%
  filter(!is.na(Year_f_SEdist))  # removes 2024 EYKT and 2025 NSEO

# Predict
pred_bridge <- predict(m.global_bridge, std_dat_bridge, type = "link", se = TRUE)
pred_bridge_2026 <- predict(m.global_bridge_2026, std_dat_bridge_2026, type = "link", se = TRUE)

alpha <- 0.05
z     <- qnorm(1 - alpha / 2)

std_dat_bridge <- std_dat_bridge %>%
  mutate(
    fit        = pred_bridge$fit,
    se         = pred_bridge$se.fit,
    lower_link = fit - z * se,
    upper_link = fit + z * se,
    lower      = m.global_bridge$family$linkinv(lower_link),
    upper      = m.global_bridge$family$linkinv(upper_link),
    bt_cpue    = exp(fit),
    bt_upper   = upper,
    bt_lower   = lower,
    bt_se      = (bt_upper - bt_cpue) / 2,
    bt_cv      = bt_se / bt_cpue
  )

std_dat_bridge_2026 <- std_dat_bridge_2026 %>%
  mutate(
    fit        = pred_bridge_2026$fit,
    se         = pred_bridge_2026$se.fit,
    lower_link = fit - z * se,
    upper_link = fit + z * se,
    lower      = m.global_bridge_2026$family$linkinv(lower_link),
    upper      = m.global_bridge_2026$family$linkinv(upper_link),
    bt_cpue    = exp(fit),
    bt_upper   = upper,
    bt_lower   = lower,
    bt_se      = (bt_upper - bt_cpue) / 2,
    bt_cv      = bt_se / bt_cpue
  )

# Export for bridging analysis
tweedie_bridge <- std_dat_bridge %>%
  select(Year, SEdist, cpue = bt_cpue, upper = bt_upper,
         lower = bt_lower, cv = bt_cv) %>%
  mutate(CPUE = "Tweedie index")

tweedie_bridge_2026 <- std_dat_bridge_2026 %>%
  select(Year, SEdist, cpue = bt_cpue, upper = bt_upper,
         lower = bt_lower, cv = bt_cv) %>%
  mutate(CPUE = "Tweedie index 2026")

tweedie_bridge2 <- rbind(tweedie_bridge,tweedie_bridge_2026)

# write.csv(tweedie_bridge2,
#           file = here("outputs/bridging analysis/IPHC_tweedie_bridge_16July2026.csv"))


plot_dat_bridge <- std_dat_bridge %>% 
  select(Year, cpue = bt_cpue, upper = bt_upper, lower = bt_lower, SEdist) %>% 
  mutate(CPUE = "Tweedie index") %>% data.frame() 

plot_dat_bridge %>% data.frame() %>%
  ggplot(aes(group = 1)) + #geom_ribbon(aes(Year, ymin = lower, ymax = upper)) +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper), fill = wes_palette("Darjeeling1")[c(5)],
              colour = NA, alpha = 0.5) +
  geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2, show.legend = F) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_wrap(~ SEdist, scales = "free") +
  scale_colour_manual(values =  wes_palette("Darjeeling1")[c(5,4)],
                      aesthetics = c("colour","fill"), name = "IPHC CPUE") +
  labs(x = "Year", y = "FISS CPUE (round kg / hook)\n") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + #c(0.85, 0.9)) +
  expand_limits(y = 0)


tweedie_bridge2 %>% 
  ggplot(aes(x = Year, y = cpue, col = CPUE, group = CPUE)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = CPUE), 
              alpha = 0.2, colour = NA) +
  facet_wrap(~SEdist) +
  theme_bw()

# ggsave(paste0("figures/IPHC_cpue_tweedie_bridge_allstations",YEAR,".png"), dpi=300, height=6, width=6, units="in")

#-------------------------------------------------------------------------------
# Compare Tweedie and bootstrap index if you want... 
boot_index <- IPHC.index %>% 
  mutate(SEdist = mngmt.area, fsh_cpue = WCPUE.mean, upper = WCPUE.hi95ci,
         lower = WCPUE.lo95ci, cv = (upper-lower)/1.96) %>%
  select(Year, SEdist, fsh_cpue, lower, upper, cv)

# Compare predicted cpue from gam to nominal cpue

cbpalette <- c("#009E73", "#0072B2", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7","#F0E442", "black", "grey")

index.compare <- boot_index %>%
  select(Year, cpue = fsh_cpue, upper, lower, SEdist) %>% 
  mutate(CPUE = "Raw index",
         Year = as.factor(Year)) %>% 
  bind_rows(std_dat_tweed %>% 
              select(Year, cpue = bt_cpue, upper = bt_upper, lower = bt_lower, SEdist) %>% 
              mutate(CPUE = "Standardized index")) %>% 
  mutate(Year = as.numeric(as.character(Year))) %>%
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = CPUE), 
              colour = NA, alpha = 0.3) +
  #geom_point(aes(Year, cpue, colour = CPUE), show.legend = F) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 0.5) +
  facet_wrap(~ SEdist, scales = "free_x") +
  scale_colour_manual(values =  cbpalette[c(1,2)],
                      aesthetics = c("colour","fill"), name = "") +
  labs(x = "Year", y = "CPUE (round kg / hook)\n") +
  theme(legend.position = "bottom") + #c(0.85, 0.9)) +
  expand_limits(y = 0) +
  theme_bw()

ggsave(paste0("REMA/Figures/",YEAR, "/IPHC_cpue_tweedie_boot_comp_cas_",YEAR,".png"), index.compare, dpi=300, height=5, width=6.5, units="in")


# compare with Rhea's index

rindex <- read.csv(paste0(here::here(), "/Data_processing/Data/IPHC.cpue.SEO_tweed_boot_rke_102824.csv"))

index.compare.r <- rindex %>%
  select(Year, cpue, CPUE, upper, lower, SEdist) %>% 
  mutate(CPUE = case_when(
    CPUE == "Bootstrap index" ~ "Raw index",
    CPUE == "Tweedie index" ~ "Standardized index"
  )) %>%
  mutate(Year = as.factor(Year)) %>%
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = CPUE, group = CPUE), 
              colour = NA, alpha = 0.3) +
  #geom_point(aes(Year, cpue, colour = CPUE), show.legend = F) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 0.5) +
  scale_x_discrete(breaks = seq(1998, 2025, by = 5)) +
  facet_wrap(~ SEdist, scales = "free_x") +
  scale_colour_manual(values =  cbpalette[c(1,2)],
                      aesthetics = c("colour","fill"), name = "") +
  labs(x = "Year", y = "CPUE (round kg / hook)\n") +
  theme(legend.position = "bottom") + #c(0.85, 0.9)) +
  expand_limits(y = 0) +
  theme_bw() + 
  theme(panel.spacing.x = unit(2, "lines"))

ggsave(paste0("REMA/Figures/",YEAR, "/IPHC_cpue_tweedie_boot_comp_rke_",YEAR,".png"), index.compare.r, dpi=300, height=5, width=6.5, units="in")

# Save these values to combine with other indices we'll generate: 
boot_index %>%
  select(Year, SEdist, cpue = fsh_cpue, upper, lower, cv) %>% 
  mutate(CPUE = "Bootstrap index",
         Year = as.factor(Year)) %>% 
  bind_rows(std_dat_bridge %>% 
              select(Year, SEdist, cpue = bt_cpue, 
                     upper = bt_upper, lower = bt_lower, cv = bt_cv) %>% 
              mutate(CPUE = "Tweedie index")) %>% 
  #mutate(Year = as.numeric(as.character(Year)))) %>% 
  data.frame() -> IPHC_cpue_indices

IPHC_cpue_indices<-left_join(IPHC_cpue_indices,IPHC.index %>% mutate(Year = as.factor(Year)) %>%
                               select(Year,SEdist = mngmt.area,no.stations),by=c("Year","SEdist"))

# SAVE all the indices for used in assessment models
write.csv(IPHC_cpue_indices,paste0("./Data_processing/Data/IPHC.cpue.SEO_tweed_boot_cas_",YEAR,".csv"))

#write.csv(IPHC_cpue_indices,paste0("./Data_processing/Data/IPHC.cpue.SEO_boot_non0_cas_",YEAR,".csv"))

# compare Tweedie index with nominal CPUE

tweedie.nom.compare <- IPHC_cpue_indices %>%
  rbind(nom.CPUE.mean)

ggplot(tweedie.nom.compare) +
  geom_point(aes(x = Year, y = cpue, color = CPUE)) +
  facet_wrap(~SEdist, ncol = 2)




tweedie.index <- std_dat_tweed %>% 
  select(Year, SEdist, cpue = bt_cpue, 
         upper = bt_upper, lower = bt_lower, cv = bt_cv) %>% 
  mutate(CPUE = "Tweedie index")


write.csv(tweedie.index, file = here("outputs/bridging analysis/IPHC_tweedie_10June2026.csv"))
