#*********************************************************************************
# Author: Aaron Lambert
# Date:   6/2/2026
# 
# Version: V3
#         This version of the sdmTMB standardization is for running the models without
#         strata (management district) as yearx strata factor.
#         This would simplify the models (less parameters 26 year coeficients vs 109
#          yearxstrata coeficients) and adress the SSC comment about combining strata for 
#          NSEO and adressing the reduced IPHC sampling...
#          
#          
# sdmTMB Spatio-Temporal CPUE Standardization
# IPHC Longline Survey – SEO Yelloweye Rockfish 
# 
# Purpose: Replaces the GAM (Tweedie) approach used in the 2024 DSR SAFE.The sdmTMB 
#          spatio-temporal index standardization was requested by the NPFMC SSC
#          in 2024. This approach is necessary because the IHPC is reducing sampling
#          due to budget constraints.
#
#*********************************************************************************

# Year that the assessment is run
Year <- 2026

library(sdmTMB)
library(sdmTMBextra)   
library(INLA)
library(tidyverse)
library(here)
library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(marmap)
library(RANN)
library(DHARMa)
library(gstat)
library(ggthemes)

cbpalette <- c("#009E73", "#0072B2", "#E69F00", "#56B4E9",
               "#D55E00", "#CC79A7", "#F0E442", "black", "grey")


# Load in CPUE (this was processed in the script IPHC_Sruvey_CPUE_index_cas.R)
# dat_raw <- read_csv(here("outputs",
#                          "IPHC_cpue_notstandardized_awl27May2026.csv"))
# dat_raw <- read_csv(here("outputs",
#                          "IPHC_cpue_notstandardized_2026.csv"))

# This one has the correct stations witht the 2025 EYKT ommitted
# dat_raw <- read_csv(here("outputs",
#                          "IPHC_cpue_notstandardized_fixed_2026.csv"))

# This one has the correct stations witht all stations included
dat_raw <- read_csv(here("outputs",
                         "IPHC_cpue_notstandardized_fixed_no25p_2026.csv"))

# Filter to get observations from less than 250 fathoms (same as 2024 GAM)
dat <- dat_raw %>%
  filter(Depth <= 250) %>%
  mutate(
    Year = as.integer(Year),
    Year_f = factor(Year),
    SEdist = factor(SEdist, levels = c("EYKT", "NSEO", "CSEO", "SSEO"))
  ) %>%
  drop_na(WCPUE, Depth, Soak, Lon, Lat)

# # Add UTM columns in km (needed for sdmTMB model)
dat <- add_utm_columns(dat, ll_names = c("Lon", "Lat"),
                       utm_crs = 32608,   # UTM zone 8N — SE Alaska
                       units   = "km")

# Take a look at number of observations, years, and how many zero catch
nrow(dat)
range(dat$Year)
sum(dat$WCPUE == 0)

# Build SPDE Mesh with Land Barrier 

# Initial mesh over observations
mesh <- make_mesh(dat,
                  xy_cols = c("X", "Y"),
                  cutoff  = 15)    # in km

# Number of initial mesh nodes
mesh$mesh$n

# Get Alaska coastline polygon, transform to UTM zone 8N
ak_sf <- ne_states(country = "united states of america", returnclass = "sf") %>%
  filter(name == "Alaska") %>%
  st_transform(32608) %>%
  st_union() %>%
  st_sf()

# Scale polygon from meters to km to match mesh units
ak_km <- ak_sf %>%
  mutate(geometry = geometry * 0.001) %>%
  st_set_crs(NA) %>%       # drop CRS (units are now km, not m)
  st_set_crs(32608)        # reattach so sdmTMBextra recognizes the object

# Add barrier to mesh
mesh_barrier <- sdmTMBextra::add_barrier_mesh(
  spde_obj = mesh,
  barrier_sf = ak_km,
  range_fraction = 0.1,    # land triangles get 10% of ocean range
  plot = TRUE   
)


# Plot barrier mesh map
# Extract triangle edges from the mesh and convert back to lon/lat for plotting
idx <- mesh$mesh$graph$tv        # triangle vertex indices (n_triangles × 3)
loc <- mesh$mesh$loc[, 1:2] * 1000  # km → m for CRS conversion

# Build edge endpoint dataframe (3 edges per triangle)
edges <- rbind(
  data.frame(x0 = loc[idx[,1],1], y0 = loc[idx[,1],2],
             x1 = loc[idx[,2],1], y1 = loc[idx[,2],2]),
  data.frame(x0 = loc[idx[,2],1], y0 = loc[idx[,2],2],
             x1 = loc[idx[,3],1], y1 = loc[idx[,3],2]),
  data.frame(x0 = loc[idx[,3],1], y0 = loc[idx[,3],2],
             x1 = loc[idx[,1],1], y1 = loc[idx[,1],2])
)

# function to convert UTM (m) coordinate pairs to lon/lat
to_ll <- function(x, y) {
  st_sfc(mapply(function(a, b) st_point(c(a, b)), x, y,
                SIMPLIFY = FALSE), crs = 32608) %>%
    st_transform(4326) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    setNames(c("lon", "lat"))
}

# Convert lat/lon
edges_ll <- bind_cols(
  to_ll(edges$x0, edges$y0) %>% rename(lon0 = lon, lat0 = lat),
  to_ll(edges$x1, edges$y1) %>% rename(lon1 = lon, lat1 = lat)
)

# Identify barrier (land) triangles from mesh_barrier object
barrier_idx <- mesh_barrier$barrier_triangles
barrier_edges <- rbind(
  data.frame(x0 = loc[idx[barrier_idx,1],1], y0 = loc[idx[barrier_idx,1],2],
             x1 = loc[idx[barrier_idx,2],1], y1 = loc[idx[barrier_idx,2],2]),
  data.frame(x0 = loc[idx[barrier_idx,2],1], y0 = loc[idx[barrier_idx,2],2],
             x1 = loc[idx[barrier_idx,3],1], y1 = loc[idx[barrier_idx,3],2]),
  data.frame(x0 = loc[idx[barrier_idx,3],1], y0 = loc[idx[barrier_idx,3],2],
             x1 = loc[idx[barrier_idx,1],1], y1 = loc[idx[barrier_idx,1],2])
)

barrier_edges_ll <- bind_cols(
  to_ll(barrier_edges$x0, barrier_edges$y0) %>% rename(lon0 = lon, lat0 = lat),
  to_ll(barrier_edges$x1, barrier_edges$y1) %>% rename(lon1 = lon, lat1 = lat)
)

# Basemap
ak_map <- ne_states(country = "united states of america", returnclass = "sf") %>%
  filter(name == "Alaska")

# Use unique station locations instead of all observations
unique_stations <- dat %>%
  distinct(Station, .keep_all = TRUE) %>%
  select(Station, Lon, Lat, SEdist)

p_mesh <- ggplot() +
  geom_sf(data = ak_map, fill = "grey85", colour = "grey50", linewidth = 0.3) +
  geom_segment(data = edges_ll,
               aes(x = lon0, y = lat0, xend = lon1, yend = lat1),
               colour = "steelblue", linewidth = 0.2, alpha = 0.5) +
  geom_segment(data = barrier_edges_ll,
               aes(x = lon0, y = lat0, xend = lon1, yend = lat1),
               colour = "darkorange", linewidth = 0.3, alpha = 0.7) +
  geom_point(data = unique_stations,  #unique locations only
             aes(Lon, Lat),
             colour = "red", size = 1.5, alpha = 0.8) +
  coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
  labs(
    # title    = "SPDE barrier mesh",
    # subtitle = "Blue = ocean triangles,\nOrange = barrier (land) triangles,\nRed = IPHC stations",
    x = "Longitude", y = "Latitude") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

p_mesh

ggsave(here("figures", "sdmTMB","sdmTMB_barrier_mesh.png"),
       p_mesh, width = 8, height = 9, dpi = 300)

# Number of barrier mesh nodes
mesh$mesh$n

# Number of triangles
nrow(mesh$mesh$graph$tv)

# Look at how close the stations are to each other
# 
# Get the distinct stations
# Need to average because each year the same station is slightly off by
# small distances...
coords <- dat %>%
  group_by(Station) %>%
  summarise(
    X = mean(X),
    Y = mean(Y),
    .groups = "drop"
  )

# Use the dist() function to get euclidean distance between points
dists   <- as.matrix(dist(coords[, c("X", "Y")]))

# The diagnal is the distance from a point to itself (i.e, = 0)
# Change to NA so it isnt picked up as the min distance
diag(dists) <- NA

# Get the quantiles
nn_dist <- apply(dists, 1, min, na.rm = TRUE)

summary(nn_dist)

hist(nn_dist, breaks = 30,
     main = "Nearest-neighbor distances between stations",
     xlab = "Distance (km)")

# Fit Candidate Models -------------------------------------------------------

# Function to fit the model
fit_sdmtmb <- function(formula, 
                       spatial = "on", 
                       family = tweedie(link = "log"),
                       spatiotemporal = "iid", ...) {
  sdmTMB(
    formula = formula,
    data = dat,
    mesh = mesh_barrier,
    family = family,
    time = "Year",
    spatial = spatial,          
    spatiotemporal = spatiotemporal,
    control = sdmTMBcontrol(newton_loops = 1L),
    ...
  )
}


# Check for years with no sampling in specific stat areas. This will be an 
# issue for fitting the model if there are yearxstrata combinations that are NA 
# in the data.
unsampled <- expand_grid(
  Year   = sort(unique(dat$Year)),
  SEdist = levels(dat$SEdist)
) %>%
  left_join(
    dat %>% group_by(Year, SEdist) %>% summarise(n = n(), .groups = "drop"),
    by = c("Year", "SEdist")
  ) %>%
  filter(is.na(n)) %>%
  mutate(n = 0)

unsampled

# Stations sampled by region and year
sampled <- expand_grid(
  Year   = sort(unique(dat$Year)),
  SEdist = levels(dat$SEdist)
) %>%
  left_join(
    dat %>% group_by(Year, SEdist) %>% summarise(n = n(), .groups = "drop"),
    by = c("Year", "SEdist")
  ) %>% 
  pivot_wider(names_from  = SEdist,values_from = n)

print(sampled, n = 30)

# Get the mean and sd of depth and soak to z-score for model fitting
depth_mean <- mean(dat$Depth)
depth_sd <- sd(dat$Depth)

soak_mean <- mean(dat$Soak)
soak_sd <- sd(dat$Soak)

# Drop empty year × stratum combinations
# IPHC did not sample some districts in 2024, 2025, 
# and likely more in the future..
dat <- dat %>%
  mutate(
    Year_f_SEdist = droplevels(interaction(Year_f, SEdist)),
    Depth_s = (Depth - depth_mean)/depth_sd,
    Soak_s = (Soak - soak_mean)/soak_sd
  )

# Lets fit some models!
# M1: depth + soak + IID spatiotemporal  
m1 <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
                 priors = sdmTMBpriors(
                   matern_st = pc_matern(range_gt = 10, sigma_lt = 5)
                 ))
sanity(m1)

set.seed(42)
sims_m1 <- simulate(m1, nsim = 500, type = "mle-mvn")
fitted_vals_m1 <- predict(m1)$est

# Get the dharma residuals
dharma_res_m1 <- DHARMa::createDHARMa(
  simulatedResponse = sims_m1,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals_m1
)

# Save plot
png(here("figures", "sdmTMB", "dharma_residuals_m1.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(dharma_res_m1)
dev.off()

# M2: depth + soak + AR(1) spatiotemporal
m2 <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
                 spatiotemporal = "ar1",
                 priors = sdmTMBpriors(
                   matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))
sanity(m2)

set.seed(42)
sims_m2 <- simulate(m2, nsim = 500, type = "mle-mvn")
fitted_vals_m2 <- predict(m2)$est

# Get the dharma residuals
dharma_res_m2 <- DHARMa::createDHARMa(
  simulatedResponse = sims_m2,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals_m2
)

png(here("figures", "sdmTMB", "dharma_residuals_m2.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(dharma_res_m2)
dev.off()

# M3: depth + soak, spatial only (no spatiotemporal field)
m3 <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
                 spatiotemporal = "off",
                 priors = sdmTMBpriors(
                   matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))
sanity(m3)

set.seed(42)
sims_m3 <- simulate(m3, nsim = 500, type = "mle-mvn")
fitted_vals_m3 <- predict(m3)$est

# Get the dharma residuals
dharma_res_m3 <- DHARMa::createDHARMa(
  simulatedResponse = sims_m3,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals_m3
)

png(here("figures", "sdmTMB", "dharma_residuals_m3.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(dharma_res_m3)
dev.off()

# M4: depth only + IID spatiotemporal
m4 <- fit_sdmtmb(WCPUE ~ 0 + Year_f  + s(Depth_s, k = 4),
                 priors = sdmTMBpriors(
                   matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))
sanity(m4)

set.seed(42)
sims_m4 <- simulate(m4, nsim = 500, type = "mle-mvn")
fitted_vals_m4 <- predict(m4)$est

# Get the dharma residuals
dharma_res_m4 <- DHARMa::createDHARMa(
  simulatedResponse = sims_m4,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals_m4
)

png(here("figures", "sdmTMB", "dharma_residuals_m4.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(dharma_res_m4)
dev.off()

# M5: soak only + IID spatiotemporal
m5 <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Soak_s, k = 4),
                 priors = sdmTMBpriors(
                   matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))
sanity(m5)

set.seed(42)
sims_m5 <- simulate(m5, nsim = 500, type = "mle-mvn")
fitted_vals_m5 <- predict(m5)$est

# Get the dharma residuals
dharma_res_m5 <- DHARMa::createDHARMa(
  simulatedResponse = sims_m5,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals_m5
)

png(here("figures", "sdmTMB", "dharma_residuals_m5.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(dharma_res_m5)
dev.off()

# M6: no covariates + IID spatiotemporal
m6 <- fit_sdmtmb(WCPUE ~ 0 + Year_f,
                 spatiotemporal = "ar1",
                 priors = sdmTMBpriors(
                   matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))
sanity(m6)

set.seed(42)
sims_m6 <- simulate(m6, nsim = 500, type = "mle-mvn")
fitted_vals_m6 <- predict(m6)$est

# Get the dharma residuals
dharma_res_m6 <- DHARMa::createDHARMa(
  simulatedResponse = sims_m6,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals_m6
)

png(here("figures", "sdmTMB", "dharma_residuals_m6.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(dharma_res_m6)
dev.off()

# M7: m2 with soak time removed (doesnt look like soak contributes a lot to the models)
m7 <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4),
                 spatiotemporal = "ar1",
                 priors = sdmTMBpriors(
                   matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))

sanity(m7)

set.seed(42)
sims_m7 <- simulate(m7, nsim = 500, type = "mle-mvn")
fitted_vals_m7 <- predict(m7)$est

# Get the dharma residuals
dharma_res_m7 <- DHARMa::createDHARMa(
  simulatedResponse = sims_m7,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals_m7
)

png(here("figures", "sdmTMB", "dharma_residuals_m7.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(dharma_res_m7)
dev.off()

# dat <- dat %>% filter(!(Year == 2025 & SEdist == "EYKT"))

# Delta models do not run unless you remove the 2025 EYKT stations (all 0 catch)
# m8 <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
#                        spatiotemporal = "iid",
#                        family = delta_lognormal(),
#                        priors = sdmTMBpriors(
#                          matern_st = pc_matern(range_gt = 10, sigma_lt = 5),
#                          matern_s = pc_matern(range_gt = 10, sigma_lt = 5)))
# 
# sanity(m8)



# DHARMa for delta model
# set.seed(42)
# sims_m8 <- simulate(m8, nsim = 500, type = "mle-mvn")
# fitted_m8 <- predict(m8)$est
# 
# dharma_res_m8 <- DHARMa::createDHARMa(
#   simulatedResponse    = sims_m8,
#   observedResponse     = dat$WCPUE,
#   fittedPredictedResponse = fitted_m8
# )
# 
# png(here("figures", "sdmTMB", "dharma_residuals_m8.png"),
#     width = 10, height = 6, units = "in", res = 300)
# plot(dharma_res_m8)
# dev.off()


# 
# str(dat)
# # Delta-gamma depth only
# m8_nodist <- fit_sdmtmb_delta(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
#                               spatiotemporal = "iid",
#                               family = delta_gamma(),
#                               priors = sdmTMBpriors(
#                                 matern_st = pc_matern(range_gt = 10, sigma_lt = 5),
#                                 matern_s = pc_matern(range_gt = 10, sigma_lt = 5)))
# 
# sanity(m8_nodist)
# 
# 
# 
# m2_nodist <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
#                         spatiotemporal = "ar1",
#                         priors = sdmTMBpriors(
#                           matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))
# 
# sanity(m2_nodist)
# 
# m1_nodist <- fit_sdmtmb(WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
#                         spatiotemporal = "iid",
#                         priors = sdmTMBpriors(
#                           matern_st = pc_matern(range_gt = 10, sigma_lt = 5)))
# 
# sanity(m1_nodist)
# 
# AIC(m_8_nodist,m2_nodist)

# # m_no_spatial: depth + soak + AR(1) spatiotemporal
# m_no_spatial <- fit_sdmtmb(WCPUE ~ Year_f * SEdist + s(Depth, k = 4) + s(Soak, k = 4),
#                  spatiotemporal = "ar1", spatial = "off")

# m_rw <- fit_sdmtmb(
#   WCPUE ~ Year_f * SEdist + s(Depth, k = 4) + s(Soak, k = 4),
#   spatial        = "off",
#   spatiotemporal = "rw"    # explicitly random walk instead of AR(1) near boundary (~1,0.99)
# )
# sanity(m_rw)
# AIC(m2, m_no_spatial, m_rw)
# The rw has a AIC that is 350 higher....no go.

# Model comparison and diagnostics ---------------------------------------------
# AIC table
aic_tbl <- AIC(m1, m2, m3, m4, m5, m6, m7) %>%
  rownames_to_column("model") %>%
  mutate(
    description = c(
      "spatial + IID + depth + soak",
      "spatial + AR(1) + depth + soak",
      "spatial only + depth + soak",
      "spatial + IID + depth",
      "spatial + IID + soak",
      "spatial + AR(1), no covariates",
      "spatial + AR(1) + depth"
      # "spatial + IID + depth + soak; delta_gamma"
      # "no spatial + AR(1) + depth + soak"
    ),
    delta_AIC = AIC - min(AIC)
  ) %>%
  arrange(AIC)

aic_tbl

# save for the SAFE report
# write.csv(x = aic_tbl, file = paste0(here(),"/outputs/sdmtmb aic table 17July2026.csv"))


# Cross validation -------------------------------------------------------------

# Create spatial blocks
# Using leave-one-block-out CV — each management unit is a block

# CV dataset — exclude 2024 and 2025
dat_cv <- dat %>%
  filter(Year <= 2023) %>%
  mutate(
    Year_f = droplevels(factor(Year)),
    Year_f_SEdist = droplevels(interaction(Year_f, SEdist, sep = "."))
  )

# # Add spatial blocks
# mutate(
#   block_geo = cut(Y,
#                   breaks = quantile(Y, probs = seq(0, 1, length.out = 6)),
#                   include.lowest = TRUE,
#                   labels = FALSE),
#   block_dist = as.integer(SEdist)


# cv_block_plot <- dat_cv %>%
#   ggplot(aes(Lon, Lat, colour = factor(block_geo))) +
#   geom_point(size = 1.5, alpha = 0.7) +
#   scale_colour_manual(values = cbpalette, name = "Block") +
#   coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
#   labs(title = "Spatial CV blocks (latitude only)",
#        x = "Longitude", y = "Latitude") +
#   theme_bw()
# 
# ggsave(here("figures", "sdmTMB","sdmTMB_cv_blocks.png"),
#        cv_block_plot, width = 8, height = 9, dpi = 300)

# Blocks by management district
# table(dat_cv$block_dist)

# Blocks across space
# table(dat_cv$block_geo)
# nrow(dat_cv)
# range(dat_cv$Year)

expand_grid(
  Year_f = levels(dat_cv$Year_f),
  SEdist = levels(dat_cv$SEdist)
) %>%
  mutate(combo = paste(Year_f, SEdist, sep = ".")) %>%
  filter(!combo %in% levels(dat_cv$Year_f_SEdist)) %>%
  print()

# CV excluding NSEO due to insufficient station coverage for spatial blocking
# dat_cv_noNSEO <- dat_cv %>%
#   filter(SEdist != "NSEO") %>%
#   mutate(
#     Year_f_SEdist = droplevels(interaction(Year_f, SEdist, sep = ".")),
#     Depth_s = (Depth - depth_mean)/depth_sd,
#     Soak_s = (Soak - soak_mean)/soak_sd
#   )
# 
# # Rebuild mesh without NSEO
# mesh_cv_noNSEO <- make_mesh(dat_cv_noNSEO,
#                             xy_cols = c("X", "Y"),
#                             cutoff  = 15)
# 
# mesh_barrier_cv_noNSEO <- sdmTMBextra::add_barrier_mesh(
#   spde_obj       = mesh_cv_noNSEO,
#   barrier_sf     = ak_km,
#   range_fraction = 0.1,
#   plot           = FALSE
# )
# 
# CV excluding NSEO due to insufficient station coverage for spatial blocking
dat_cv <- dat_cv %>%
  # filter(SEdist != "NSEO") %>%
  mutate(
    Year_f_SEdist = droplevels(interaction(Year_f, SEdist, sep = ".")),
    Depth_s = (Depth - depth_mean)/depth_sd,
    Soak_s = (Soak - soak_mean)/soak_sd,
    # Recompute spatial blocks on the filtered (NSEO-excluded) data
    block_geo = cut(Y,
                    breaks = quantile(Y, probs = seq(0, 1, length.out = 6)),
                    include.lowest = TRUE,
                    labels = FALSE)
  )

# Check new block balance
# table(dat_cv_noNSEO$block_geo)
# table(dat_cv_noNSEO$block_geo, dat_cv_noNSEO$SEdist)


cv_block_plot <- dat_cv %>%
  ggplot(aes(Lon, Lat, colour = factor(block_geo))) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = cbpalette, name = "Block") +
  coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
  labs(title = "Spatial CV blocks (latitude only)",
       x = "Longitude", y = "Latitude") +
  theme_bw()

ggsave(here("figures", "sdmTMB","sdmTMB_cv_blocks.png"),
       cv_block_plot, width = 8, height = 9, dpi = 300)


# Rebuild mesh without NSEO
mesh_cv <- make_mesh(dat_cv,
                            xy_cols = c("X", "Y"),
                            cutoff  = 15)

mesh_barrier_cv <- sdmTMBextra::add_barrier_mesh(
  spde_obj = mesh_cv,
  barrier_sf = ak_km,
  range_fraction = 0.1,
  plot = FALSE
)


# Get a set of folds to compare across all the cv runs
set.seed(42)
random_folds <- sample(rep(1:5, length.out = nrow(dat_cv)))

cv_m1_random <- sdmTMB_cv(
  formula = WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
  data = dat_cv,
  mesh = mesh_barrier_cv,
  family = tweedie(link = "log"),
  time = "Year",
  spatial = "on",
  spatiotemporal = "iid",
  k_folds = 5,
  fold_ids = random_folds,
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
)

# set.seed(42)

cv_m2_random <- sdmTMB_cv(
  formula = WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
  data = dat_cv,
  mesh = mesh_barrier_cv,
  family = tweedie(link = "log"),
  time = "Year",
  spatial = "on",
  spatiotemporal = "ar1",
  k_folds = 5,
  fold_ids = random_folds,
  control = sdmTMBcontrol(
    newton_loops = 2L,
    start = list(ar1_phi = 0.5)
  ),
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
)


cv_m3_random <- sdmTMB_cv(
  formula = WCPUE ~ 0 + Year_f + s(Depth_s, k = 4) + s(Soak_s, k = 4),
  data = dat_cv,
  mesh = mesh_barrier_cv,
  family = tweedie(link = "log"),
  time = "Year",
  spatial = "on",
  spatiotemporal = "off",
  k_folds = 5,
  fold_ids = random_folds,
  control = sdmTMBcontrol(
    newton_loops = 2L,
    start = list(ar1_phi = 0.5)
  ),
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
)
# set.seed(42)

cv_m4_random <- sdmTMB_cv(
  formula = WCPUE ~ 0 + Year_f + s(Depth_s, k = 4),
  data = dat_cv,
  mesh = mesh_barrier_cv,
  family = tweedie(link = "log"),
  time = "Year",
  spatial = "on",
  spatiotemporal = "iid",
  k_folds = 5,
  fold_ids = random_folds,
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
)

cv_m5_random <- sdmTMB_cv(
  formula = WCPUE ~ 0 + Year_f +  s(Soak_s, k = 4),
  data = dat_cv,
  mesh = mesh_barrier_cv,
  family = tweedie(link = "log"),
  time = "Year",
  spatial = "on",
  spatiotemporal = "iid",
  k_folds = 5,
  fold_ids = random_folds,
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
)


cv_m6_random <- sdmTMB_cv(
  formula = WCPUE ~ 0 + Year_f,
  data = dat_cv,
  mesh = mesh_barrier_cv,
  family = tweedie(link = "log"),
  time = "Year",
  spatial = "on",
  spatiotemporal = "ar1",
  k_folds = 5,
  fold_ids = random_folds,
  control = sdmTMBcontrol(
    newton_loops = 2L,
    start = list(ar1_phi = 0.5)
  ),
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
)

set.seed(42)

cv_m7_random <- sdmTMB_cv(
  formula = WCPUE ~ 0 + Year_f + s(Depth_s, k = 4),
  data = dat_cv,
  mesh = mesh_barrier_cv,
  family = tweedie(link = "log"),
  time = "Year",
  spatial = "on",
  spatiotemporal = "ar1",
  k_folds = 5,
  fold_ids = random_folds
  # control        = sdmTMBcontrol(
  #   start = list(ar1_phi = 0.5)
  # ),
  # priors         = sdmTMBpriors(
  #   matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
)


rmse_fn <- function(cv_model, model_name) {
  cv_model$data %>%
    filter(!is.na(cv_predicted)) %>%
    summarise(
      model = model_name,
      RMSE = sqrt(mean((WCPUE - cv_predicted)^2)),
      MAE = mean(abs(WCPUE - cv_predicted)),
      n = n()
    )
}

mae_fn <- function(cv_model, model_name) {
  cv_model$data %>%
    filter(!is.na(cv_predicted)) %>%
    summarise(
      model = model_name,
      MAE = mean(abs(WCPUE - cv_predicted)),
      n = n()
    )
}


cv_summary <- data.frame(
  model      = c("M1 IID + depth + soak",
                 "M2 AR1 + depth + soak",
                 "M3 Spatial only + depth + soak",
                 "M4 IID + depth",
                 "M5 IID + soak",
                 "M6 AR1 + no covariates",
                 "M7 AR1 + depth"
                 # "spatial + IID + depth + soak; delta_gamma"
  ),
  sum_loglik = c(cv_m1_random$sum_loglik,
                 cv_m2_random$sum_loglik,
                 cv_m3_random$sum_loglik,
                 cv_m4_random$sum_loglik,
                 cv_m5_random$sum_loglik,
                 cv_m6_random$sum_loglik,
                 cv_m7_random$sum_loglik),
  RMSE       = c(rmse_fn(cv_m1_random, "M1")$RMSE,
                 rmse_fn(cv_m2_random, "M2")$RMSE,
                 rmse_fn(cv_m3_random, "M3")$RMSE,
                 rmse_fn(cv_m4_random, "M4")$RMSE,
                 rmse_fn(cv_m5_random, "M5")$RMSE,
                 rmse_fn(cv_m6_random, "M6")$RMSE,
                 rmse_fn(cv_m7_random, "M7")$RMSE),
  MAE        = c(rmse_fn(cv_m1_random, "M1")$MAE,
                 rmse_fn(cv_m2_random, "M2")$MAE,
                 rmse_fn(cv_m3_random, "M3")$MAE,
                 rmse_fn(cv_m4_random, "M4")$MAE,
                 rmse_fn(cv_m5_random, "M5")$MAE,
                 rmse_fn(cv_m6_random, "M6")$MAE,
                 rmse_fn(cv_m7_random, "M7")$MAE),
  AIC        = c(AIC(m1), 
                 AIC(m2), 
                 AIC(m3), 
                 AIC(m4),
                 AIC(m5),
                 AIC(m6),
                 AIC(m7))
) %>%
  mutate(
    delta_loglik = sum_loglik - max(sum_loglik),
    delta_RMSE   = RMSE - min(RMSE),
    delta_AIC    = AIC - min(AIC)
  ) %>%
  arrange(desc(sum_loglik))


cv_summary

# Save the cv summary
# write.csv(x = cv_summary,
#           file = paste0(here(),"/Outputs/model cv and aic summary 17July2026.csv"))


# Select best model
best_model <- get(aic_tbl$model[1])
# best_model <- m2

# Randomized quantile residuals
set.seed(42)
sims <- simulate(best_model, nsim = 500, type = "mle-mvn")
fitted_vals <- predict(best_model)$est

# Get the dharma residuals
dharma_res <- DHARMa::createDHARMa(
  simulatedResponse = sims,
  observedResponse = dat$WCPUE,
  fittedPredictedResponse = fitted_vals
)

plot(dharma_res)

# See which stations are being underpredicted
# Identify outlier observations
outlier_idx <- which(dharma_res$scaledResiduals > 0.99)

print(dat[outlier_idx, ] %>% 
        select(Year, Station, SEdist, WCPUE, Depth, Soak),n= 50)

# Get Moran I results to see if there is still spatial correaltion not
# absorbed by the spatial field
moran_results <- map_dfr(sort(unique(dat$Year)), function(yr) {
  yr_idx <- which(dat$Year == yr)
  
  if (length(yr_idx) < 4) return(NULL)
  
  dharma_yr <- DHARMa::createDHARMa(
    simulatedResponse = sims[yr_idx, ],
    observedResponse = dat$WCPUE[yr_idx],
    fittedPredictedResponse = fitted_vals[yr_idx]
  )
  
  test <- DHARMa::testSpatialAutocorrelation(
    dharma_yr,
    x    = dat$X[yr_idx],
    y    = dat$Y[yr_idx],
    plot = FALSE
  )
  
  tibble(
    Year      = yr,
    n         = length(yr_idx),
    statistic = test$statistic["observed"], 
    p_value   = test$p.value
  )
}
)

print(moran_results,n=30)

# Look for significant moran results
moran_results %>%
  filter(p_value < 0.05) 

# # Look at 2017 to see if there is a reason for the correlation (even if it is small)
# dat %>% 
#   filter(Year == 2017) %>%
#   group_by(SEdist) %>%
#   summarise(n = n(), mean_cpue = mean(WCPUE))

sample_size_tbl <- dat %>%
  group_by(Year, SEdist) %>%
  summarise(n = n(), .groups = "drop") %>%
  # filter(SEdist == "NSEO") %>%
  pivot_wider(values_from = n, names_from = SEdist) %>% 
  arrange(Year) %>%
  print(n = 30)

# write.csv(x = sample_size_tbl, file = here("Outputs/Sample size table corrected bounds 25June2026.csv"))


# Residuals vs. covariates
dat$resid <- residuals(best_model, type = "mle-mvn")

# Plots residuals
p_resid <- (
  ggplot(dat, aes(Depth, resid)) +
    geom_point(alpha = 0.3) + geom_smooth() +
    labs(title = "Residuals vs. depth", 
         x = "Depth (fathoms)", 
         y = "Residual")
) + (
  ggplot(dat, aes(Soak, resid)) +
    geom_point(alpha = 0.3) + geom_smooth() +
    labs(title = "Residuals vs. soak time", 
         x = "Soak time (hours)",
         y = "Residual")
) + (
  ggplot(dat, aes(Year, resid, group = Year)) +
    geom_boxplot() +
    labs(title = "Residuals by year",
         x = "Year",
         y = "Residual")
) + (
  ggplot(dat, aes(SEdist, resid)) +
    geom_boxplot() +
    labs(title = "Residuals by management unit", 
         x = "Management unit", 
         y = "Residual")
)
p_resid


ggsave(here("figures", "sdmTMB","sdmTMB_redids_plots.png"),
       p_resid, width = 8, height = 8, dpi = 300)

# Variogram 
# convert to sf with residuals already attached
resid_sf <- dat %>%
  st_as_sf(coords = c("X", "Y"), crs = 32608)

vario <- gstat::variogram(resid ~ 1, data = resid_sf)

# Plot the variogram
plot(vario,
     main = "Semivariogram of sdmTMB residuals",
     xlab = "Distance (km)",
     ylab = "Semivariance")


# Check estimated spatial range from sdmTMB (in km)
tidy(best_model, "ran_pars")

# Extract fixed effects (Year_f_SEdist coefficients)
fixed_effects <- tidy(best_model, effects = "fixed", conf.int = TRUE)

# Extract random effect parameters (spatial, spatiotemporal, etc.)
ran_pars <- tidy(best_model, effects = "ran_pars", conf.int = TRUE)

# Combine into one clean table
model_params <- bind_rows(
  fixed_effects %>% mutate(type = "Fixed effect"),
  ran_pars      %>% mutate(type = "Random effect / variance parameter")
)

# Save full parameter table
write.csv(model_params,
          here("outputs", "m2_parameter_estimates.csv"),
          row.names = FALSE)

# For the SAFE report — just the variance parameters (not all 110 fixed effects)
safe_params <- ran_pars %>%
  mutate(
    term = case_when(
      term == "range"  ~ "Spatial range (km)",
      term == "phi" ~ "Tweedie dispersion (phi)",
      term == "sigma_O" ~ "Spatial SD (sigma_O)",
      term == "sigma_E" ~ "Spatiotemporal SD (sigma_E)",
      term == "tweedie_p" ~ "Tweedie power (p)",
      term == "rho"  ~ "AR(1) correlation (rho)",
      term == "sd__s(Depth_s)"~ "Depth smooth SD",
      term == "sd__s(Soak_s)" ~ "Soak smooth SD",
      TRUE ~ term
    )
  ) %>%
  select(
    Parameter = term,
    Estimate = estimate,
    SE = std.error,
    Lower_95 = conf.low,
    Upper_95 = conf.high
  ) %>%
  mutate(across(c(Estimate, SE, Lower_95), ~round(., 3))) %>%
  mutate(Upper_95 = case_when(
    Upper_95 > 1000 ~ NA_real_,  # replace astronomically large values with NA
    TRUE            ~ round(Upper_95, 3)
  ))

safe_params

# write.csv(safe_params,
#           here("outputs", "m2_variance_parameters_SAFE_17July2026.csv"),
#           row.names = FALSE)

# Predict to grid -------------------------------------------------------------

# Get the management area shapfiles
seak_mgmt <- st_read(here("data", "shapefiles",
                          "SEAK_Groundfish_Mgmt_Areas_Pseudo.shp"))

# Check contents
print(seak_mgmt %>% st_drop_geometry() %>% select(Name, Code))

# Filter to SEO strata and transform to UTM zone 8N (meters)
seo_bounds_m <- seak_mgmt %>%
  filter(Code %in% c("EYKT", "NSEO", "CSEO", "SSEO")) %>%
  st_transform(32608)

# Also create km version for mesh 
seo_bounds_km <- seo_bounds_m %>%
  mutate(geometry = geometry / 1000) %>%
  st_set_crs(NA) %>%
  st_set_crs(32608)

# Get bounding box in km for grid construction
bbox <- st_bbox(seo_bounds_km)

# Plot to confirm
ggplot() +
  geom_sf(data = seo_bounds_m %>% st_transform(4326),
          aes(fill = Code), alpha = 0.3) +
  scale_fill_manual(values = cbpalette, name = "Management unit") +
  coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
  labs(title = "SEO management unit boundaries",
       x = "Longitude", y = "Latitude") +
  theme_bw()

# Add land for fun....
ggplot() +
  geom_sf(data = ak_sf %>% st_transform(4326),
          fill = "grey85", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = seo_bounds_m %>% st_transform(4326),
          aes(fill = Code), alpha = 0.3, colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = cbpalette, name = "Management \nunit") +
  coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
  labs(
    x = "Longitude", y = "Latitude") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

station_p <- ggplot() +
  geom_sf(data = ak_sf %>% st_transform(4326),
          fill = "grey85", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = seo_bounds_m %>% st_transform(4326),
          aes(fill = Code), alpha = 0.3, colour = "grey30", linewidth = 0.3) +
  geom_point(data = unique_stations,
             aes(x = Lon, y = Lat))+
  scale_fill_manual(values = cbpalette, name = "Management \nunit") +
  coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
  labs(
    x = "Longitude", y = "Latitude") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

ggsave(here("figures", "sdmTMB","iphc stations and seo districs.png"),
       station_p, width = 8, height = 9, dpi = 300)

# Download bathymetry data for prediction data
# ETOPO bathymetry at 1 arc-minute resolution
# Downloaded from NOAA via marmap package

# Download bathymetry (cached after first download)
bathy <- getNOAA.bathy(
  lon1       = -141, lon2 = -132,
  lat1       = 53,   lat2 = 61,
  resolution = 1     # 1 arc-minute
)

# Convert to dataframe and filter to ocean only
bathy_df <- as.xyz(bathy) %>%
  setNames(c("Lon", "Lat", "Depth_m")) %>%
  filter(Depth_m < 0) %>%                      # ocean only (negative = below sea level)
  mutate(Depth_fath = abs(Depth_m) * 0.546807) # meters to fathoms

nrow(bathy_df)

# Depth range in area
round(range(bathy_df$Depth_fath), 1)


# Make prediction grids
# Trying three approaches to present to plan team
# 1. The full management area
# 2. The area with depth less than 250 fathoms (x2; one at 2km resolution)
# 3. The IPHC survey stations (accounts for years with no observations)

# Function to convert grid UTM km to lon/lat for plotting
to_lonlat <- function(grid_df) {
  grid_df %>%
    mutate(X_m = X * 1000, Y_m = Y * 1000) %>%
    st_as_sf(coords = c("X_m", "Y_m"), crs = 32608) %>%
    st_transform(4326) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    setNames(c("Lon", "Lat")) %>%
    bind_cols(grid_df %>% select(SEdist, grid))
}

# Full management area 
# 15 km resolution
grid_points_15km <- expand.grid(
  X = seq(bbox["xmin"], bbox["xmax"], by = 15),
  Y = seq(bbox["ymin"], bbox["ymax"], by = 15)
)

# 2 km resolution
grid_points_2km <- expand.grid(
  X = seq(bbox["xmin"], bbox["xmax"], by = 2),
  Y = seq(bbox["ymin"], bbox["ymax"], by = 2)
)

# As sf object
grid_sf_15km <- grid_points_15km %>%
  st_as_sf(coords = c("X", "Y"), crs = 32608)

# As sf object
grid_sf_2km <- grid_points_2km %>%
  st_as_sf(coords = c("X", "Y"), crs = 32608)

# Get the full grid with land border
grid_1_full <- st_join(grid_sf_15km, seo_bounds_km,
                       join = st_within) %>%
  filter(!is.na(Code)) %>%
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(X, Y, SEdist = Code) %>%
  mutate(grid = "Full management area")

# Get the full grid with land border(2km resolution)
grid_1_full_2km <- st_join(grid_sf_2km, seo_bounds_km,
                           join = st_within) %>%
  filter(!is.na(Code)) %>%
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(X, Y, SEdist = Code) %>%
  mutate(grid = "Full management area 2km")

# Cells per district
table(grid_1_full$SEdist)
table(grid_1_full_2km$SEdist)

# Cells for whole grid
nrow(grid_1_full)
nrow(grid_1_full_2km)

# Convert grid to lon/lat for depth lookup
grid_1_ll <- grid_1_full %>%
  mutate(X_m = X * 1000, Y_m = Y * 1000) %>%
  st_as_sf(coords = c("X_m", "Y_m"), crs = 32608) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  setNames(c("Lon", "Lat")) %>%
  bind_cols(grid_1_full %>% select(X, Y, SEdist, grid))

grid_1_ll_2km <- grid_1_full_2km %>%
  mutate(X_m = X * 1000, Y_m = Y * 1000) %>%
  st_as_sf(coords = c("X_m", "Y_m"), crs = 32608) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  setNames(c("Lon", "Lat")) %>%
  bind_cols(grid_1_full_2km %>% select(X, Y, SEdist, grid))

# Nearest neighbor depth lookup
nn_depth <- RANN::nn2(
  data  = bathy_df[, c("Lon", "Lat")],
  query = grid_1_ll[, c("Lon", "Lat")],
  k     = 1
)

nn_depth_2km <- RANN::nn2(
  data  = bathy_df[, c("Lon", "Lat")],
  query = grid_1_ll_2km[, c("Lon", "Lat")],
  k     = 1
)

grid_1_ll$Depth_fath <- bathy_df$Depth_fath[nn_depth$nn.idx]
grid_1_ll_2km$Depth_fath <- bathy_df$Depth_fath[nn_depth_2km$nn.idx]

# # For each grid cell, find all bathymetry points within 7.5 km (half cell width)
# # and average their depths
# nn_all <- RANN::nn2(
#   data  = bathy_df[, c("Lon", "Lat")],
#   query = grid_1_ll[, c("Lon", "Lat")],
#   k     = 20    # find 20 nearest neighbors
# )
# 
# # Average depth across neighbors within reasonable distance
# grid_1_ll$Depth_fath <- apply(nn_all$nn.idx, 1, function(idx) {
#   mean(bathy_df$Depth_fath[idx])
# })


# Grid that is depth restricted (15 km, 14-250 fm)
# Apply depth filter
grid_2_depth <- grid_1_ll %>%
  filter(Depth_fath >= 14 & Depth_fath <= 250) %>%
  select(X, Y, SEdist, Depth_fath) %>%
  mutate(grid = "Depth restricted (14-250 fm)")

grid_2_depth_2km <- grid_1_ll_2km %>%
  filter(Depth_fath >= 14 & Depth_fath <= 250) %>%
  select(X, Y, SEdist, Depth_fath) %>%
  mutate(grid = "Depth restricted (14-250 fm)")

# Number of cells
table(grid_2_depth$SEdist)
nrow(grid_2_depth)
table(grid_2_depth_2km$SEdist)
nrow(grid_2_depth_2km)

# Grid of fixed historical survey stations
# The 136 unique IPHC station locations sampled at least once 1998-present
# Assigned to their most frequently observed management unit
# Boundary stations (n=7) assigned to most common stratum

fixed_stations <- dat %>%
  group_by(Station, SEdist) %>%
  summarise(
    X     = mean(X),
    Y     = mean(Y),
    Depth = mean(Depth),
    n     = n(),
    .groups = "drop"
  ) %>%
  group_by(Station) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-n)

grid_3_stations <- fixed_stations %>%
  select(X, Y, SEdist) %>%
  mutate(grid = "Fixed survey stations")

# Cells
print(table(grid_3_stations$SEdist))
nrow(grid_3_stations)

# Plot the three grids
all_grids <- bind_rows(
  to_lonlat(grid_1_full),
  to_lonlat(grid_2_depth %>% select(X, Y, SEdist, grid)),
  to_lonlat(grid_3_stations)
) %>%
  mutate(grid = factor(grid, levels = c(
    "Full management area",
    "Depth restricted (14-250 fm)",
    "Fixed survey stations"
  )))

p_grids <- ggplot() +
  geom_sf(data = ak_map, fill = "grey85", colour = "grey50", linewidth = 0.3) +
  geom_point(data = all_grids,
             aes(Lon, Lat, colour = SEdist),
             size = 1, alpha = 0.7) +
  scale_colour_manual(values = cbpalette, name = "Management unit") +
  facet_wrap(~grid, ncol = 3) +
  coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
  labs(title    = "Prediction grid comparison – SEO yelloweye rockfish",
       subtitle = "15 km resolution",
       x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "top",
        strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 90))

p_grids

ggsave(here("figures","sdmTMB", "prediction_grid_comparison.png"),
       p_grids, width = 14, height = 7, dpi = 300)

# Plot the 15 km prediction grid
ggplot() +
  geom_sf(data = ak_map, fill = "grey85", colour = "grey50", linewidth = 0.3) +
  geom_point(data = to_lonlat(grid_2_depth),
             aes(Lon, Lat, colour = SEdist),
             size = 1, alpha = 0.7) +
  scale_colour_manual(values = cbpalette, name = "Management unit") +
  # facet_wrap(~grid, ncol = 3) +
  coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
  labs(
    # title    = "Prediction grid comparison – SEO yelloweye rockfish",
    #    subtitle = "2 km resolution",
    x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "top",
        strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 90)) -> grid_15km_p

# ggsave(here("figures","sdmTMB", "prediction_grid_15km.png"),
#        grid_15km_p, width = 7, height = 7, dpi = 300)

# Grid 1: Full management area 
grid_1_full_depth <- grid_1_ll %>%   
  select(X, Y, SEdist, grid, Depth_fath) %>%
  rename(Depth = Depth_fath)

# Grid 2: Depth restricted 
grid_2_depth_d <- grid_2_depth %>%
  rename(Depth = Depth_fath)

grid_2_depth_d_2km <- grid_2_depth_2km %>%
  rename(Depth = Depth_fath)

# Grid 3: Fixed stations — use mean observed depth per station
grid_3_stations_depth <- fixed_stations %>%
  select(X, Y, SEdist, Depth) %>%    
  mutate(grid = "Fixed survey stations")


# Predict_index without semi_join and Year_f_SEdist
predict_index <- function(pred_grid_in, grid_name, best_model = best_model, se_fit = TRUE) {
  
  stopifnot(all(c("X", "Y", "SEdist", "Depth") %in% names(pred_grid_in)))
  
  pred_full <- expand_grid(
    Year = sort(unique(dat$Year))
  ) %>%
    cross_join(pred_grid_in %>% select(X, Y, SEdist, Depth)) %>%
    mutate(
      Year_f  = factor(Year, levels = levels(dat$Year_f)),
      Soak_s  = (median(dat$Soak) - soak_mean) / soak_sd,
      Depth_s = (Depth - depth_mean) / depth_sd
    )
  
  cat("\nRunning predictions for:", grid_name,
      "(", nrow(pred_full), "rows)\n")
  
  pred_out <- predict(best_model,
                      newdata = pred_full,
                      se_fit  = se_fit)
  
  index <- pred_out %>%
    mutate(cpue_pred = exp(est)) %>%
    group_by(Year, SEdist) %>%
    # group_by(Year) %>% # uncomment to get overall SEO trend
    summarise(
      mean_cpue = mean(cpue_pred),
      sd_cpue   = sd(cpue_pred),
      n_cells   = n(),
      .groups   = "drop"
    ) %>%
    mutate(
      se_cpue  = sd_cpue / sqrt(n_cells),
      lower_95 = mean_cpue - 1.96 * se_cpue,
      upper_95 = mean_cpue + 1.96 * se_cpue,
      CV       = se_cpue / mean_cpue,
      grid     = grid_name
    ) %>%
    select(Year, SEdist, mean_cpue, CV, lower_95, upper_95, n_cells, grid)
    # select(Year, mean_cpue, CV, lower_95, upper_95, n_cells, grid) # Uncomment to get SEO index
  
  
  return(list(index    = index,
              raw_pred = pred_out))
}


# Get predictions
result_full <- predict_index(grid_1_full_depth,
                             "Full management area",
                             se_fit = FALSE, best_model = m2)

result_depth  <- predict_index(grid_2_depth_d,
                               "Depth restricted (14-250 fm)",
                               se_fit = FALSE, best_model = m2)

result_depth_2km <- predict_index(grid_2_depth_d_2km,
                                  "Depth restricted (2km: 14-250 fm)",
                                  se_fit = FALSE, best_model = m2)

result_stations  <- predict_index(grid_3_stations_depth,
                                  "Fixed survey stations",
                                  se_fit = FALSE, best_model = m2)

# # Extract index for comparison
index_full <- result_full$index
index_depth <- result_depth$index
index_depth_2km<- result_depth_2km$index
index_stations <- result_stations$index

# # Save raw predictions for spatial maps
pred_raw_full <- result_full$raw_pred
pred_raw_depth <- result_depth$raw_pred
pred_raw_depth_2km <- result_depth_2km$raw_pred
pred_raw_stations <- result_stations$raw_pred



# Index for other models not selected
result_depth_mod1 <- predict_index(grid_2_depth_d,
                                       "Depth restricted (14-250 fm)",
                                       se_fit = FALSE, best_model = m1)

result_depth_mod3 <- predict_index(grid_2_depth_d,
                                   "Depth restricted (14-250 fm)",
                                   se_fit = FALSE, best_model = m3)

index_m1 <- result_depth_mod1$index

index_m3 <- result_depth_mod3$index
# result_depth_2km_mod2 <- predict_index(grid_2_depth_d_2km,
#                                        "Depth restricted (14-250 fm)",
#                                        se_fit = FALSE, best_model = m2)

index_m1$model <- "M.1"

index_m3$model <- "M.3"

index_m2 <- result_depth$index

index_m2$model <- "M.2"


cor(index_m1$mean_cpue, index_m2$mean_cpue)
cor(index_m2$mean_cpue, index_m3$mean_cpue)


comp_df <- rbind(index_m1, index_m2, index_m3)


# Correlation table between M1, M2, M3 by stratum
cor_tbl <- comp_df %>%
  select(Year, SEdist, mean_cpue, model) %>%
  pivot_wider(names_from = model, values_from = mean_cpue) %>%
  group_by(SEdist) %>%
  summarise(
    r_M1_M2 = round(cor(M.1, M.2, use = "complete.obs"), 2),
    r_M1_M3 = round(cor(M.1, M.3, use = "complete.obs"), 2),
    r_M2_M3 = round(cor(M.2, M.3, use = "complete.obs"), 2),
    .groups = "drop"
  ) %>%
  bind_rows(
    # Overall correlation across all strata
    comp_df %>%
      select(Year, SEdist, mean_cpue, model) %>%
      pivot_wider(names_from = model, values_from = mean_cpue) %>%
      summarise(
        SEdist  = "All strata",
        r_M1_M2 = round(cor(M.1, M.2, use = "complete.obs"), 2),
        r_M1_M3 = round(cor(M.1, M.3, use = "complete.obs"), 2),
        r_M2_M3 = round(cor(M.2, M.3, use = "complete.obs"), 2)
      )
  )



# write.csv(comp_df,
#           here("outputs",
#                "bridging analysis",
#                "IPHC.cpue.SEO_sdmTMB_m1_m2_m3_17July2026.csv"),
#           row.names = FALSE)

index_comp_plot <- comp_df %>% 
  ggplot(
    aes(Year, 
        mean_cpue,
        colour = model, 
        linetype = model,
        fill = model)) +
  # geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
  #             alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_wrap(~SEdist, scales = "free_y") +
  scale_fill_colorblind()+
  scale_color_colorblind()+
  labs(
    # title = "Comparison of M1, M2, and M3",
    # subtitle = "M1 = IID, M2 = AR(1) Tweedie, M3 = spatial only",
    x = "Year", y = "Standardized CPUE (kg/hook)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

index_comp_plot

# ggsave(here("figures","sdmTMB", "m1_m2_m3_comp_plot_SAFE.png"),
#        index_comp_plot, width = 7, height = 7, dpi = 300)

# Combine all indices 
index_compare <- bind_rows(
  index_full,
  index_depth,
  index_depth_2km,
  index_stations
) %>%
  mutate(grid = factor(grid, levels = c(
    "Full management area",
    "Depth restricted (14-250 fm)",
    "Depth restricted (2km: 14-250 fm)",
    "Fixed survey stations"
  ))) 



# Plot index comparison
p_index_compare <- ggplot(index_compare,
                          aes(Year, 
                              mean_cpue,
                              colour = grid, 
                              linetype = grid,
                              fill = grid)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_wrap(~SEdist, scales = "free_y") +
  scale_colour_manual(
    values = c("Full management area"              = "#E69F00",
               "Depth restricted (14-250 fm)"      = "#0072B2",
               "Depth restricted (2km: 14-250 fm)" = "black",
               "Fixed survey stations"             = "#009E73"),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Full management area"              = "#E69F00",
               "Depth restricted (14-250 fm)"      = "#0072B2",
               "Depth restricted (2km: 14-250 fm)" = "black",
               "Fixed survey stations"             = "#009E73"),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Full management area"              = "dotted",
               "Depth restricted (14-250 fm)"      = "dashed",
               "Depth restricted (2km: 14-250 fm)" = "dotdash",
               "Fixed survey stations"             = "solid"),
    name = NULL
  ) +
  labs(title    = "sdmTMB CPUE index – prediction grid sensitivity",
       subtitle = "IPHC SEO yelloweye rockfish",
       x = "Year", y = "Standardized CPUE (kg/hook)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom")

p_index_compare

ggsave(paste0(here(),"/figures/sdmTMB/sdmTMB_index_grid_sensitivity_delta",Year,".png"),
       p_index_compare, width = 10, height = 7, dpi = 300)

# Correlation Between Grid Approaches 
index_wide <- index_compare %>%
  select(Year, SEdist, mean_cpue, grid) %>%
  pivot_wider(names_from  = grid,
              values_from = mean_cpue,
              names_repair = "universal")

cor_tbl <- index_wide %>%
  group_by(SEdist) %>%
  summarise(
    cor_full_vs_stations  = cor(`Full.management.area`,
                                `Fixed.survey.stations`,
                                method = "pearson"),
    cor_depth_vs_stations = cor(`Depth.restricted..14.250.fm.`,
                                `Fixed.survey.stations`,
                                method = "pearson"),
    cor_depth_2km_vs_stations = cor(`Depth.restricted..2km..14.250.fm.`,
                                    `Fixed.survey.stations`,
                                    method = "pearson"),
    cor_depth_2km_vs_15km = cor(`Depth.restricted..2km..14.250.fm.`,
                                `Depth.restricted..14.250.fm.`,
                                    method = "pearson"),
    
    .groups = "drop"
  )

# Correlations with fixed station approach
cor_tbl

# Export the four indices for REMA model
rema_all <- index_compare %>%
  mutate(method = grid) %>%
  select(method, Year, SEdist, CPUE = mean_cpue, CV,
         lower = lower_95, upper = upper_95)

write.csv(rema_all,
          here("outputs", "bridging analysis",
               "IPHC.cpue.SEO_sdmTMB_grid_comparison_delta_15July2026.csv"),
          row.names = FALSE)

# Compare sdmtmb to gam index
gam_index <- read.csv(here("outputs", "bridging analysis",
                           "IPHC_tweedie_bridge_16July2026.csv")) %>%
  filter(CPUE == "Tweedie index") %>%
  select(Year, SEdist, cpue, cv) %>%
  mutate(method = "GAM (Tweedie)")

compare <- bind_rows(
  index_compare %>%
    select(Year, SEdist, cpue = mean_cpue, cv = CV, type = grid) %>%
    mutate(method = "sdmTMB (Tweedie)"),
  gam_index %>% mutate(type = "GAM")
)

# Look at correlations with GAM
all_index_wide <- compare %>%
  select(Year, SEdist, cpue, type) %>%
  filter(Year<=2023) %>% 
  pivot_wider(names_from  = type,
              values_from = cpue,
              names_repair = "universal")

cor_tbl_gam <- all_index_wide %>%
  group_by(SEdist) %>%
  summarise(
    cor_GAM_vs_Full  = cor(`Full.management.area`,
                           `GAM`,
                           method = "pearson"),
    cor_GAM_vs_depth = cor(`Depth.restricted..2km..14.250.fm.`,
                           `GAM`,
                           method = "pearson"),
    cor_GAM_vs_stations = cor(`Fixed.survey.stations`,
                              `GAM`,
                              method = "pearson"),
    .groups = "drop"
  )

cor_tbl_gam

p_compare <- compare %>% 
  ggplot(aes(Year, cpue, colour = type, linetype = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_wrap(~SEdist, scales = "free_y") +
  # scale_colour_manual(values = c("GAM (Tweedie)"    = "#009E73",
  # "sdmTMB (Tweedie)" = "#0072B2")) +
  scale_color_manual(values = cbpalette)+
  labs(title    = "GAM vs. sdmTMB standardized CPUE",
       subtitle = "IPHC SEO yelloweye rockfish",
       x = "Year", y = "Standardized CPUE (kg/hook)",
       colour = NULL, linetype = NULL) +
  theme_bw(base_size = 16) +
  theme(legend.position = "top")

# p_compare <- compare %>% 
#   filter(type %in% c("GAM", "Fixed survey stations")) %>% 
#   ggplot(
#     aes(Year, cpue, colour = type, linetype = type)) +
#   geom_line(linewidth = 0.7) +
#   geom_point(size = 1.2) +
#   facet_wrap(~SEdist, scales = "free_y") +
#   # scale_colour_manual(values = c("GAM (Tweedie)"    = "#009E73",
#   # "sdmTMB (Tweedie)" = "#0072B2")) +
#   scale_color_colorblind()+
#   labs(title    = "GAM vs. sdmTMB standardized CPUE",
#        subtitle = "IPHC SEO yelloweye rockfish",
#        x = "Year", y = "Standardized CPUE (kg/hook)",
#        colour = NULL, linetype = NULL) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "bottom")

p_compare

# ggsave(here("figures", "sdmTMB","sdmTMB_vs_GAM_IPHC_cpue.png"),
#        p_compare, width = 14, height = 7, dpi = 300)


# Plot comparing 2km and GAM for DSR appendix
p22_compare <- compare %>%
  filter(type %in% c("GAM","Depth restricted (2km: 14-250 fm)")) %>% 
  mutate(type = recode(type,
                       "GAM" = "GAM (Tweedie)",
                       "Depth restricted (2km: 14-250 fm)" = "Spatiotemporal"
  )) %>% 
  ggplot(aes(Year, cpue, colour = type, linetype = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_wrap(~SEdist, scales = "free_y") +
  # scale_colour_manual(values = c("GAM (Tweedie)"    = "#009E73",
  # "sdmTMB (Tweedie)" = "#0072B2")) +
  scale_color_manual(values = cbpalette[c(3,4)])+
  labs(x = "Year", y = "Standardized CPUE (kg/hook)",
       colour = NULL, linetype = NULL) +
  scale_x_continuous(breaks = seq(1995,2025,5))+
  theme_bw(base_size = 18) +
  theme(legend.position = "top")

p22_compare

ggsave(here("figures", "sdmTMB","sdmTMB_2km_vs_GAM_IPHC_cpue.png"),
       p22_compare,  dpi = 300)


# Best model (M2) index plots by strata
result_depth_M2  <- predict_index(grid_2_depth_d,
                                   "Depth restricted (14-250 fm)",
                                   se_fit = FALSE, best_model = m2)

M2_index <- result_depth_M2$index

M2_index %>% 
  ggplot(aes(Year, mean_cpue)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.5)+
  facet_wrap(~SEdist, scales = "free_y") +
  # scale_colour_manual(values = c("GAM (Tweedie)"    = "#009E73",
  # "sdmTMB (Tweedie)" = "#0072B2")) +
  scale_color_manual(values = cbpalette[c(3,4)])+
  labs(x = "Year", y = "Standardized CPUE (kg/hook)",
       colour = NULL, linetype = NULL) +
  scale_x_continuous(breaks = seq(1995,2025,5))+
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = .9)) -> M2_index_p


M2_index_p

# ggsave(here("figures","sdmTMB", "M2_index_SAFE.png"),
#        M2_index_p, width = 7, height = 7, dpi = 300)


# Spatial Field Plots ------------------------------------------------------ 
# Visualise the estimated spatial (omega) and spatiotemporal (epsilon) random
# fields for the most recent year.

plot.years <- sort(unique(dat$Year))

pdf(file = here("figures", "sdmTMB", "spatialfield_plots.pdf"))
for (y in plot.years) {
  
  # y <- 2020
  
  YEAR <- y
  
  p_omega <- pred_raw_depth_2km %>%
    filter(Year == YEAR) %>%
    mutate(
      X_m = X * 1000,
      Y_m = Y * 1000
    ) %>%
    st_as_sf(coords = c("X_m", "Y_m"), crs = 32608) %>%
    st_transform(4326) %>%
    mutate(
      Lon = st_coordinates(.)[, 1],
      Lat = st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry() %>%
    ggplot() +
    geom_sf(data = ak_map, fill = "grey85", colour = "grey50", linewidth = 0.3) +
    geom_point(aes(Lon, Lat, colour = omega_s), size = 3) +
    scale_colour_gradient2(low  = "blue", 
                           mid  = "white", 
                           high = "red",
                           midpoint = 0, 
                           name = "ω") +
    coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
    labs(title = paste("Persistent spatial field (omega) –", YEAR),
         x = "Longitude", y = "Latitude") +
    theme_bw()
  
  # Plot epsilon (spatiotemporal field) for most recent year
  p_eps <- pred_raw_depth_2km %>%
    filter(Year == YEAR) %>%
    mutate(
      X_m = X * 1000,
      Y_m = Y * 1000
    ) %>%
    st_as_sf(coords = c("X_m", "Y_m"), crs = 32608) %>%
    st_transform(4326) %>%
    mutate(
      Lon = st_coordinates(.)[, 1],
      Lat = st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry() %>%
    ggplot() +
    geom_sf(data = ak_map, fill = "grey85", colour = "grey50", linewidth = 0.3) +
    geom_point(aes(Lon, Lat, colour = epsilon_st), size = 3) +
    scale_colour_gradient2(low  = "blue",
                           mid  = "white",
                           high = "red",
                           midpoint = 0,
                           name = "ε") +
    coord_sf(xlim = c(-140, -133), ylim = c(54, 61)) +
    labs(title = paste("Spatiotemporal field (epsilon) –", YEAR),
         x = "Longitude", y = "Latitude") +
    theme_bw()
  
  print(p_omega + p_eps) 
  
  
  
}

dev.off()



# Plot the predicted CPUE as densitys on map (Spatial maps)---------------------
# Function to add lat and lon for plotting
add_lonlat <- function(pred_df, grid_name) {
  pred_df %>%
    mutate(X_m = X * 1000, Y_m = Y * 1000) %>%
    st_as_sf(coords = c("X_m", "Y_m"), crs = 32608) %>%
    st_transform(4326) %>%
    mutate(
      Lon  = st_coordinates(.)[, 1],
      Lat  = st_coordinates(.)[, 2],
      cpue = exp(est),
      grid = grid_name
    ) %>%
    st_drop_geometry() %>%
    select(Year, SEdist, Lon, Lat, cpue, est, grid)
}

# Add lat and lon to prediction outputs
spatial_full <- add_lonlat(pred_raw_full,     "Full management area")
spatial_depth <- add_lonlat(pred_raw_depth,    "Depth restricted (14-250 fm)")
spatial_depth_2km <- add_lonlat(pred_raw_depth_2km, "Depth restricted (2km; 14-250 fm)")
spatial_stations <- add_lonlat(pred_raw_stations, "Fixed survey stations")

# Years to show in plots (restricint number of plots for visual purposes)
key_years <- c(1998, 2002, 2006, 2010, 2014, 2016, 2020, 2021, 2022, 2023, 2024, 2025)
all_years <- sort(unique(dat$Year))

# Consistent colour scale based on 99th percentile of station grid
cpue_max <- quantile(spatial_stations$cpue, 0.99, na.rm = TRUE)
cpue_min <- 0

# Ocean background polygon (This is just the box behind the plot...)
ocean_bbox <- st_bbox(c(xmin = -141, xmax = -132,
                        ymin = 53.5, ymax = 61.5),
                      crs = 4326) %>%
  st_as_sfc() %>%
  st_sf()

# Crop basemap
ak_map_crop <- ak_map %>%
  st_crop(xmin = -141, xmax = -132,
          ymin = 53.5, ymax = 61.5)

# 
# Tile size for fixed stations 
# use smaller tiles since they are points not grid cells
tile_w_stn <- 0.12
tile_h_stn <- 0.08
# 
# Tile size for 15km predictions grids
tile_w_15km <- 0.28   
tile_h_15km <- 0.145  

# Tile size for 2km
tile_w_2km <- 0.033   
tile_h_2km <- 0.018  

# Function to plot the prediction spatially

plot_cpue_maps <- function(spatial_df, grid_name, tile_w, tile_h,
                           ncols = 5, years = key_years) {
  
  grid_cpue_max <- quantile(spatial_df$cpue, 0.99, na.rm = TRUE)
  
  ggplot() +
    geom_sf(data = ocean_bbox,
            fill = "#cfe2f3",
            colour = NA) +
    geom_sf(data  = ak_map_crop,
            fill = "grey80",
            colour = "grey50",
            linewidth = 0.2) +
    geom_sf(data = seo_bounds_m %>% st_transform(4326),
            fill = NA,
            colour = "white",
            linewidth = 0.4) +
    geom_tile(data = spatial_df %>% filter(Year %in% years),
              aes(x = Lon,
                  y = Lat,
                  fill = cpue),
              width  = tile_w,
              height = tile_h) +
    scale_fill_viridis_c(
      option = "magma",
      name   = "CPUE\n(kg/hook)",
      trans  = "sqrt",
      limits = c(0, grid_cpue_max),   # grid-specific scale
      oob    = scales::squish,
      labels = scales::number_format(accuracy = 0.01)
    ) +
    facet_wrap(~Year, ncol = ncols) +
    coord_sf(xlim = c(-140, -133),
             ylim = c(54, 61),
             expand = FALSE) +
    labs(
      # title    = paste("Predicted yelloweye rockfish CPUE –", grid_name),
      # subtitle = paste0("sdmTMB (m2) | sqrt scale | max = ",
      # round(grid_cpue_max, 3), " kg/hook"),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position   = "right",
      legend.key.height = unit(2, "cm"),
      legend.key.width  = unit(0.5, "cm"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      strip.text  = element_text(size = 12, face = "bold"),
      strip.background  = element_rect(fill = "grey95"),
      axis.text.x = element_text(angle = 90),
      # axis.text = element_blank(),
      # axis.ticks  = element_blank(),
      panel.spacing  = unit(0.1, "lines"),
      panel.background  = element_rect(fill = "#cfe2f3"),
      panel.border  = element_rect(colour = "grey70")
    )
}

p_full <- plot_cpue_maps(spatial_full,
                         "Full management area (15 km grid)",
                         tile_w_15km, tile_h_15km)

p_depth <- plot_cpue_maps(spatial_depth,
                          "Depth restricted 14-250 fm (15 km grid)",
                          tile_w_15km, tile_h_15km)

p_depth_2km <- plot_cpue_maps(spatial_depth_2km,
                              "Depth restricted 14-250 fm (2 km grid)",
                              tile_w_2km, tile_h_2km, 
                              years = key_years, ncols = 3)

p_stations  <- plot_cpue_maps(spatial_stations,
                              "Fixed survey stations",
                              tile_w_stn, tile_h_stn)

p_full
p_depth
p_depth_2km
p_stations



ggsave(here("figures", "sdmTMB","sdmTMB_cpue_maps_full.png"),
       p_full, width = 14, height = 10, dpi = 300)
ggsave(here("figures", "sdmTMB","sdmTMB_cpue_maps_depth.png"),
       p_depth, width = 14, height = 10, dpi = 300)
ggsave(here("figures", "sdmTMB","sdmTMB_cpue_maps_stations.png"),
       p_stations, width = 14, height = 10, dpi = 300)
ggsave(here("figures", "sdmTMB","sdmTMB_cpue_maps_stations.png"),
       p_depth_2km, width = 8, height = 10, dpi = 300)


# Single year comparison across all three grids
yr_compare <- 2023

spatial_compare <- bind_rows(
  spatial_full     %>% filter(Year == yr_compare),
  spatial_depth    %>% filter(Year == yr_compare),
  spatial_depth_2km    %>% filter(Year == yr_compare),
  spatial_stations %>% filter(Year == yr_compare)
) %>%
  mutate(
    grid = factor(grid, levels = c(
      "Full management area",
      "Depth restricted (14-250 fm)",
      "Depth restricted (2km; 14-250 fm)",
      "Fixed survey stations"
    )),
    # Assign tile size per grid type
    tile_w = case_when(
      grid == "Fixed survey stations"             ~ tile_w_stn,
      grid == "Depth restricted (2km; 14-250 fm)" ~ tile_w_2km,
      TRUE                                        ~ tile_w_15km
    ),
    tile_h = case_when(
      grid == "Fixed survey stations"             ~ tile_h_stn,
      grid == "Depth restricted (2km; 14-250 fm)" ~ tile_h_2km,
      TRUE                                        ~ tile_h_15km
    )
  )

p_compare_yr <- ggplot() +
  geom_sf(data = ocean_bbox,
          fill = "#cfe2f3",
          colour = NA) +
  geom_tile(data = spatial_compare,
            aes(x = Lon,
                y = Lat,
                fill = cpue,
                width = tile_w,
                height = tile_h)) +
  geom_sf(data = ak_map_crop,
          fill = "grey80",
          colour = "grey50",
          linewidth = 0.2) +
  geom_sf(data = seo_bounds_m %>% st_transform(4326),
          fill = NA,
          colour = "white",
          linewidth = 0.4) +
  scale_fill_viridis_c(
    option = "magma",
    name   = "CPUE\n(kg/hook)",
    trans  = "sqrt",
    limits = c(0, quantile(spatial_compare$cpue, 0.99, na.rm = TRUE)),
    oob    = scales::squish,
    labels = scales::number_format(accuracy = 0.01)
  ) +
  facet_wrap(~grid, ncol = 4) +
  coord_sf(xlim = c(-140, -133),
           ylim = c(54, 61),
           expand = FALSE) +
  labs(
    title    = paste("Predicted yelloweye CPUE by prediction grid –", yr_compare),
    subtitle = "sdmTMB spatio-temporal model (m2)",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position   = "right",
    legend.key.height = unit(2, "cm"),
    strip.text = element_text(size = 9, face = "bold"),
    strip.background  = element_rect(fill = "grey95"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing  = unit(0.2, "lines"),
    panel.background = element_rect(fill = "#cfe2f3")
  )

p_compare_yr

ggsave(here("figures","sdmTMB", paste0("sdmTMB_cpue_grid_comparison_", yr_compare, ".png")),
       p_compare_yr, width = 14, height = 6, dpi = 300)


# Combined index across SEO
# Note that you need to comment out the groupby in predict_index() to run this
result_depth_SEO  <- predict_index(grid_2_depth_d,
                               "Depth restricted (14-250 fm)",
                               se_fit = FALSE, best_model = m2)

SEO_index <- result_depth_SEO$index

SEO_index %>% 
ggplot(aes(Year, mean_cpue)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.5)+
  # facet_wrap(~SEdist, scales = "free_y") +
  # scale_colour_manual(values = c("GAM (Tweedie)"    = "#009E73",
  # "sdmTMB (Tweedie)" = "#0072B2")) +
  scale_color_manual(values = cbpalette[c(3,4)])+
  labs(x = "Year", y = "Standardized CPUE (kg/hook)",
       colour = NULL, linetype = NULL) +
  scale_x_continuous(breaks = seq(1995,2025,5))+
  theme_bw(base_size = 14) +
  theme(legend.position = "top") -> SEO_comb_index_p

# ggsave(here("figures","sdmTMB", "SEOcombined_index_SAFE.png"),
#        SEO_comb_index_p, width = 7, height = 7, dpi = 300)


# Compare boot strap index to sdmTMB index
sdm_index <- read.csv(paste0(here::here(), "/outputs/bridging analysis/IPHC.cpue.SEO_sdmTMB_m1_m2_m3_17July2026.csv"))

boot_index <- read.csv(paste0(here::here(), "/outputs/2026/REMA/IPHC.cpue.SEO_bootindex_2025.csv"))

# Model M2 is the AR1 sdmTMB model with depth and soak smooths
ind_sdm <- sdm_index %>% 
  filter(model == "M.2") %>% 
  select(strata = SEdist, year = Year, cpue = mean_cpue, cv = CV, lwr = lower_95, upr = upper_95) %>% 
  mutate(model = "sdmTMB (Current Method)")

ind_boot <- boot_index %>% 
  select(strata = mngmt.area, year = Year, cpue = WCPUE.bootmean, cv = WCPUE.cv, lwr = WCPUE.lo95ci, upr = WCPUE.hi95ci) %>% 
  mutate(model = "Bootstrapped (2022 Method)")

gam.cpue.df <- read.csv(here("outputs/bridging analysis/IPHC_tweedie_bridge_24June2026.csv"), row.names = NULL) %>% 
  filter(CPUE == "Tweedie index 2026") %>% 
  select(strata = SEdist, year = Year, cpue, cv, lwr = lower, upr = upper) %>% 
  mutate(model = "GAM (2024 Method)")


index_comb <- rbind(ind_sdm, ind_boot, gam.cpue.df)

index_comb %>% 
  ggplot(aes(x = year, y = cpue, color = model, fill = model))+
  # geom_point()+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = lwr, ymax = upr, alpha = model), color = NA)+
  scale_color_colorblind(name = "")+
  scale_fill_colorblind(name = "")+
  scale_alpha_manual(values  = c(0.2,0.2,.6), name = "")+
  labs(x = "Year",
       y = "CPUE (round kg / hook)")+
  facet_wrap(~strata, scales = "free_y")+
  theme_bw(base_size = 14)+
  theme(legend.position = "top") -> boot_sdm_index_p

boot_sdm_index_p

# ggsave(here("figures","sdmTMB", "boot_vs_sdm_index_SAFE.png"),
#        boot_sdm_index_p, width = 7, height = 7, dpi = 300)
