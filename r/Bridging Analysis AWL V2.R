# Aaron Lambert
# 6/24/2026
# 2026 SEO Yelloweye Rockfish - REMA Bridging Analysis
# 
#
# Purpose:
#   Stepwise bridging from the 2022 accepted model to the 2024 proposed base model.
#   Required by SSC (Dec. 2024): show the stepwise effect of each change to
#   the IPHC index separately
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
#  mod.7 "22.2 + std  CPUE + correct stations"
#         CPUE in kg per hook standardized using a spatiotemporal model with the
#         Tweedie distribution in order to accommodate zero inflation, 
#         adjusted for hook saturation and competition, 
#         with the correct EYKT boundary & 
#         the correct lat for EYKT and omitting SEI stations.
#
## OFL/ABC computed at each step using M = 0.02 (SSC-endorsed for 2025).
## M = 0.044 shown as sensitivity on B3 (SSC request, Dec. 2024).

#*****************************************************************************

# Setup ----------------------------------------------------------------

# Install/update rema if needed:
# devtools::install_github("afsc-assessments/rema", dependencies = TRUE,
#                          build_vignettes = TRUE, force = TRUE)

# File path is needed to load rema package
(rema_path <- find.package('rema'))
# (rema_examples <- file.path(rema_path, 'example_scripts'))
# list.files(rema_examples)

library(rema)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(cowplot)
library(here)

YEAR <- 2026


# Palette (colorblind-friendly, consistent with REMA_run.R)
cbpalette <- c("#009E73", "#0072B2", "#E69F00", "#56B4E9",
               "#D55E00", "#CC79A7", "gold", "black", "grey")

# One color per bridging step, named to match model_name strings below
bridge_cols <- c(
  "22.2"                            = cbpalette[1],
  "22.2 + hook adj."                = cbpalette[2],
  "22.2 + lon"                      = cbpalette[3],
  "22.2 + kg/hook"                  = cbpalette[4],
  "22.2 + correct stations"         = cbpalette[5],
  "22.2 + GAM CPUE"                 = cbpalette[6],
  "22.2 + GAM CPUE + correct stations"  = cbpalette[7],
  "22.2 + sdmTMB CPUE"              = cbpalette[8]
)

step_levels <- names(bridge_cols)

# Tier 5 harvest control rule
M_base   <- 0.02    # SSC-endorsed for 2025 specifications
M_alt    <- 0.044   # sensitivity requested by SSC
nonYE_ofl <- 26     # non-YE DSR Tier 6 contributions
nonYE_abc <- 20

# Load data ------------------------------------------------------------

# ROV biomass 

# biomass data from 2022
# bio_2022 <- read.csv(here("Data_processing/Data/SEO_YE_Biomass_subdistrict_2022.csv")) %>%
#   select(strata = Subdistrict, year = Year, biomass = Biomass.mt, cv = Biomass.cv)

# biomass to 2023
bio_2023 <- read.csv(here("data/SEO_YE_Biomass_subdistrict_2024-10-28.csv")) %>%
  select(strata = Subdistrict, year = Year, biomass = Biomass.mt, cv = Biomass.cv)

# IPHC CPUE indices ------------------------------------------------------------

# 2022 Model index
# cpue_2022 <- read.csv(here("Data_processing/Data/IPHC.cpue.SEO_non0_2022.csv")) %>%
#   mutate(cv = CV) %>%                        
#   select(strata = Stratum, year = Year, cpue = Mean, cv)

# Boot strapped indices
boot.cpue.df <- read.csv(here("outputs/bridging analysis/bridge.cpue.boot 24June2026.csv"), row.names = NULL)

# GAM standardized CPUE
gam.cpue.df <- read.csv(here("outputs/bridging analysis/IPHC_tweedie_bridge_24June2026.csv"), row.names = NULL)

sdm.cpue.df <- read.csv(here("outputs/bridging analysis/IPHC.cpue.SEO_sdmTMB_m1_m2_m3_17July2026.csv"), 
                        row.names = NULL) %>% 
  filter(model == "M.2")


# sdmTMB standardized CPUE

# bootstrap cpue
cpue_bootstrap_2024 <-  boot.cpue.df %>%
  filter(name == "22.2") %>% 
  select(strata = mngmt.area, year = Year, cpue = CPUE.bootmean, cv = CPUE.cv)

# bootstrap cpue +adj for hook saturation and competition
cpue_hook_adj <- boot.cpue.df %>%
  filter(name == "22.2_hhokadj_kghook") %>% 
  select(strata = mngmt.area, year = Year, cpue = CPUE.bootmean, cv = CPUE.cv)

# bootstrap cpue + lon for hook saturation and competition with incorrect longitude
cpue_lon <- boot.cpue.df %>%
  filter(name == "22.2_lon") %>% 
  select(strata = mngmt.area, year = Year, cpue = CPUE.bootmean, cv = CPUE.cv)

# bootstrap kg/hook cpue + kg/hook for hook saturation and competition, correct boundary
cpue_kg_hook <- boot.cpue.df %>%
  filter(name == "22.2_hhokadj_kghook") %>% 
  select(strata = mngmt.area, year = Year, cpue = WCPUE.bootmean, cv = WCPUE.cv)

cpue_2026 <- boot.cpue.df %>%
  filter(name == "22.2 + correct stations") %>% 
  select(strata = mngmt.area, year = Year, cpue = WCPUE.bootmean, cv = WCPUE.cv)

# GAM tweedie std CPUE for hook saturation and competition and kg/hook
cpue_GAM <- gam.cpue.df %>% 
  filter(CPUE=="Tweedie index") %>%
  select(strata = SEdist, year=Year, cpue = cpue, cv = cv) 

cpue_GAM_2026 <- gam.cpue.df %>% 
  filter(CPUE=="Tweedie index 2026") %>%
  select(strata = SEdist, year=Year, cpue = cpue, cv = cv) 

# sdmTMB std CPUE for hook saturation and competition and kg/hook
cpue_SDM <- sdm.cpue.df%>% 
  select(strata = SEdist, year=Year, cpue = mean_cpue, cv = CV) 

# Fit models --------------------------------------------------

fit_rema_bridge <- function(model_name, biomass_dat, cpue_dat) {
  inp <- prepare_rema_input(
    model_name = model_name,
    multi_survey = TRUE,
    biomass_dat = biomass_dat,
    cpue_dat = cpue_dat,
    extra_biomass_cv = list(assumption = "extra_cv",
                            pointer_extra_biomass_cv = c(1, 1, 1, 1)),
    extra_cpue_cv = NULL,
    PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
    q_options = list(pointer_q_cpue    = c(1, 2, 3, 4))
  )
  mod  <- fit_rema(inp)
  conv <- check_convergence(mod)
  list(input = inp, 
       model = mod, 
       output = tidy_rema(mod))
}

# mod.0
mod.0 <- fit_rema_bridge("22.2",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_bootstrap_2024)

# mod.1
mod.1 <- fit_rema_bridge("22.2 + hook adj.",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_hook_adj)

# mod.2
mod.2 <- fit_rema_bridge("22.2 + lon",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_lon)

# mod.3
mod.3 <- fit_rema_bridge("22.2 + kg/hook",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_kg_hook)

# mod.4
mod.4 <- fit_rema_bridge("22.2 + correct stations",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_2026)

# mod.5
mod.5 <- fit_rema_bridge("22.2 + GAM CPUE",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_GAM)

# mod.6
mod.6 <- fit_rema_bridge("22.2 + GAM CPUE + correct stations",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_GAM_2026)

# mod.7
mod.7 <- fit_rema_bridge("22.2 + sdmTMB CPUE",
                         biomass_dat = bio_2023,
                         cpue_dat    = cpue_SDM)

all_steps <- list(mod.0, mod.1, mod.2, mod.3, mod.4, mod.5, mod.6, mod.7)

# Harvest specifications at each step ----------------------------------

calc_specs <- function(step, M = M_base) {
  yr <- max(step$output$total_predicted_biomass$year)
  step$output$total_predicted_biomass %>%
    filter(year == yr) %>%
    transmute(
      Step = step$input$model_name,
      Ref_year = yr,
      M = M,
      Biomass_mt = round(pred),
      Bio_lci = round(pred_lci),
      Bio_uci = round(pred_uci),
      YE_OFL = round(M * pred),
      YE_ABCmax = round(0.75 * M * pred),
      DSR_OFL = round(M * pred) + nonYE_ofl,
      DSR_ABCmax = round(0.75 * M * pred) + nonYE_abc
    )
}

# Use the above function to get reference points and biomass est
specs_M002 <- bind_rows(lapply(all_steps, calc_specs, M = M_base))
specs_M044 <- bind_rows(lapply(all_steps, calc_specs, M = M_alt))  

# Look at the tables
specs_M002 %>% select(Step, Biomass_mt, YE_OFL, YE_ABCmax, DSR_OFL, DSR_ABCmax)
print(specs_M044 %>% select(Step, M, Biomass_mt, YE_OFL, YE_ABCmax, DSR_OFL, DSR_ABCmax))

# Save the tables
# write.csv(specs_M002,
#           here(paste0("outputs/bridging analysis/harvest_specs_bridging_M002_", YEAR, ".csv")),
#           row.names = FALSE)

# write.csv(specs_M044,
#           here(paste0("REMA/Output/bridging/harvest_specs_M044_sensitivity_", YEAR, ".csv")),
#           row.names = FALSE)

# Outputs for plotting -----------------------------------------

# Functions to ger plotting values
get_total  <- function(s) s$output$total_predicted_biomass %>%
  mutate(model_name = s$input$model_name)
get_strata <- function(s) s$output$biomass_by_strata      %>%
  mutate(model_name = s$input$model_name)
get_cpue   <- function(s) s$output$cpue_by_strata         %>%
  mutate(model_name = s$input$model_name)

# Use lapply to get the values from the step models
all_total  <- bind_rows(lapply(all_steps, get_total))  %>%
  mutate(model_name = factor(model_name, levels = step_levels))
all_strata <- bind_rows(lapply(all_steps, get_strata)) %>%
  mutate(model_name = factor(model_name, levels = step_levels))
all_cpue   <- bind_rows(lapply(all_steps, get_cpue))   %>%
  mutate(model_name = factor(model_name, levels = step_levels))

# Plots ----------------------------------------------------------------

linetypes <- c(1,2,4,6,8, 5,10,3)

# Total biomass trajectory 
p_total <- ggplot(all_total,
                  aes(x = year, color = model_name, fill = model_name,
                      linetype = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), alpha = 0.15, color = NA) +
  geom_line(aes(y = pred), linewidth = 0.8) +
  scale_color_manual(values = bridge_cols, name = NULL) +
  scale_fill_manual( values = bridge_cols, name = NULL) +
  scale_linetype_manual(values = linetypes, name = NULL) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(1994, YEAR, by = 1)) +
  xlab("Year") + ylab("Total estimated biomass (t)") +
  # ggtitle("SEO yelloweye rockfish — bridging analysis",
  #         subtitle = "Total YE biomass") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_total

ggsave(p_total,
       filename = here(paste0("figures/bridging analysis/bridge_total_biomass_", YEAR, ".png")),
       width = 8, height = 5, units = "in", bg = "white", dpi = 300)

# Total biomass trajectory no error
p_total_ne <- ggplot(all_total,
                  aes(x = year, color = model_name, fill = model_name)) +
  # geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), alpha = 0.15, color = NA) +
  geom_line(aes(y = pred), linewidth = 1.5) +
  scale_color_manual(values = bridge_cols, name = NULL) +
  scale_fill_manual( values = bridge_cols, name = NULL) +
  # scale_linetype_manual(values = linetypes, name = NULL) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(1994, YEAR, by = 1)) +
  xlab("Year") + ylab("Total estimated biomass (t)") +
  # ggtitle("SEO yelloweye rockfish — bridging analysis",
  #         subtitle = "Total YE biomass") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_total_ne

ggsave(p_total_ne,
       filename = here(paste0("figures/bridging analysis/bridge_total_biomass_ne_", YEAR, ".png")),
       width = 8, height = 5, units = "in", bg = "white", dpi = 300)

# Biomass fits by strata (observations plotted for std. CPUE)
obs_strata <- all_strata %>% filter(model_name == "22.2")

p_strata <- ggplot(all_strata,
                   aes(x = year, color = model_name, fill = model_name,
                       linetype = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), alpha = 0.15, color = NA) +
  geom_line(aes(y = pred), linewidth = 0.6) +
  geom_point(data = obs_strata, aes(y = obs), color = "black", size = 1.2) +
  geom_errorbar(data = obs_strata, aes(ymin = obs_lci, ymax = obs_uci),
                color = "black", linewidth = 0.35, width = 0) +
  facet_wrap(~strata, scales = "free_y", ncol = 2) +
  scale_color_manual(values = bridge_cols, name = NULL) +
  scale_fill_manual( values = bridge_cols, name = NULL) +
  scale_linetype_manual(values = linetypes, name = NULL) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(1994, YEAR, by = 6)) +
  xlab("Year") + ylab("Biomass (t)") +
  # ggtitle("ADF&G ROV survey fits by management unit") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_strata

ggsave(p_strata,
       filename = here(paste0("figures/bridging analysis/bridge_biomass_strata_", YEAR, ".png")),
       width = 10, height = 5, units = "in", bg = "white", dpi = 300)

p_strata_ne <- ggplot(all_strata,
                   aes(x = year, color = model_name, fill = model_name)) +
  # geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), alpha = 0.15, color = NA) +
  geom_line(aes(y = pred), linewidth = 1) +
  geom_point(data = obs_strata, aes(y = obs), color = "black", size = 1.2) +
  geom_errorbar(data = obs_strata, aes(ymin = obs_lci, ymax = obs_uci),
                color = "black", linewidth = 0.35, width = 0) +
  facet_wrap(~strata, scales = "free_y", ncol = 1) +
  scale_color_manual(values = bridge_cols, name = NULL) +
  scale_fill_manual( values = bridge_cols, name = NULL) +
  # scale_linetype_manual(values = linetypes, name = NULL) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(1994, YEAR, by = 6)) +
  xlab("Year") + ylab("Biomass (t)") +
  # ggtitle("ADF&G ROV survey fits by management unit") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_strata_ne

ggsave(p_strata_ne,
       filename = here(paste0("figures/bridging analysis/bridge_biomass_strata_ne_", YEAR, ".png")),
       width = 6, height = 8, units = "in", bg = "white", dpi = 300)

# IPHC CPUE fits by strata

p_cpue <- ggplot(all_cpue,
                 aes(x = year, color = model_name, fill = model_name,
                     linetype = model_name)) +
  # geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), alpha = 0.15, color = NA) +
  geom_line(aes(y = pred), linewidth = 0.6) +
  geom_point(aes(y = obs), size = 1) +
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), linewidth = 0.3, width = 0) +
  facet_wrap(~strata, scales = "free_y", ncol = 2) +
  scale_color_manual(values = bridge_cols, name = NULL) +
  scale_fill_manual( values = bridge_cols, name = NULL) +
  scale_linetype_manual(values = linetypes, name = NULL) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(1998, YEAR, by = 6)) +
  xlab("Year") + ylab("IPHC setline survey CPUE (kg/hook)") +
  ggtitle("IPHC CPUE index fits by management unit") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_cpue

ggsave(p_cpue,
       filename = here(paste0("figures/bridging analysis/bridge_cpue_strata_", YEAR, ".png")),
       width = 10, height = 5, units = "in", bg = "white", dpi = 300)

# Harvest specifications dot-plot across steps
specs_long <- specs_M002 %>%
  select(Step, YE_OFL, YE_ABCmax) %>%
  tidyr::pivot_longer(c(YE_OFL, YE_ABCmax), names_to = "Metric", values_to = "value") %>%
  mutate(Metric = recode(Metric, "YE_OFL" = "OFL", "YE_ABCmax" = "maxABC"),
         Step   = factor(Step, levels = step_levels))

p_specs <- ggplot(specs_long,
                  aes(x = Step, y = value, color = Metric, shape = Metric, group = Metric)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("OFL" = cbpalette[3], "maxABC" = cbpalette[8])) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 14)) +
  xlab(NULL) + ylab("Yelloweye rockfish (t)") +
  # ggtitle("Harvest specifications across bridging steps",
  #         subtitle = paste0("Tier 5, M = ", M_base)) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_specs

ggsave(p_specs,
       filename = here(paste0("figures/bridging analysis/bridge_harvest_specs_", YEAR, ".png")),
       width = 7, height = 5, units = "in", bg = "white", dpi = 300)

# Parameter estimates --------------------------------------------------

all_params <- bind_rows(lapply(all_steps, function(s) {
  s$output$parameter_estimates %>%
    mutate(model_name = s$input$model_name,
           across(where(is.numeric), ~round(., 5)))
}))

print(all_params)
# write.csv(all_params,
#           here(paste0("REMA/Output/bridging/parameter_estimates_bridging_", YEAR, ".csv")),
#           row.names = FALSE)

# SAFE summary table ---------------------------------------------------

safe_table <- specs_M002 %>%
  mutate(`Biomass 95% CI` = paste0("(", Bio_lci, "-", Bio_uci, ")")) %>%
  select(Step, Biomass_mt, `Biomass 95% CI`,
         YE_OFL, YE_ABCmax, DSR_OFL, DSR_ABCmax)

print(safe_table)

write.csv(safe_table,
          here(paste0("outputs/bridging analysis/SAFE_bridging_summary_", YEAR, ".csv")),
          row.names = FALSE)


cpue_SDM %>%
  left_join(cpue_GAM, by = c("year", "strata")) %>%
  filter(year >= 2018) %>%
  select(year, strata, sdmTMB = cpue.x, GAM = cpue.y) %>%
  mutate(ratio = sdmTMB / GAM) %>%
  arrange(strata, year)



# Plot the indices
# mod.0
cpue_bootstrap_2024$model <- "22.2"

# mod.1
cpue_hook_adj$model <- "22.2 + hook adj."

# mod.2
cpue_lon$model <- "22.2 + lon"

# mod.3
cpue_kg_hook$model <- "22.2 + kg/hook"

# mod.4
cpue_2026$model <- "22.2 + correct stations"

# mod.5
cpue_GAM$model <- "22.2 + GAM CPUE"

# mod.6
cpue_GAM_2026$model <- "22.2 + GAM CPUE + correct stations"

# mod.7
cpue_SDM$model <- "22.2 + sdmTMB CPUE"

Index_all <- rbind(cpue_bootstrap_2024,
                   cpue_hook_adj,
                   cpue_lon,
                   cpue_kg_hook,
                   cpue_2026,
                   cpue_GAM,
                   cpue_GAM_2026,
                   cpue_SDM)

Index_all$model <- factor(Index_all$model, levels = c("22.2",
                                     "22.2 + hook adj.",
                                     "22.2 + lon",
                                     "22.2 + kg/hook",
                                     "22.2 + correct stations",
                                     "22.2 + GAM CPUE",
                                     "22.2 + GAM CPUE + correct stations",
                                     "22.2 + sdmTMB CPUE"))

# Compare the different indices for the SAFE appendix
Index_all %>% 
  ggplot(aes(x = year, y = cpue, color = model))+
  geom_line(size = 1.2, alpha = 0.75)+
  # geom_line(data = df, aes(x = Year, y = cpue), linewidth = 2, inherit.aes = F)+
  scale_x_continuous(breaks = seq(1998,2026,4))+
  scale_color_manual(values = cbpalette, name = "")+
  facet_wrap(~strata, scales = "free_y", ncol = 1)+
  labs(x = "Year", y = "CPUE")+
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))-> index_compare_p

index_compare_p

ggsave(index_compare_p,
       filename = here(paste0("figures/bridging analysis/bridge_index_compare", YEAR, ".png")),
       width = 6, height = 8, units = "in", bg = "white", dpi = 300)

# Compare the different indice cvs for the SAFE appendix
Index_all %>% 
  ggplot(aes(x = year, y = cv, color = model))+
  geom_line(size = 1.2, alpha = 0.75)+
  scale_color_manual(values = cbpalette, name = "")+
  scale_x_continuous(breaks = seq(1998,2026,4))+
  facet_wrap(~strata, scales = "free_y", ncol = 1)+
  labs(x = "Year", y = "CV")+
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))-> cv_compare_p

cv_compare_p

ggsave(cv_compare_p,
       filename = here(paste0("figures/bridging analysis/bridge_cv_compare", YEAR, ".png")),
       width = 6, height = 8, units = "in", bg = "white", dpi = 300)



# CV sensitivity analysis ------------------------------------------------------

# # Mean CV across all non-sdmTMB indices by stratum
ref_cv <- Index_all %>%
  filter(model != "22.2 + sdmTMB CPUE") %>%
  group_by(strata) %>%
  summarise(overall_mean_cv = mean(cv, na.rm = TRUE), .groups = "drop")
 
# Mean by year and stratum
ref_cv_yr <- Index_all %>%
  filter(model != "22.2 + sdmTMB CPUE") %>%
  group_by(strata, year) %>%
  summarise(mean_cv = mean(cv, na.rm = TRUE), .groups = "drop")

# Join back to sdmTMB index by year and stratum
cpue_SDM_refcv <- cpue_SDM %>%
  left_join(ref_cv_yr, by = c("strata", "year")) %>%
  left_join(ref_cv, by = c("strata")) %>% 
  mutate(cv = coalesce(mean_cv, overall_mean_cv)) %>%  # fall back to sdmTMB cv if no ref available
  select(strata, year, cpue, cv)

# CV multipliers to test (applied to sdmTMB CVs)
multipliers <- c(1, 2, 3, 4, 5)

# Fit multiplier sensitivity models
fit_sdm_cv_sens <- function(multiplier, biomass_dat, cpue_sdm_base) {
  cpue_scaled <- cpue_sdm_base %>%
    mutate(cv = cv * multiplier)
  fit_rema_bridge(
    model_name = paste0("sdmTMB CV x", multiplier),
    biomass_dat = biomass_dat,
    cpue_dat = cpue_scaled
  )
}

sens_models <- lapply(multipliers, fit_sdm_cv_sens,
                      biomass_dat   = bio_2023,
                      cpue_sdm_base = cpue_SDM)

names(sens_models) <- paste0("CV_x", multipliers)

# Also fit the reference CV model
sens_ref <- fit_rema_bridge("sdmTMB CV = mean other indices",
                            biomass_dat = bio_2023,
                            cpue_dat    = cpue_SDM_refcv)

all_sens <- c(sens_models, list(ref_cv_model = sens_ref))

# Harvest specs for each sensitivity
sens_step_levels <- c(paste0("sdmTMB CV x", multipliers),
                      "sdmTMB CV = mean other indices")

specs_sens <- bind_rows(lapply(all_sens, calc_specs, M = M_base)) %>%
  mutate(Step = factor(Step, levels = sens_step_levels))

specs_sens %>% select(Step, Biomass_mt, Bio_lci, Bio_uci, YE_OFL, YE_ABCmax)

# Total biomass trajectories across CV sensitivities
sens_cols <- c(
  "sdmTMB CV x1" = cbpalette[8],
  "sdmTMB CV x2" = cbpalette[2],
  "sdmTMB CV x3" = cbpalette[3],
  "sdmTMB CV x4" = cbpalette[4],
  "sdmTMB CV x5" = cbpalette[5],
  "sdmTMB CV = mean other indices" = cbpalette[1]
)

all_sens_total <- bind_rows(lapply(all_sens, get_total)) %>%
  mutate(model_name = factor(model_name, levels = sens_step_levels))

p_sens_total <- ggplot(all_sens_total,
                       aes(x = year, color = model_name, fill = model_name)) +
  # geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), alpha = 0.15, color = NA) +
  geom_line(aes(y = pred), linewidth = 2) +
  # geom_line(data = get_total(mod.6),
  #           aes(y = pred), color = "black", linewidth = 1,
  #           linetype = "solid", inherit.aes = TRUE) +
  scale_color_manual(values = sens_cols, name = NULL) +
  scale_fill_manual(values  = sens_cols, name = NULL) +
  # scale_linetype_manual(values = c(2, 4, 5, 6, 8, 1), name = NULL) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(1994, YEAR, by = 4)) +
  xlab("Year") + ylab("Total estimated biomass (t)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_sens_total

ggsave(p_sens_total,
       filename = here(paste0("figures/bridging analysis/sdmTMB_CV_sensitivity_biomass_", YEAR, ".png")),
       width = 8, height = 5, units = "in", bg = "white", dpi = 300)

# Harvest specs sensitivity plot
specs_sens_long <- specs_sens %>%
  select(Step, YE_OFL, YE_ABCmax) %>%
  tidyr::pivot_longer(c(YE_OFL, YE_ABCmax), names_to = "Metric", values_to = "value") %>%
  mutate(Metric = recode(Metric, "YE_OFL" = "OFL", "YE_ABCmax" = "maxABC"))



ref_specs <- calc_specs(mod.7, M = M_base)  

p_sens_specs <- ggplot(specs_sens_long,
                       aes(x = Step, y = value, color = Metric,
                           shape = Metric, group = Metric)) +
  geom_hline(yintercept = ref_specs$YE_OFL,
             linetype = "dashed", color = cbpalette[3], linewidth = 0.7) +
  geom_hline(yintercept = ref_specs$YE_ABCmax,
             linetype = "dashed", color = cbpalette[8], linewidth = 0.7) +
  annotate("text", x = 0.6, y = ref_specs$YE_OFL + 10,
           label = "sdmTMB base OFL", size = 3, color = cbpalette[3], hjust = 0) +
  annotate("text", x = 0.6, y = ref_specs$YE_ABCmax + 10,
           label = "sdmTMB base maxABC", size = 3, color = cbpalette[8], hjust = 0) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("OFL" = cbpalette[3], "maxABC" = cbpalette[8])) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
  xlab(NULL) + ylab("Yelloweye rockfish (t)") +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_sens_specs

ggsave(p_sens_specs,
       filename = here(paste0("figures/bridging analysis/sdmTMB_CV_sensitivity_specs_", YEAR, ".png")),
       width = 7, height = 5, units = "in", bg = "white", dpi = 300)

# Summary table
write.csv(specs_sens %>%
            select(Step, Biomass_mt, Bio_lci, Bio_uci, YE_OFL, YE_ABCmax, DSR_OFL, DSR_ABCmax),
          here(paste0("outputs/bridging analysis/sdmTMB_CV_sensitivity_", YEAR, ".csv")),
          row.names = FALSE)




