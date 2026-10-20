################################################################################
#                                                                              #
# The following script creates Figure 7 in Wolfe (2026), Monotonic Growth      #
# Models. All paleoanthropological data is derived from their intial           #
# publications (See main text). Comparisons were made against the previously   #
# fit US SVAD model.                                                           #
#                                                                              #
################################################################################

# Package Load
library(magrittr)
library(tidyverse)
library(brms)
library(posterior)
library(patchwork)

# Load US Model
us_fit_m1_mo <- readRDS("models/us_mo.rds")

# Import Stature Data
dat2 <- readRDS("/Users/christopherwolfe/Downloads/US_0-21_FINAL.rds")
dat3 <- inner_join(us_sub, dat2, by = c("SVAD_identifier" = "SVAD_medrec"))

# Femur / Stature Model
us_fit_fs <- brm(bf((FDL_L/stature) ~ mo(man_M1_L), sigma ~ mo(man_M1_L)), data = dat3, 
                  backend = "cmdstanr")
# Prep Plots
## Grid of ordinal levels
newdat <- data.frame(man_M1_L = sort(unique(dat3$man_M1_L)))

# Prediction intervals (include sigma — show the heteroskedasticity)
pred <- predict(us_fit_fs, newdata = newdat, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  rename(pi_lo = Q2.5, pi_hi = Q97.5)

# Fitted intervals (uncertainty in the mean only)
fit <- fitted(us_fit_fs, newdata = newdat, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  rename(ci_lo = Q2.5, ci_hi = Q97.5, mu = Estimate)

plot_df <- bind_cols(newdat, fit[, c("mu", "ci_lo", "ci_hi")],
                     pred[, c("pi_lo", "pi_hi")])

# Plot Results
fossils <- tibble::tribble(
  ~specimen,      ~x,  ~fdl,  ~stat_lo, ~stat_hi,
  "Dederiyeh 1",     5,  159.7,  792,      817,
  "Lagar Velho 1",   8,  199,    910,      978,
  "Qafzeh 10",     9,  233,    981,      1063,
  "El Sidrón J1",    11,  250,    1049.8,   1200,
  "KNM-WT 15000",  13,  400,    1410,     1693
) |>
  dplyr::mutate(ymin = fdl / stat_hi, ymax = fdl / stat_lo,
                specimen = factor(specimen, levels = specimen))

fossil_cols <- c("Dederiyeh 1"   = "darkgreen",
                 "Lagar Velho 1" = "firebrick",
                 "Qafzeh 10"   = "goldenrod",
                 "El Sidrón J1"   = "antiquewhite3",
                 "KNM-WT 15000" = "darkorange")

# Shared fossil layer + scale, added to both plots
fossil_layer <- list(
  geom_errorbar(data = fossils,
                aes(x = x, ymin = ymin, ymax = ymax, color = specimen),
                linewidth = 1.5, inherit.aes = FALSE),
  scale_color_manual(values = fossil_cols, name = NULL)
)

fig3a <- ggplot(plot_df, aes(x = man_M1_L, y = mu)) +
  geom_linerange(aes(ymin = pi_lo, ymax = pi_hi),
                 color = "steelblue", linewidth = 0.6) +
  geom_linerange(aes(ymin = ci_lo, ymax = ci_hi),
                 color = "steelblue4", linewidth = 2) +
  geom_point(color = "steelblue4", size = 2.5) +
  fossil_layer +
  geom_jitter(data = dat3, aes(x = man_M1_L, y = FDL_L / stature),
              width = 0.08, height = 0, alpha = 0.5,
              size = 0.5, color = "grey25", inherit.aes = FALSE) +
  scale_x_continuous(breaks = sort(unique(dat3$man_M1_L))) +
  labs(x = "Mandibular M1 stage", y = "FDL / stature") +
  theme_classic(base_size = 12)

fig3b <- ggplot(plot_df, aes(x = man_M1_L, y = mu)) +
  geom_ribbon(aes(ymin = pi_lo, ymax = pi_hi),
              fill = "steelblue", alpha = 0.15) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = "steelblue", alpha = 0.35) +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  geom_point(color = "steelblue4", size = 2.2) +
  fossil_layer +
  geom_jitter(data = dat3, aes(x = man_M1_L, y = FDL_L / stature),
              width = 0.12, height = 0, alpha = 0.5,
              size = 1.4, color = "grey25", inherit.aes = FALSE) +
  scale_x_continuous(breaks = sort(unique(dat3$man_M1_L))) +
  labs(x = "Mandibular M1 stage", y = "FDL / stature") +
  theme_classic(base_size = 12)

fig <- (fig3a | fig3b) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Fig7_V2.pdf", fig, width = 9, height = 7, device = "pdf", units = "in")

###############################END##############################################
