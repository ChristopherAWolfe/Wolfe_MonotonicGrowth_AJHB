################################################################################
#                                                                              #
# The following script creates Figure 6 in Wolfe (2026), Monotonic Growth      #
# Models. The bioarchaeological data derives from Wolfe and Herrmann 2022.     #
# The US SVAD monotonic effects model is imported for comparison.              #
#                                                                              #
################################################################################

# Package Load
library(magrittr)
library(tidyverse)
library(brms)
library(posterior)
library(patchwork)

# Data Import + Cleaning
hunter <- read.csv("hunter_dentalDev_FDL.csv")
freed_teeth <- read.csv("41DL316_RECODE_tooth.scores_ANALYSIS.csv")
freed_all <- read.csv("freedmans_REAL_CSV.csv")
freed <- inner_join(freed_teeth, freed_all, by = c("ID" = "PIN"))
freed %<>% select(ID, M1, FLEN)
hunter %<>% select(ID, M1, FDL)
freed %<>% mutate(M1_recode = ifelse(M1 >= 8, M1 - 1, M1))
hunter %<>% mutate(M1_recode = ifelse(M1 >= 8, M1 - 1, M1))
colnames(freed)[4] <- "man_M1_L"
colnames(hunter)[4] <- "man_M1_L"

colnames(freed)[3] <- "FDL_L"
colnames(hunter)[3] <- "FDL_L"
hunter <- hunter[-c(14, 28),]

# Model Import and Creation

## SVAD US Model (Fit previously)
us_fit_m1_mo <- readRDS("models/us_mo.rds")

## Freedman's Model
freed_mo <- brm(bf(FDL_L ~ mo(man_M1_L), sigma ~ mo(man_M1_L)), data = freed, 
                backend = "cmdstanr")

## Georgia Model
huunter_mo <- brm(bf(FDL_L ~ mo(man_M1_L), sigma ~ mo(man_M1_L)), data = hunter, 
                  backend = "cmdstanr", control = list(adapt_delta = 0.99999, max_treedepth = 12))

# Plot Prep

## Helper to extract predictions for one model + sample
get_pred_df <- function(model, data, sample_label) {
  newdat <- data.frame(man_M1_L = sort(unique(data$man_M1_L)))
  
  pred <- predict(model, newdata = newdat, probs = c(0.025, 0.975)) |>
    as.data.frame() |>
    rename(pi_lo = Q2.5, pi_hi = Q97.5)
  
  fit <- fitted(model, newdata = newdat, probs = c(0.025, 0.975)) |>
    as.data.frame() |>
    rename(ci_lo = Q2.5, ci_hi = Q97.5, mu = Estimate)
  
  bind_cols(newdat, fit[, c("mu", "ci_lo", "ci_hi")],
            pred[, c("pi_lo", "pi_hi")]) |>
    mutate(sample = sample_label)
}

## Build combined prediction frame
plot_df <- bind_rows(
  get_pred_df(us_fit_m1_mo, us_sub,  "SVAD"),
  get_pred_df(freed_mo,    freed,  "Freedman's"),
  get_pred_df(huunter_mo,   hunter, "Hunter")
)

## Combined raw data
raw_df <- bind_rows(
  us_sub  |> select(man_M1_L, FDL_L) |> mutate(sample = "SVAD"),
  freed |> select(man_M1_L, FDL_L) |> mutate(sample = "Freedman's"),
  hunter |> select(man_M1_L, FDL_L) |> mutate(sample = "Hunter")
)

## Consistent factor ordering for legend
plot_df$sample <- factor(plot_df$sample, levels = c("SVAD", "Freedman's", "Hunter"))
raw_df$sample  <- factor(raw_df$sample,  levels = c("SVAD", "Freedman's", "Hunter"))

## Palette
pal <- c("SVAD" = "#1b6ca8", "Freedman's" = "goldenrod", "Hunter" = "coral")

# Plot Growth Profiles
bioarch <- ggplot(plot_df, aes(x = man_M1_L, y = mu,
                               color = sample, fill = sample, group = sample)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.28, color = NA) +
  # Growth trajectory
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  # Raw data
  geom_jitter(data = raw_df, aes(x = man_M1_L, y = FDL_L, color = sample),
              width = 0.10, height = 0, alpha = 0.35, size = 1.2,
              inherit.aes = FALSE) +
  scale_color_manual(values = pal, name = "Sample") +
  scale_fill_manual(values = pal,  name = "Sample") +
  scale_x_continuous(breaks = sort(unique(raw_df$man_M1_L))) +
  labs(x = "Maxillary M1 stage",
       y = "FDL (mm)") +
  theme_classic(base_size = 12) + theme(legend.position = "bottom")

# Prep Z-Scores
us  <- us_sub |> mutate(sample = "SVAD")
fd  <- freed |> mutate(sample = "Freedman's")
ht <- hunter |> mutate(sample = "Hunter")

## ---- STAGE: model-based mu and sigma at each stage ----------------------
stage_grid <- data.frame(man_M1_L = sort(unique(us$man_M1_L)))

mu_stage    <- fitted(us_fit_m1_mo, newdata = stage_grid)[, "Estimate"]
sigma_stage <- fitted(us_fit_m1_mo, newdata = stage_grid, dpar = "sigma")[, "Estimate"]

ref_stage <- data.frame(
  man_M1_L = stage_grid$man_M1_L,
  mu_stage = mu_stage,
  sigma_stage = sigma_stage
)

## ---- Compute z-scores for SA + COL --------------------------------------
score <- function(df) {
  df |>
    left_join(ref_stage, by = "man_M1_L") |>
    mutate(z_stage = (FDL_L - mu_stage) / sigma_stage)
}

freed_z  <- score(fd)
hunter_z <- score(ht)
zs    <- bind_rows(freed_z, hunter_z)

FD_COLOR  <- "goldenrod"
HT_COLOR <- "coral"
GREY_REF  <- "#9aa0a6"
pal       <- c("Freedman's" = FD_COLOR, "Hunter" = HT_COLOR)

zs %<>% select(man_M1_L, z_stage, sample) %>% na.omit() %>% droplevels()

# Plot Z-Scores
p_stage <- ggplot(zs, aes(x = factor(man_M1_L), y = z_stage, color = sample)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = 2,
           fill = GREY_REF, alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = 1,
           fill = GREY_REF, alpha = 0.10) +
  geom_hline(yintercept = 0,  color = GREY_REF, linewidth = 0.4) +
  geom_hline(yintercept = c(-2, 2), color = GREY_REF,
             linewidth = 0.3, linetype = "dashed") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.45),
             size = 2.4, alpha = 0.85,
             stroke = 0.6, shape = 21,
             aes(fill = sample), color = "white") +
  scale_fill_manual(values = pal, name = "Sample") +
  scale_color_manual(values = pal, name = "Sample") +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(x = "Mandibular M1 stage",
       y = "Z-Score") +
  theme_classic(base_size = 11) + theme(legend.position = "none") + guides(fill = "none", color = "none")
p_stage

fig <- (bioarch | p_stage) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Fig6_V2.pdf", fig, width = 9, height = 7, device = "pdf", units = "in")

###############################END##############################################
