#####################################
## @Description: Figure 5 - GLMM analysis for pertussis incidence vs vaccination program
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-22 11:08:31
## @LastEditors: Li Kangguo
## @LastEditTime: 2025-11-03 10:52:17
#####################################

library(tidyverse)
library(glmmTMB)
library(gtsummary)
library(broom.mixed)
library(openxlsx)
library(flextable)
library(gt)
library(brms)

rm(list = ls())

# data --------------------------------------------------------------------

region_names <- c("Global", 'WHO region',
                  "African Region", "Eastern Mediterranean Region", "European Region",
                  "Region of the Americas", "South-East Asia Region", "Western Pacific Region")

## Vaccine data ------------------------------------------------------------

load('./Output/DataAll.RData')

## bind --------------------------------------------------------------------

DataAll <- DataAll |>
     select(location_id, location_name, Location, Location_ID, ISO2, 
            VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant,
            contains('WHO_Inci_'), contains('GBD_Inci_'), contains('Population_'))

# attach Region mapping: https://cdn.who.int/media/docs/default-source/air-pollution-documents/air-quality-and-health/un-agencies-region-classification-for-country.xlsx?sfvrsn=289af35f_3
unpop <- read.xlsx("./Data/CLASS_2025_07_02.xlsx")
unpop_map <- unpop |>
     select(Code, Income.group) |>
     rename(ISO3 = Code, Region = Income.group)

DataAll <- DataAll |>
     left_join(unpop_map, by = c(Location_ID = 'ISO3'))

rm(unpop, unpop_map)

# attach UHC service coverage index data: https://data.who.int/indicators/i/3805B1E/9A706FD
uhc_data <- read.csv("./Data/UHC service coverage index.csv") |>
     select(Year = DIM_TIME, GEO_NAME_SHORT, uhc_index = INDEX_N)

# using 2021 data to replace 2024 value
uhc_2024 <- uhc_data |>
     filter(Year == 2021) |>
     mutate(Year = 2024)

uhc_data <- bind_rows(uhc_data, uhc_2024) |> 
     # cut index to very high (80+), high (60-79), medium (40-59), low (20-39), very low (<20)
     # https://www.who.int/news-room/questions-and-answers/item/tracking-universal-health-coverage
     mutate(uhc_cat = cut(uhc_index,
                          breaks = c(-Inf, 20, 40, 60, 80, Inf),
                          labels = c('Very low', 'Low', 'Medium', 'High', 'Very high')))

rm(uhc_2024)

## data clean -------------------------------------------------------------------

## Build long dataset stacking WHO and GBD Inci data with a `source` column
# pivot WHO and GBD incidence columns separately, then bind_rows into a long table
who_long <- DataAll |>
     select(location_id, Location_ID, location_name, starts_with('WHO_Inci_'), VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant) |>
     pivot_longer(cols = matches('WHO_Inci_\\d{4}'), names_to = 'InciVar', values_to = 'Incidence') |>
     mutate(Year = as.integer(str_extract(InciVar, '\\d{4}')), source = 'WHO') |>
     select(-InciVar) |> 
     # add UHC index by year and location
     left_join(uhc_data, by = c('location_name' = 'GEO_NAME_SHORT', 'Year' = 'Year'))

gbd_long <- DataAll |>
     select(location_id, Location_ID, location_name, starts_with('GBD_Inci_'), VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant) |>
     pivot_longer(cols = matches('GBD_Inci_\\d{4}'), names_to = 'InciVar', values_to = 'Incidence') |>
     mutate(Year = as.integer(str_extract(InciVar, '\\d{4}')), source = 'GBD') |>
     select(-InciVar) |> 
     # add UHC index by year and location
     left_join(uhc_data, by = c('location_name' = 'GEO_NAME_SHORT', 'Year' = 'Year'))

population_long <- DataAll |>
     select(location_id, Location_ID, starts_with('Population_')) |>
     pivot_longer(cols = matches('Population_\\d{4}'), names_to = 'PopVar', values_to = 'Population') |>
     mutate(Year = as.integer(str_extract(PopVar, '\\d{4}'))) |>
     select(-PopVar)

DataLongCounts <- bind_rows(who_long, gbd_long) |>
     left_join(population_long, by = c('location_id', 'Location_ID', 'Year')) |>
     # keep only years of interest and ensure predictors are factors
     filter(Year %in% c(2019, 2021, 2024)) |>
     mutate(
          VaccineDose = factor(VaccineDose, levels = levels(DataAll$VaccineDose)),
          TimeFirstShotG = factor(TimeFirstShotG, levels = levels(DataAll$TimeFirstShotG)),
          TimeLastShotG = factor(TimeLastShotG, levels = levels(DataAll$TimeLastShotG)),
          VaccinePregnant = factor(VaccinePregnant, levels = rev(levels(DataAll$VaccinePregnant)))
     ) |>
     filter(!is.na(Incidence) & !is.na(Population) & Population > 0)

# Prepare a log-transformed incidence outcome for Bayesian hierarchical Gaussian models
# we add a small offset to Cases to avoid log(0); compute rate per person and take log
DataLongRates <- DataLongCounts |>
     mutate(
          rate = (Incidence * Population + 0.5) / Population,  # incidence per person with 0.5 offset
          log_rate = log(rate)
     ) |>
     # ensure region is available and use as factor; Region was joined into DataAll earlier
     left_join(DataAll |> select(location_id, Location_ID, Region), by = c('location_id', 'Location_ID')) |>
     mutate(
          Region = ifelse(is.na(Region), 'Unknown', Region),
          Region = factor(Region),
          uhc_cat = ifelse(is.na(uhc_cat), 'Unknown', uhc_cat),
          uhc_cat = factor(uhc_cat, levels = c('Very low', 'Low', 'Medium', 'High', 'Very high', 'Unknown')),
          location_id = factor(location_id)
     )

# Fit ---------------------------------------------------------------------

car::vif(lm(log_rate ~ VaccineDose + TimeFirstShotG + TimeLastShotG +
                 VaccinePregnant + uhc_index, data = DataLongRates))

# Bayesian hierarchical Gaussian model (log-rate outcome) using brms
fit_bayes_hier_gaussian <- function(df, year, source = NULL, min_obs = 10) {
     dsub <- df |> filter(Year == year)
     if(!is.null(source)) dsub <- dsub |> filter(source == !!source)
     if(nrow(dsub) < min_obs) { message('Too few obs for ', year, ' source=', source); return(NULL) }

     # covariates: VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant
     # nested random intercepts: (1 | Region/location_id)
     fmla <- as.formula('log_rate ~ VaccineDose + TimeFirstShotG + TimeLastShotG + VaccinePregnant + uhc_cat + (1|Region/location_id)')

     # weakly informative priors as specified
     priors <- c(
          prior(normal(0, 0.35), class = "b"),
          prior(normal(0, 2), class = "Intercept"),
          prior(normal(0, 2), class = "sigma"),
          prior(normal(0, 2), class = "sd")
     )

     # increase cores for sampling
     options(mc.cores = parallel::detectCores())

     fit <- tryCatch({
          brm(formula = fmla,
              data = dsub,
              family = gaussian(),
              prior = priors,
              chains = 4,
              iter = 4000,
              warmup = 1000,
              control = list(adapt_delta = 0.95, max_treedepth = 15),
              seed = 20251029,
              refresh = 0)
     }, error = function(e) { message('brms fit error ', year, ' source=', source, ': ', e$message); return(NULL) })

     if(is.null(fit)) return(NULL)

     # extract posterior draws for fixed effects
     post <- as.data.frame(posterior_samples(fit))
     b_names <- names(post)[grepl('^b_', names(post))]

     summaries <- map_dfr(b_names, function(nm) {
          draws <- post[[nm]]
          tibble(
               term = sub('^b_', '', nm),
               median = median(draws),
               lower = quantile(draws, 0.025),
               upper = quantile(draws, 0.975),
               pr_gt_zero = mean(draws > 0)
          )
     }) |>
     mutate(year = year, source = ifelse(is.null(source), 'both', source))

     list(fit = fit, summary = summaries)
}

# Fit Bayesian hierarchical models for each available source-year combo
combos_bayes <- DataLongRates |>
     group_by(source, Year) |>
     summarise(n = n(), .groups = 'drop')

bayes_models <- list()
bayes_results <- list()
for(r in seq_len(nrow(combos_bayes))) {
     src <- combos_bayes$source[r]
     yr  <- combos_bayes$Year[r]
     lab <- paste(src, yr, sep = '_')
     res <- fit_bayes_hier_gaussian(DataLongRates, as.integer(yr), source = src)
     if(!is.null(res)) {
          bayes_models[[lab]] <- res$fit
          bayes_results[[lab]] <- res$summary
          # save model object
          saveRDS(res$fit, file = file.path('./Output', paste0('brms_hier_gaussian_', lab, '.rds')))
     }
}

## combine summaries and export
results_bayes_all <- imap_dfr(bayes_results, ~ mutate(.x, model = .y)) |>
     separate(model, into = c('source', 'year'), sep = '_', remove = FALSE) |>
     mutate(year = as.integer(year)) |>
     select(model, source, year, term, median, lower, upper, pr_gt_zero)

# also provide rate-ratio scale (exp of log-scale coefficients)
results_bayes_all <- results_bayes_all |>
     mutate(rr_median = exp(median), rr_lower = exp(lower), rr_upper = exp(upper))

# save combined models+summary object for reproducibility
saveRDS(list(models = bayes_models, summary = results_bayes_all),
        file = './Output/bayes_models_and_summary.rds')

# create a forest plot of rate ratios for non-intercept terms
plot_df <- results_bayes_all |> filter(!term %in% c('Intercept')) |> 
     mutate(term = factor(term, levels = unique(term)))

p_forest <- ggplot(plot_df, aes(x = term, y = rr_median, ymin = rr_lower, ymax = rr_upper, color = source)) +
     geom_pointrange(position = position_dodge(width = 0.6), size = 0.6) +
     geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray50') +
     coord_flip() +
     facet_wrap(~ year, scales = 'fixed') +
     scale_y_log10() +
     scale_x_discrete(limits = rev(levels(plot_df$term))) +
     labs(x = NULL, y = 'Rate ratio (exp(coef))', title = 'Posterior rate ratios (median and 95% CI)') +
     theme_minimal() +
     theme(legend.position = 'bottom')

ggsave('./Output/Figure 5.png',
       p_forest,
       width = 10, height = 6, dpi = 300)

# write numeric summary
write.csv(results_bayes_all,
          file = './Output/Figure 5.xlsx')
