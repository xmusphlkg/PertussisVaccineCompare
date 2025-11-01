#####################################
## @Description: Figure 5 - GLMM analysis for pertussis incidence vs vaccination program
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-22 11:08:31
## @LastEditors: Li Kangguo
## @LastEditTime: 2025-11-01 10:32:55
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

## WHO data ----------------------------------------------------------------

decode_html_entities <- function(x) {
     x <- str_replace_all(x, "&#39;", "'")
     x <- str_replace_all(x, "&amp;", "&")
     x
}

# mapping from WHO-case names (after HTML decode) to UN population names
case_to_pop_map <- c(
     "United Kingdom of Great Britain and Northern Ireland" = "United Kingdom",
     "Netherlands (Kingdom of the)" = "Netherlands",
     "Micronesia (Federated States of)" = "Micronesia (Fed. States of)",
     "Democratic Republic of the Congo" = "Dem. Rep. of the Congo",
     "Democratic People's Republic of Korea" = "Dem. People's Rep. of Korea",
     "Lao People's Democratic Republic" = "Lao People's Dem. Republic",
     "Bonaire" = "Bonaire, Sint Eustatius and Saba",
     "Saba" = "Bonaire, Sint Eustatius and Saba",
     "Sint Eustatius" = "Bonaire, Sint Eustatius and Saba",
     "occupied Palestinian territory, including east Jerusalem" = "State of Palestine",
     "Wallis and Futuna" = "Wallis and Futuna Islands",
     "Kosovo (in accordance with UN Security Council resolution 1244 (1999))" = "Kosovo (under UNSC res. 1244)",
     "Côte d'Ivoire" = "Côte d'Ivoire" # kept for clarity after decode
)

standardize_location <- function(x, map = case_to_pop_map) {
     x2 <- decode_html_entities(x)
     x2 <- str_trim(x2)
     # apply explicit mapping when present
     x2 <- ifelse(x2 %in% names(map), unname(map[x2]), x2)
     x2
}

# read cases data from WHO data
DataWHOCases <- read.xlsx('./Data/Pertussis reported cases and incidence.xlsx')
names(DataWHOCases)[1] <- 'Location'

# pivot data to long format
DataWHOCases <- DataWHOCases |>
     filter(!is.na(Disease)) |>
     select(-Disease) |>
     pivot_longer(cols = -Location, names_to = 'Year', values_to = 'Cases') |>
     mutate(Year = as.integer(Year),
            Cases = as.numeric(str_replace_all(Cases, ',', '')),
            Location = standardize_location(Location)) |>
     group_by(Location, Year) |>
     summarise(Cases = sum(Cases), .groups = 'drop') |>
     filter(!Location %in% region_names, Year %in% c(2019, 2021, 2024)) |> 
     pivot_wider(names_from = Year, values_from = Cases, names_prefix = 'WHO_Cases_')

## GBD data ---------------------------------------------------------------

# read GBD incidence data
DataGBDCases <- read.csv('./Data/IHME-GBD_2021_DATA-7facc03b-1.csv')

DataGBDCases <- DataGBDCases|>
     # drop regional data
     filter(!location_name %in% region_names, measure_name != "Deaths",
            age_name == 'All ages', metric_name == 'Number',
            year %in% c(2019, 2021, 2024)) |>
     select(location_id, location_name, year, val)

DataGBDCases <- DataGBDCases |>
     pivot_wider(names_from = year, values_from = val, names_prefix = 'GBD_Cases_')

rm(region_names, case_to_pop_map, decode_html_entities, standardize_location)

## bind --------------------------------------------------------------------

DataAll <- DataAll |>
     select(location_id, location_name, Location, Location_ID, ISO2, Population_2019, Population_2021,
            Population_2024, VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant) |>
     left_join(DataWHOCases, by = c('location_name' = 'Location')) |>
     left_join(DataGBDCases, by = c('location_id', 'location_name'))

# attach Region mapping: https://cdn.who.int/media/docs/default-source/air-pollution-documents/air-quality-and-health/un-agencies-region-classification-for-country.xlsx?sfvrsn=289af35f_3
unpop <- read.xlsx("./Data/un-agencies-region-classification-for-country.xlsx")
unpop_map <- unpop |>
     select(ISO.Country.code, WHO.Region.name2) |>
     rename(ISO3 = ISO.Country.code, Region = WHO.Region.name2)

DataAll <- DataAll |>
     left_join(unpop_map, by = c(Location_ID = 'ISO3'))

# model -------------------------------------------------------------------

## Build long dataset stacking WHO and GBD cases with a `source` column
# pivot WHO and GBD case columns separately, then bind_rows into a long table
who_long <- DataAll |>
     select(location_id, starts_with('WHO_Cases_'), VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant, starts_with('Population_')) |>
     pivot_longer(cols = matches('WHO_Cases_\\d{4}'), names_to = 'CasesVar', values_to = 'Cases') |>
     mutate(Year = as.integer(str_extract(CasesVar, '\\d{4}')), source = 'WHO') |>
     select(-CasesVar)

gbd_long <- DataAll |>
     select(location_id, starts_with('GBD_Cases_'), VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant, starts_with('Population_')) |>
     pivot_longer(cols = matches('GBD_Cases_\\d{4}'), names_to = 'CasesVar', values_to = 'Cases') |>
     mutate(Year = as.integer(str_extract(CasesVar, '\\d{4}')), source = 'GBD') |>
     select(-CasesVar)

DataLongCounts <- bind_rows(who_long, gbd_long) |>
     # attach population by year
     mutate(PopVar = paste0('Population_', Year)) |>
     rowwise() |>
     mutate(Population = as.numeric(cur_data_all()[[PopVar]])) |>
     ungroup() |>
     select(-PopVar) |>
     # keep only years of interest and ensure predictors are factors
     filter(Year %in% c(2019, 2021, 2024)) |>
     mutate(
          VaccineDose = factor(VaccineDose),
          TimeFirstShotG = factor(TimeFirstShotG),
          TimeLastShotG = factor(TimeLastShotG),
          VaccinePregnant = factor(VaccinePregnant),
          source = factor(source)
     ) |>
     filter(!is.na(Cases) & !is.na(Population) & Population > 0)

# Prepare a log-transformed incidence outcome for Bayesian hierarchical Gaussian models
# we add a small offset to Cases to avoid log(0); compute rate per person and take log
DataLongRates <- DataLongCounts |>
     mutate(
          Cases_adj = Cases + 0.5,
          rate = Cases_adj / Population,
          log_rate = log(rate)
     ) |>
     # ensure region is available and use as factor; Region was joined into DataAll earlier
     left_join(DataAll |> select(location_id, Region), by = 'location_id') |>
     mutate(
          Region = ifelse(is.na(Region), 'Unknown', Region),
          Region = factor(Region),
          location_id = factor(location_id)
     )

# Bayesian hierarchical Gaussian model (log-rate outcome) using brms
# formula: log_rate ~ covariates + (1 | Region/location_id)
fit_bayes_hier_gaussian <- function(df, year, source = NULL, min_obs = 10) {
     dsub <- df |> filter(Year == year)
     if(!is.null(source)) dsub <- dsub |> filter(source == !!source)
     if(nrow(dsub) < min_obs) { message('Too few obs for ', year, ' source=', source); return(NULL) }

     # covariates: VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant
     # nested random intercepts: (1 | Region/location_id)
     fmla <- as.formula('log_rate ~ VaccineDose + TimeFirstShotG + TimeLastShotG + VaccinePregnant + (1 | Region/location_id)')

     # weakly informative priors as specified
     priors <- c(
          prior(normal(0, 2), class = 'b'),
          prior(normal(0, 2), class = 'Intercept'),
          prior(normal(0, 2), class = 'sigma'),
          prior(normal(0, 2), class = 'sd')
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

# write numeric summary
write.csv(results_bayes_all, file = './Output/Results_Bayes_Hier_Gaussian.csv')

# save combined models+summary object for reproducibility
saveRDS(list(models = bayes_models, summary = results_bayes_all),
        file = './Output/bayes_models_and_summary.rds')

# create a forest plot of rate ratios for non-intercept terms
plot_df <- results_bayes_all |> filter(!term %in% c('Intercept')) |> 
     mutate(term = factor(term, levels = unique(term)))

if(nrow(plot_df) > 0) {
     p_forest <- ggplot(plot_df, aes(x = term, y = rr_median, ymin = rr_lower, ymax = rr_upper, color = source)) +
          geom_pointrange(position = position_dodge(width = 0.6), size = 0.6) +
          geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray50') +
          coord_flip() +
          facet_wrap(~ year, scales = 'free_y') +
          scale_y_log10() +
          labs(x = NULL, y = 'Rate ratio (exp(coef))', title = 'Posterior rate ratios (median and 95% CI)') +
          theme_minimal() +
          theme(legend.position = 'bottom')
     
     ggsave('./Output/forest_fixed_effects_bayes.png', p_forest, width = 10, height = 6, dpi = 300)
}
