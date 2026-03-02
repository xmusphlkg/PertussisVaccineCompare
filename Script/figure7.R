#####################################
## @Description: Figure 7 - Bayesian hierarchical analysis for pertussis incidence vs vaccination program
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-22 11:08:31
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-02 11:43:08
#####################################

library(tidyverse)
library(gtsummary)
library(broom.mixed)
library(openxlsx)
library(flextable)
library(gt)
library(brms)
library(patchwork)

rm(list = ls())

# data --------------------------------------------------------------------

region_names <- c("Global", 'WHO region',
                  "African Region", "Eastern Mediterranean Region", "European Region",
                  "Region of the Americas", "South-East Asia Region", "Western Pacific Region")

## Vaccine data ------------------------------------------------------------

load('./Output/DataAll.RData')

DataAll <- DataAll |>
     select(location_id, location_name, Location, Location_ID, ISO2, 
            VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant,
            contains('WHO_Inci_'), contains('GBD_Inci_'), contains('Population_'))

## DTP coverage ------------------------------------------------------------

DataCoverage <- read.xlsx('./Data/Diphtheria tetanus toxoid and pertussis (DTP) vaccination coverage 2026-25-02 14-33 UTC.xlsx')

DataDTP3 <- DataCoverage |>
     filter(ANTIGEN %in% c("DTPCV1", "DTPCV3"), 
            COVERAGE_CATEGORY == 'WUENIC', 
            GROUP == 'COUNTRIES') |> 
     transmute(Location_ID = CODE,
               Year = as.integer(YEAR),
               Antigen = ANTIGEN,
               Coverage = as.numeric(COVERAGE)) |>
     distinct(Location_ID, Year, Antigen, .keep_all = TRUE) |>
     mutate(Antigen = recode(Antigen,
                             "DTPCV1" = "DTP1",
                             "DTPCV3" = "DTP3")) |>
     pivot_wider(names_from = Antigen, values_from = Coverage)

## check all locations in DataAll have DTP coverage data for 2019, 2021, 2024
missing_dtp <- DataAll |>
     select(Location_ID, Location) |>
     distinct() |>
     anti_join(DataDTP3 |> filter(Year %in% c(2019, 2021, 2024)), by = c('Location_ID' = 'Location_ID'))

## population data ---------------------------------------------------------

# attach Region mapping: https://cdn.who.int/media/docs/default-source/air-pollution-documents/air-quality-and-health/un-agencies-region-classification-for-country.xlsx?sfvrsn=289af35f_3
unpop <- read.xlsx("./Data/CLASS_2025_07_02.xlsx")
unpop_map <- unpop |>
     select(Code, Income.group) |>
     rename(ISO3 = Code, Region = Income.group)

DataAll <- DataAll |>
     left_join(unpop_map, by = c(Location_ID = 'ISO3'))

rm(unpop, unpop_map)

## UHC data ----------------------------------------------------------------

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
     # add DTP coverage (DTP1, DTP3) by country-year
     left_join(DataDTP3 |> select(Location_ID, Year, DTP1, DTP3),
               by = c('Location_ID', 'Year')) |>
     mutate(
          Region = ifelse(is.na(Region), 'Unknown', Region),
          Region = factor(Region),
          uhc_cat = ifelse(is.na(uhc_cat), 'Unknown', uhc_cat),
          # scale DTP3 coverage per 10-percentage-point increase for interpretability
          DTP3_10 = ifelse(is.na(DTP3), NA_real_, DTP3 / 10)
     )

# Fit ---------------------------------------------------------------------

car::vif(lm(log_rate ~ VaccineDose + TimeFirstShotG + TimeLastShotG +
                 VaccinePregnant + uhc_index, data = DataLongRates))

# Bayesian hierarchical Gaussian model (log-rate outcome) using brms
fit_bayes_hier_gaussian <- function(df, year, source = NULL, min_obs = 10) {
     dsub <- df |> filter(Year == year)
     if(!is.null(source)) dsub <- dsub |> filter(source == !!source)
     # keep only observations with available DTP3 coverage
     dsub <- dsub |> filter(!is.na(DTP3_10))
     if(nrow(dsub) < min_obs) { message('Too few obs for ', year, ' source=', source); return(NULL) }

     # covariates: VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant, DTP3_10 (per 10%-point)
     # nested random intercepts: (1 | uhc_cat/location_id)
     fmla <- as.formula('log_rate ~ VaccineDose + TimeFirstShotG + TimeLastShotG + VaccinePregnant + DTP3_10 + (1|uhc_cat/location_id)')

     # weakly informative priors as specified
     priors <- c(
          prior(normal(0, 2), class = "b"),
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
     summarise(n = n(), .groups = 'drop') |> 
     arrange(desc(source), Year)

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
     mutate(group = case_when(grepl('^VaccineDose', term) ~ 'Number of scheduled doses',
                              grepl('^TimeFirstShotG', term) ~ 'Timing of first scheduled dose',
                              grepl('^TimeLastShotG', term) ~ 'Timing of last scheduled childhood dose',
                              grepl('^VaccinePregnant', term) ~ 'Maternal pertussis vaccination policy',
                              term == 'DTP3_10' ~ 'DTP3 coverage (per 10 percentage points)',
                              TRUE ~ term),
            trem_clean = case_when(grepl('^VaccineDose', term) ~ sub('^VaccineDose', '', term),
                                   term == 'TimeFirstShotG1.52.5m' ~ '[1.5,2.5)m',
                                   term == 'TimeFirstShotG2.53.0m' ~ '[2.5,3.0]m',
                                   term == 'TimeLastShotG12y' ~ '[1,2]y',
                                   term == 'TimeLastShotG35y' ~ '[3,5]y',
                                   term == 'TimeLastShotG612y' ~ '[6,12]y',
                                   term == 'TimeLastShotG1318y' ~ '[13,18]y',
                                   term == 'VaccinePregnantRecommended' ~ 'Recommended',
                                   TRUE ~ term))
            

labs_vars_y <- c('Incidence rate, 2019 (WHO)', 
                 'Incidence rate, 2021 (WHO)', 
                 'Incidence rate, 2024 (WHO)', 
                 'Incidence rate, 2019 (GBD)', 
                 'Incidence rate, 2021 (GBD)')

fill_color <- c("#2A6EBB", "#F0AB00", "#C50084", "#7D5CC6", "#E37222")
names(fill_color) <- labs_vars_y

plot_function <- function(i){
     outcome <- plot_df |> 
          filter(group == unique(plot_df$group)[i]) |> 
          mutate(source = factor(paste0(source, ' ', year),
                                 levels = c('WHO 2019', 'WHO 2021', 'WHO 2024', 'GBD 2019', 'GBD 2021'),
                                 labels = labs_vars_y))
     
     outcome |> 
          ggplot(aes(x = trem_clean, y = rr_median, ymin = rr_lower, ymax = rr_upper,
                     color = source)) +
          geom_pointrange(position = position_dodge(width = 0.6), show.legend = T) +
          geom_line(aes(group = source), position = position_dodge(width = 0.6), show.legend = F) +
          geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray50') +
          scale_color_manual(values = fill_color,
                             drop = F) +
          labs(x = unique(plot_df$group)[i],
               y = 'Adjusted rate ratio (95% CrI)',
               title = LETTERS[i],
               color = 'Data source')+
          theme_bw() +
          guides(color = guide_legend(ncol = 2))
}

fig <- lapply(seq_len(length(unique(plot_df$group))), plot_function)

fig[[6]] <- patchwork::guide_area()

# increase right margin
fig <- wrap_plots(fig, ncol = 3, guides = 'collect')&
     theme(legend.position = 'right')

ggsave('./Output/Figure 7.png',
       fig,
       width = 13, height = 6, dpi = 300)

ggsave('./Output/Figure 7.pdf',
       plot = fig,
       width = 12, height = 6, 
       device = cairo_pdf,
       family = "Helvetica")

# write numeric summary
write.xlsx(results_bayes_all,
           file = './Output/Figure 7.xlsx')
