#####################################
## @Description: Figure 4 : Global pertussis incidence vs. vaccination program
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-21 19:20:50
## @LastEditors: Li Kangguo
## @LastEditTime: 2025-10-22 11:23:43
#####################################

# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(openxlsx)
library(paletteer)
library(Cairo)

rm(list = ls())

# data ----------------------------------------------------------------

region_names <- c("Global", 'WHO region',
                  "African Region", "Eastern Mediterranean Region", "European Region",
                  "Region of the Americas", "South-East Asia Region", "Western Pacific Region")

## vaccine -----------------------------------------------------------------

# Country iso codes
DataCountry <- read.csv('./Data/ISO_CODE.csv') |> 
     select(-location_name)

DataVaccine <- read.xlsx("./data/VaccineSchedule.xlsx") 

DataVaccine <- DataVaccine |> 
     select(CODE, NAME, VaccinePregnant, VaccineDose, TimeFirstShot, TimeLastShot)

# add Greenland data: https://pmc.ncbi.nlm.nih.gov/articles/PMC7034463/
DataVaccine <- DataVaccine |> 
     rbind(data.frame(CODE = 'GRL',
                      NAME = 'Greenland',
                      VaccinePregnant = 1,
                      VaccineDose = 4,
                      TimeFirstShot = 3,
                      TimeLastShot = 6*12))

# Fix Romania CODE: ROU -> ROM
DataVaccine$CODE[DataVaccine$CODE == 'ROU'] <- 'ROM'

DataVaccine <- DataVaccine |> 
     mutate(TimeFirstShotG = cut(as.numeric(TimeFirstShot), 
                                 breaks = c(-Inf, 1.5, 2.5, Inf),
                                 labels = c('[1.0,1.5)m', 
                                            '[1.5,2.5)m',
                                            '[2.5,3.0]m'),
                                 include.lowest = TRUE),
            TimeFirstShotG = factor(TimeFirstShotG, 
                                    levels = c('[1.0,1.5)m', 
                                               '[1.5,2.5)m',
                                               '[2.5,3.0]m')),
            TimeLastShotG = cut(as.numeric(TimeLastShot), 
                                breaks = c(-Inf, 12, 30, 66, 150, Inf),
                                labels = c('[3,6]m', 
                                           '[1,2]y',
                                           '[3,5]y',
                                           '[6,12]y',
                                           '[13,18]y'),
                                include.lowest = TRUE)) |> 
     select(-NAME)

## WHO cases ---------------------------------------------------------------

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
df_cases <- read.xlsx('./Data/Pertussis reported cases and incidence.xlsx')
names(df_cases)[1] <- 'Location'

# pivot data to long format
df_cases <- df_cases |>
     filter(!is.na(Disease)) |>
     select(-Disease) |>
     pivot_longer(cols = -Location, names_to = 'Year', values_to = 'Cases') |>
     mutate(
          Year = as.integer(Year),
          Cases = as.numeric(str_replace_all(Cases, ',', '')),
          Location = standardize_location(Location)
     ) |>
     filter(!Location %in% region_names)

# read population data from UN data
df_population <- read.csv('./Data/unpopulation_dataportal_20251021120614.csv')

# pivot data to long format
df_population <- df_population |>
     select(
          Location = Category,
          Location_ID = Estimate_Type,
          Year = Estimate_Method,
          Population = Subregion
     ) |>
     mutate(Location = standardize_location(Location))

# check data complete
unique(df_cases$Location)[
     !unique(df_cases$Location) %in% unique(df_population$Location)
]

unique(df_population$Location)[
     !unique(df_population$Location) %in% unique(df_cases$Location)
]

DataIncidenceWHO <- df_cases |>
     group_by(Location, Year) |>
     summarise(Cases = sum(Cases), .groups = 'drop') |>
     left_join(df_population, by = c('Location', 'Year')) |>
     mutate(Incidence = (Cases / Population) * 1e5) |> 
     filter(Year %in% c(2019, 2021, 2024), !is.na(Cases)) |> 
     select(Location_ID, Location, Year, Incidence) |>
     pivot_wider(names_from = Year, values_from = Incidence, names_prefix = 'WHO_Inci_')

# also keep raw Cases and Population for the target years in wide format so downstream
# scripts can fit count models with offsets without re-loading raw files
DataCasesPopWide <- df_cases |>
     group_by(Location, Year) |>
     summarise(Cases = sum(Cases), .groups = 'drop') |>
     left_join(df_population, by = c('Location', 'Year')) |>
     filter(Year %in% c(2019, 2021, 2024)) |>
     select(Location_ID, Location, Year, Cases, Population) |>
     pivot_wider(names_from = Year, values_from = c(Cases, Population), names_sep = "_")

rm(df_cases, df_population, case_to_pop_map, standardize_location, decode_html_entities)

## GBD cases ---------------------------------------------------------------

# read GBD incidence data
df_incidence_gbd_country <- read.csv('./Data/IHME-GBD_2021_DATA-7facc03b-1.csv')

df_incidence_gbd_country <- df_incidence_gbd_country|>
     # drop regional data
     filter(!location_name %in% region_names, measure_name != "Deaths",
            age_name == 'Age-standardized', metric_name == 'Rate',
            year %in% c(2019, 2021, 2024)) |>
     select(location_id, location_name, year, val)

DataIncidenceGBD <- df_incidence_gbd_country |>
     pivot_wider(names_from = year, values_from = val, names_prefix = 'GBD_Inci_')

rm(df_incidence_gbd_country)

## merge data ---------------------------------------------------------------

DataInciRaw <- DataIncidenceWHO |> 
     full_join(DataCountry, by = c('Location_ID' = 'ISO3')) |> 
     full_join(DataIncidenceGBD, by = 'location_id') |> 
     # join the raw counts/population (wide) so DataAll contains Cases_2019, Population_2019, etc.
     left_join(DataCasesPopWide, by = c('Location', 'Location_ID'))

# add GBD data
DataAll <- DataInciRaw |> 
     full_join(DataVaccine, by = c('Location_ID' = 'CODE')) |> 
     mutate(TimeLastShotG = factor(TimeLastShotG, 
                                   levels = c('[3,6]m', 
                                              '[1,2]y',
                                              '[3,5]y',
                                              '[6,12]y',
                                              '[13,18]y')),
            VaccineDose = as.factor(VaccineDose),
            VaccinePregnant = factor(VaccinePregnant, 
                                     levels = c(1, 0),
                                     labels = c('Recommended', 'Not recommended')))

# visualization -----------------------------------------------------------

visual_vars <- c('VaccineDose', 'TimeFirstShotG', 'TimeLastShotG', 'VaccinePregnant')
lab_vars <- c('Vaccine Doses', 
              'Time to First Shot', 
              'Time to Last Shot', 
              'Maternal Vaccination')

visual_vars_y <- c('WHO_Inci_2019', 'WHO_Inci_2021', 'WHO_Inci_2024', 'GBD_Inci_2019', 'GBD_Inci_2021')
labs_vars_y <- c('Incidence rate, 2019 (WHO)', 
                 'Incidence rate, 2021 (WHO)', 
                 'Incidence rate, 2024 (WHO)', 
                 'Incidence rate, 2019 (GBD)', 
                 'Incidence rate, 2021 (GBD)')

scientific_10 <- function(x) { 
     parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) 
}

y <- 'WHO_Inci_2019'; i <- 1; start_index <- 0

fill_color <- c("#2A6EBB", "#F0AB00", "#C50084", "#7D5CC6", "#E37222", "#69BE28", "#00B2A9", "#CD202C", "#747678")
names(fill_color) <- labs_vars_y

# number of bootstrap resamples used for CIs and for tests (make configurable)
bootstrap_n <- 1000
# number of permutations for permutation test
permutation_n <- 5000

compute_panel <- function(y, i) {
     x <- visual_vars[i]
     type <- labs_vars_y[visual_vars_y == y]

     DataAllTemp <- DataAll |> 
          select(location_id, !!x, !!y) |> 
          rename(var = !!x, y = !!y) |>
          drop_na() |> 
          mutate(type = type,
                 type =  factor(type, levels = labs_vars_y))

     # Generate neighbor pairwise comparisons for unique levels of var
     unique_vars <- levels(DataAllTemp$var)

     # Drop not enough level
     unique_vars <- unique_vars[unique_vars %in% unique(DataAllTemp$var)]

     # generate all pairwise comparisons between levels of var
     if(length(unique_vars) > 1) {
          my_comparisons <- combn(unique_vars, 2, simplify = FALSE)
     } else {
          my_comparisons <- list()
     }

     # Compute bootstrap mean and 95% CI per group (bootstrap by resampling with replacement)
     bootstrap_summary <- DataAllTemp |>
          group_by(var, type) |> 
          summarise(n = n(), mean = mean(y, na.rm = TRUE), .groups = 'drop') |> 
          rowwise() |> 
          mutate(
               ci = list({
                    vals <- DataAllTemp$y[DataAllTemp$var == var & DataAllTemp$type == type]
                    if(length(vals) <= 1) {
                         c(lower = mean, upper = mean)
                    } else {
                         boots <- replicate(bootstrap_n, mean(sample(vals, replace = TRUE), na.rm = TRUE))
                         quantile(boots, probs = c(0.025, 0.975), na.rm = TRUE)
                    }
               })
          ) |> 
          mutate(ci_lower = ci[[1]], ci_upper = ci[[2]]) |> 
          select(-ci, -n) |> 
          ungroup()

     # Compute permutation-based pairwise tests using an inline permutation test implementation
     if(length(my_comparisons) > 0) {
          perm_test_pair <- function(vec, grp, n_perm = 5000, seed = NULL) {
               if(!is.null(seed)) set.seed(seed)
               levs <- levels(factor(grp))
               if(length(levs) != 2) return(NA_real_)
               obs <- mean(vec[grp == levs[2]], na.rm = TRUE) - mean(vec[grp == levs[1]], na.rm = TRUE)
               combined <- vec
               n <- length(vec)
               perm_diffs <- numeric(n_perm)
               for(b in seq_len(n_perm)) {
                    perm_grp <- sample(grp, size = n, replace = FALSE)
                    perm_diffs[b] <- mean(combined[perm_grp == levs[2]], na.rm = TRUE) - mean(combined[perm_grp == levs[1]], na.rm = TRUE)
               }
               mean(abs(perm_diffs) >= abs(obs), na.rm = TRUE)
          }

          p_list <- lapply(my_comparisons, function(pair) {
               g1 <- pair[1]; g2 <- pair[2]
               sub <- DataAllTemp |> filter(var %in% c(g1, g2))
               if(nrow(sub) < 2) return(tibble(group1 = g1, group2 = g2, p = NA_real_))
               pval <- tryCatch(perm_test_pair(sub$y, sub$var, n_perm = permutation_n), error = function(e) NA_real_)
               tibble(group1 = g1, group2 = g2, p = pval)
          })
          p_df <- bind_rows(p_list)

          # adjust and format
          p_df <- p_df |> mutate(p.adj = p.adjust(p, method = 'holm')) |>
               mutate(p.label = ifelse(is.na(p.adj), NA_character_, formatC(p.adj, format = 'f', digits = 2)),
                      p.signif = case_when(is.na(p.adj) ~ '',
                                           p.adj < 0.001 ~ '***',
                                           p.adj < 0.01 ~ '**',
                                           p.adj < 0.05 ~ '*',
                                           TRUE ~ 'ns'))

          # compute plotting positions using the global maximum of bootstrap CI upper bounds
          max_ci_upper <- max(bootstrap_summary$ci_upper, na.rm = TRUE)
          if(is.infinite(max_ci_upper) || is.na(max_ci_upper)) {
               max_ci_upper <- max(DataAllTemp$y, na.rm = TRUE)
          }
          spacing <- 0.25 * abs(max_ci_upper)

          p_df <- p_df |>
               mutate(xmin = group1, xmax = group2) |>
               arrange(group1, group2) |>
               mutate(y.position = max_ci_upper + spacing * (row_number() - 1), p.value = p)
     } else {
          p_df <- NULL
     }

     list(
          DataAllTemp = DataAllTemp,
          bootstrap_summary = bootstrap_summary,
          p_df = p_df,
          unique_vars = unique_vars,
          type = type
     )
}

# Compute all panels ahead of plotting and store results in a list
compute_all_panels <- function() {
     panels <- vector('list', length(visual_vars_y))
     names(panels) <- visual_vars_y
     for(j in seq_along(visual_vars_y)) {
          panels[[j]] <- vector('list', length(visual_vars))
          for(i in seq_along(visual_vars)) {
               panels[[j]][[i]] <- compute_panel(visual_vars_y[j], i)
          }
     }
     panels
}

# Precompute all panels
panels_list <- compute_all_panels()

# Plot from a precomputed panel result
plot_from_res <- function(res, i, start_index) {
     DataAllTemp <- res$DataAllTemp
     bootstrap_summary <- res$bootstrap_summary
     p_df <- res$p_df
     unique_vars <- res$unique_vars
     type <- res$type

     ggplot(DataAllTemp, aes(x = var, y = y, color = type)) +
          geom_pointrange(data = bootstrap_summary,
                          aes(x = var, y = mean, ymin = ci_lower, ymax = ci_upper, color = type),
                          position = position_dodge(width = 0.6),
                          size = 0.8,
                          fatten = 1.5,
                          show.legend = TRUE) +
          {if(!is.null(p_df)) ggpubr::stat_pvalue_manual(p_df, label = "p.signif", hide.ns = FALSE) else NULL} +
          scale_color_manual(values = fill_color,
                             breaks = labs_vars_y,
                             drop = FALSE) +
          scale_y_continuous(limits = c(0, NA),
                             expand = expansion(mult = c(0, 0.25))) +
          scale_x_discrete(limits = unique_vars) +
          theme_bw() +
          theme(panel.grid = element_blank()) +
          labs(title = LETTERS[start_index + i],
               x = lab_vars[i],
               color = NULL,
               y = "Incidence rate (per 100,000 population)")
}

# Save figures (plot from precomputed results)
fig_1 <- lapply(seq_along(visual_vars), function(i) plot_from_res(panels_list[[1]][[i]], i, 0))
fig_2 <- lapply(seq_along(visual_vars), function(i) plot_from_res(panels_list[[2]][[i]], i, 4))
fig_3 <- lapply(seq_along(visual_vars), function(i) plot_from_res(panels_list[[3]][[i]], i, 8))
fig_4 <- lapply(seq_along(visual_vars), function(i) plot_from_res(panels_list[[4]][[i]], i, 12))
fig_5 <- lapply(seq_along(visual_vars), function(i) plot_from_res(panels_list[[5]][[i]], i, 16))

figs <- c(fig_1, fig_2, fig_3, fig_4, fig_5)

fig <- wrap_plots(figs, ncol = 4, guides = 'collect', axis_titles = 'collect') & 
     theme(legend.position = 'bottom',)

# save
ggsave('./Output/Figure 4.png',
       plot = fig,
       width = 12, height = 12, dpi = 300)

ggsave('./Output/Figure 4.pdf',
       plot = fig,
       # plot = p1 / p2,
       width = 12, height = 12, 
       device = cairo_pdf,
       family = "Helvetica")

# save figure data --------------------------------------------------------

# Save computed stats to Excel: one sheet per panel/variable (avoid overlapping writes)
wb <- createWorkbook()
for(j in seq_along(panels_list)) {
     for(i in seq_along(panels_list[[j]])) {
          res <- panels_list[[j]][[i]]
          # create a safe sheet name: panel letter + short descriptor
          panel_letter <- LETTERS[j]
          type_clean <- gsub('[^A-Za-z0-9]', '_', as.character(res$type))
          var_clean <- gsub('[^A-Za-z0-9]', '_', visual_vars[i])
          sheet_name <- paste0(panel_letter, '_', substr(type_clean, 1, 20), '_', substr(var_clean, 1, 10))
          # Excel sheet names must be <= 31 chars
          sheet_name <- substr(sheet_name, 1, 31)
          # ensure sheet name is unique in workbook
          existing <- names(wb)
          suffix <- 1
          base_name <- sheet_name
          while(sheet_name %in% existing) {
               suffix <- suffix + 1
               # keep within 31 chars
               sheet_name <- substr(paste0(base_name, '_', suffix), 1, 31)
          }
          addWorksheet(wb, sheet_name)

          # write bootstrap_summary if present
          if(!is.null(res$bootstrap_summary) && nrow(res$bootstrap_summary) > 0) {
               writeData(wb, sheet = sheet_name, x = res$bootstrap_summary,
                          startCol = 1, startRow = 1, colNames = TRUE, rowNames = FALSE)
               current_row <- nrow(res$bootstrap_summary) + 3
          } else {
               current_row <- 1
          }

          # write p_df if present
          if(!is.null(res$p_df) && nrow(res$p_df) > 0) {
               writeData(wb, sheet = sheet_name, x = res$p_df,
                          startCol = 1, startRow = current_row, colNames = TRUE, rowNames = FALSE)
               current_row <- current_row + nrow(res$p_df) + 3
          }

          # write DataAllTemp if present
          if(!is.null(res$DataAllTemp) && nrow(res$DataAllTemp) > 0) {
               writeData(wb, sheet = sheet_name, x = res$DataAllTemp,
                          startCol = 1, startRow = current_row, colNames = TRUE, rowNames = FALSE)
          }
     }
}

saveWorkbook(wb,
             './Output/Figure 4.xlsx',
             overwrite = TRUE)

# save DataAll for figure5
save(DataAll, file = './Output/DataAll.RData')
