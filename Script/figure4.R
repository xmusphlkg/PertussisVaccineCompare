
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
     full_join(DataIncidenceGBD, by = 'location_id')

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
              'Time to First Shot (months)', 
              'Time to Last Shot (months)', 
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

x <- 'VaccineDose'; y <- 'WHO_Inci_2019'; i <- 1; start_index <- 0

fill_color <- paletteer_d("werpals::pan")
names(fill_color) <- labs_vars_y

plot_panel <- function(y, i, start_index) {
     
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
     
     my_comparisons <- Map(c, unique_vars[-length(unique_vars)], unique_vars[-1])
     
     # print stat compare
     DataAllTemp |> 
          group_by(var) |> 
          summarise(inci = mean(y),
                    .ungroup = T) |> 
          mutate(inci = format(round(inci, 2))) |> 
          print()
     
     ggplot(DataAllTemp,
            aes(x = var, y = y, color = type))+
          geom_boxplot(show.legend = T) +
          geom_jitter(width = 0.2, alpha = 0.5, show.legend = T) +
          scale_color_manual(values = fill_color,
                             breaks = labs_vars_y,
                             drop = FALSE) +
          # Add statistical comparison
          ggpubr::stat_compare_means(comparisons = my_comparisons,
                                     label = "p.signif",
                                     method = "t.test",
                                     p.adjust.method = "holm") +
          scale_y_continuous(trans = 'sqrt',
                             expand = expansion(mult = c(0, 0.25)))+
          scale_x_discrete(limits = unique_vars) +
          theme_bw()+
          theme(panel.grid = element_blank()) +
          labs(title = LETTERS[start_index + i],
               x = lab_vars[i],
               color = NULL,
               y = "Incidence rate (per 100,000 population)")
     
}

# Save figures
fig_1 <- lapply(1:length(visual_vars), plot_panel, y = visual_vars_y[1], start_index = 0)
fig_2 <- lapply(1:length(visual_vars), plot_panel, y = visual_vars_y[2], start_index = 4)
fig_3 <- lapply(1:length(visual_vars), plot_panel, y = visual_vars_y[3], start_index = 8)
fig_4 <- lapply(1:length(visual_vars), plot_panel, y = visual_vars_y[4], start_index = 12)
fig_5 <- lapply(1:length(visual_vars), plot_panel, y = visual_vars_y[5], start_index = 16)

figs <- c(fig_1, fig_2, fig_3, fig_4, fig_5)

fig <- wrap_plots(figs, ncol = 4, guides = 'collect', axis_titles = 'collect') & 
     theme(legend.position = 'bottom',)

# save
ggsave('./Output/Figure 4.png',
       plot = fig,
       width = 12, height = 9, dpi = 300)

ggsave('./Output/Figure 4.pdf',
       plot = fig,
       # plot = p1 / p2,
       width = 12, height = 9, 
       device = cairo_pdf,
       family = "Helvetica")

# save figure data --------------------------------------------------------

write.xlsx(DataAll,
           './Output/Figure 4.xlsx')
