#####################################
## @Description: Figure 1 for pertussis incidence trend analysis
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-22 09:01:50
## @LastEditors: Li Kangguo
## @LastEditTime: 2025-11-01 11:29:24
#####################################

# packages ----------------------------------------------------------------

library(tidyverse)
library(openxlsx)
library(stringr)
library(patchwork)
library(paletteer)


library(nih.joinpoint)

rm(list = ls())

decode_html_entities <- function(x) {
     x <- str_replace_all(x, "&#39;", "'")
     x <- str_replace_all(x, "&amp;", "&")
     x
}

source('./script/joinpoint_setting.R')

fill_color_source <- c("#2A6EBB", "#F0AB00", "#C50084", "#7D5CC6", "#E37222", "#69BE28", "#00B2A9", "#CD202C", "#747678")[1:2]
fill_color_heatmap <- paletteer_d("awtools::a_palette")

# data --------------------------------------------------------------------

## WHO --------------------------------------------------------------------

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

region_names <- c(
     "African Region", "Eastern Mediterranean Region", "European Region",
     "Region of the Americas", "South-East Asia Region", "Western Pacific Region",
     "Global", 'WHO'
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
          Location_ID = Category_ID,
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

df_incidence_who_country <- df_cases |>
     left_join(df_population, by = c('Location', 'Year')) |>
     mutate(Incidence = (Cases / Population) * 1e5)

rm(df_cases, df_population, case_to_pop_map, standardize_location, decode_html_entities)

# read regional data from WHO data
df_incidence_who_region <- read.xlsx('./Data/Pertussis reported cases and incidence region.xlsx')
names(df_incidence_who_region)[1] <- 'Region'

df_incidence_who_region <- df_incidence_who_region |>
     filter(!is.na(Disease)) |>
     select(-Disease, -Denominator) |>
     pivot_longer(cols = -Region, names_to = 'Year', values_to = 'Incidence') |>
     mutate(Year = as.integer(Year),
            Source = 'WHO',
            Incidence = as.numeric(str_replace_all(Incidence, ',', '')),
            lower = NA,
            upper = NA) |> 
     rename(Location = Region) |> 
     arrange(Location, Year)

## GBD ---------------------------------------------------------------------

# read GBD incidence data
df_incidence_gbd_country <- read.csv('./Data/IHME-GBD_2021_DATA-7facc03b-1.csv')

df_incidence_gbd_country <- df_incidence_gbd_country |>
     select(location_id, location_name, year, measure_name, age_name, metric_name, val, lower, upper) |>
     # drop regional data
     filter(!location_name %in% region_names)

# read GBD regional data
df_incidence_gbd_region <- read.csv('./Data/IHME-GBD_2021_DATA-629ff673-1.csv')

df_incidence_gbd_region <- df_incidence_gbd_region |> 
     filter(!location_name %in% "WHO region", measure_name != "Deaths", age_name == 'Age-standardized', metric_name == 'Rate') |>
     select(Location = location_name, Year = year, Incidence = val, lower, upper) |> 
     mutate(Source = 'GBD') |> 
     arrange(Location, Year)

## bind data ---------------------------------------------------------------

df_incidence_region <- rbind(df_incidence_who_region, df_incidence_gbd_region) |> 
     filter(Year >= 2000)

region_names <- c("Global",
                  "African Region", "Eastern Mediterranean Region", "European Region",
                  "Region of the Americas", "South-East Asia Region", "Western Pacific Region")

# JP analysis -------------------------------------------------------------

## WHO ---------------------------------------------------------------------

## build joinpoint model for WHO data
model_incidence_who_global <- joinpoint(df_incidence_who_region |> filter(Location == 'Global', Year %in% 2000:2019),
                                        Year,
                                        Incidence,
                                        run_opt = run_opt_who,
                                        export_opt = export_opt_who)

## build joinpoint model for WHO data
model_incidence_who_region <- joinpoint(df_incidence_who_region |> filter(Location != 'Global', Year %in% 2000:2019),
                                        Year,
                                        Incidence,
                                        by = Location,
                                        run_opt = run_opt_who,
                                        export_opt = export_opt_who)

## build joinpoint model for GBD data
model_incidence_gbd_global <- joinpoint(df_incidence_gbd_region |> filter(Location == 'Global', Year %in% 2000:2019),
                                        Year,
                                        Incidence,
                                        run_opt = run_opt_gbd,
                                        export_opt = export_opt_gbd)

## build joinpoint model for GBD data
model_incidence_gbd_region <- joinpoint(df_incidence_gbd_region |> filter(Location != 'Global', Year %in% 2000:2019),
                                        Year,
                                        Incidence,
                                        by = Location,
                                        run_opt = run_opt_gbd,
                                        export_opt = export_opt_gbd)

# clean data
df_aapc_who <- rbind(get_aapc(model_incidence_who_global) |> mutate(location = 'Global'),
                     get_aapc(model_incidence_who_region)) |>
     mutate(`AAPC (95%CI)` = paste0(aapc, p_value_label),
            aapc = as.numeric(aapc),
            location = factor(location, levels = region_names))

df_aapc_gbd <- rbind(get_aapc(model_incidence_gbd_global) |> mutate(location = 'Global'),
                     get_aapc(model_incidence_gbd_region)) |>
     mutate(`AAPC (95%CI)` = paste0(aapc, p_value_label),
            aapc = as.numeric(aapc),
            location = factor(location, levels = region_names))

save.image('./Output/jp_model.RData')

load('./Output/jp_model.RData')

# figure ------------------------------------------------------------------

## panel line --------------------------------------------------------------

i <- 2

panel_line_function <- function(i){
     data_region <- df_incidence_region |> 
          filter(Location == region_names[i])
     
     breaks <- pretty(c(0, data_region$Incidence, data_region$upper))
     
     if (region_names[i] == 'Global'){
          data_region_jp <- rbind(
               model_incidence_who_global$data_export |> mutate(Source = 'WHO'),
               model_incidence_gbd_global$data_export |> mutate(Source = 'GBD')
          )
     } else {
          data_region_jp <- rbind(
               model_incidence_who_region$data_export |> mutate(Source = 'WHO'),
               model_incidence_gbd_region$data_export |> mutate(Source = 'GBD')
          ) |>
               filter(location == region_names[i]) |>
               select(-location)
     }
     
     ggplot(data_region, aes(x = Year, y = Incidence, color = Source)) +
          geom_point() +
          geom_linerange(aes(ymin = lower, ymax = upper)) +
          geom_line(data = data_region_jp,
                    aes(x = year, y = model, color = Source)) +
          scale_y_continuous(breaks = breaks, limits = range(breaks), expand = c(0,0)) +
          scale_x_continuous(breaks = seq(2000, 2025, by = 5), limits = c(2000, 2025),
                             expand = expansion(add = c(1, 0))) +
          scale_color_manual(values = fill_color_source) +
          scale_fill_manual(values = fill_color_source) +
          labs(title = paste(LETTERS[i], region_names[i], sep = ': '),
               x = NULL,
               y = "Incidence rate per 100,000 population") +
          theme_bw()+
          theme(panel.grid = element_blank(),
                legend.position = ifelse(i == 1, 'inside', 'none'),
                legend.background = element_rect(fill = 'transparent'),
                legend.justification.inside = c(1, 1),
                legend.position.inside = c(0.99, 0.99))+
          guides(color = guide_legend(title = "Data Source", nrow = 1),
                 fill = guide_legend(title = "Data Source", nrow = 1))
}

# create panels
fig1 <- lapply(1:length(region_names), panel_line_function)

ggsave('./Output/Figure1_1.pdf',
       wrap_plots(fig1, ncol = 3, axis_titles = 'collect'),
       width = 12,
       height = 8,
       device = cairo_pdf,
       family = 'Helvetica')

# panel heatmap -----------------------------------------------------------

heatmap_breaks <- pretty(c(df_aapc_gbd$aapc, df_aapc_who$aapc))

panel_heatmap_function <- function(i){
     data_aapc <- list(df_aapc_who, df_aapc_gbd)[[i]] |> 
          # replace ~20 n Year
          mutate(Year = str_replace(Year, '~20', '~'))
     
     ggplot(data_aapc, aes(x = Year, y = location, fill = aapc)) +
          geom_tile(color = "white") +
          geom_text(aes(label = `AAPC (95%CI)`), color = "black", size = 3) +
          scale_fill_gradientn(colors = fill_color_heatmap,
                               breaks = heatmap_breaks,
                               limits = range(heatmap_breaks),
                               name = "AAPC") +
          scale_y_discrete(limits = rev(levels(data_aapc$location)),
                           expand = c(0,0)) +
          scale_x_discrete(expand = c(0,0)) +
          labs(title = paste(LETTERS[i + length(region_names)], 
                             c("WHO", "GBD")[i], 
                             sep = ': '),
               x = "Period",
               y = NULL) +
          theme_bw()+
          theme(legend.position = 'right',
                plot.title.position = "plot",
                panel.grid = element_blank())+
          guides(fill = guide_colorbar(barwidth = 0.7, barheight = 8, title.position = "top"))
}

# create panels
fig2 <- lapply(1:2, panel_heatmap_function)

ggsave('./Output/Figure1_2.pdf',
       wrap_plots(fig2, ncol = 2, axes = 'collect_y', guides = 'collect')&
            theme(legend.position = 'right'),
       width = 8,
       height = 8/3,
       device = cairo_pdf,
       family = 'Helvetica')

fig <- wrap_plots(c(fig1, fig2),
                  nrow = 3)

# save figure
ggsave('./Output/Figure 1.png',
       fig,
       width = 12,
       height = 8,
       dpi = 300)
