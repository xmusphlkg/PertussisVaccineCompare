#####################################
## @Description: Figure 3: Global maternal pertussis vaccination program
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-21 18:48:11
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-07 11:54:15
#####################################

# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(openxlsx)
library(paletteer)
library(Cairo)
library(sf)

rm(list = ls())

# data ----------------------------------------------------------------

## Map data ---------------------------------------------------

# Country iso codes
DataCountry <- read.csv('./Data/ISO_CODE.csv') |> 
     select(-location_name)

# Load map data
DataMap <- st_read('./Data/Map GS(2021)648 - geojson/globalmap.shp',
                   quiet = TRUE)

DataMapBorder <- st_read('./Data/Map GS(2021)648 - geojson/china_border.shp',
                         quiet = TRUE)

DataVaccine <- read.xlsx("./Data/VaccineSchedule.xlsx") 

# add Greenland data: https://pmc.ncbi.nlm.nih.gov/articles/PMC7034463/
DataVaccine <- DataVaccine |> 
     select(CODE, VaccinePregnant, VaccinePregnantTime) |> 
     rbind(data.frame(CODE = 'GRL',
                      VaccinePregnant = 1,
                      VaccinePregnantTime = NA))

DataVaccine <- DataVaccine |> 
     mutate(VaccinePregnant = factor(VaccinePregnant, 
                                     levels = c(1 ,0),
                                     labels = c('Recommended', 'Not recommended')))|>
     left_join(DataCountry, by = c('CODE' = 'ISO3'))


DataMapPlot <- DataMap |> 
     left_join(DataVaccine, by = c('SOC' = 'CODE'))

DataCoverageRaw <- read.xlsx("./Data/maternal immunization coverage.xlsx")

DataCoverage <- DataCoverageRaw |> 
     select(CODE, NAME, YEAR, COVERAGE_CATEGORY_DESCRIPTION, Coverage, Drop) |> 
     filter(is.na(Drop), !is.na(Coverage)) |> 
     select(-Drop) |> 
     # Clean Coverage values: extract numeric value from strings like
     # "Article 63.51 [46.95, 78.87]" or numeric entries
     mutate(Coverage = as.character(Coverage),
            Coverage = readr::parse_number(Coverage)) |> 
     # join country names / iso info for plotting
     left_join(DataCountry, by = c('CODE' = 'ISO3'))

# Extract latest Coverage for each country: keep rows with max YEAR,
# and if multiple entries in same year, take the mean Coverage.
DataCoverageLatest <- DataCoverage |> 
     group_by(CODE) |> 
     filter(YEAR == max(YEAR, na.rm = TRUE)) |> 
     summarize(
         YEAR = max(YEAR, na.rm = TRUE),
         Coverage = if (all(is.na(Coverage))) NA_real_ else mean(Coverage, na.rm = TRUE),
         COVERAGE_CATEGORY_DESCRIPTION = paste(unique(na.omit(COVERAGE_CATEGORY_DESCRIPTION)), collapse = "; "),
         NAME = dplyr::first(NAME)
     ) |> 
     left_join(DataCountry, by = c('CODE' = 'ISO3')) |> 
     filter(!is.na(ISO2))

# visualization: booster vaccination program ------------------------------------

## panel A ---------------------------------------

fill_color <- c("#2A6EBB", "#F0AB00")

p1 <- ggplot(DataMapPlot) +
     geom_sf(data = DataMapBorder, color = 'grey', fill = NA) +
     geom_sf(aes(fill = VaccinePregnant)) +
     coord_sf(xlim = c(-180, 180), ylim = c(-60, 75), expand = FALSE) +
     scale_fill_manual(values = fill_color,
                       na.translate = FALSE,
                       na.value = 'white') +
     labs(title = 'A',
          fill = NULL) +
     ggthemes::theme_map()+
     theme(legend.position = 'inside',
           plot.background = element_rect(fill = 'white', color = NA),
           legend.position.inside = c(0.05, 0.1),
           legend.justification.inside = c(0, 0),
           legend.key.size = unit(0.4, "cm"),
           legend.spacing.y = unit(0.1, "cm"),
           plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
           panel.border = element_rect(fill = NA, color = 'black'))

## panel B------------------------------------

DataTime <- DataVaccine |> 
     select(CODE, VaccinePregnantTime) |> 
     filter(!is.na(VaccinePregnantTime),
            VaccinePregnantTime != 'invisible',
            VaccinePregnantTime != 'No timing specified') |> 
     # split into two columns
     separate(VaccinePregnantTime, into = c('VaccinePregnantTime1', 'VaccinePregnantTime2'), sep = '-', remove = F) |> 
     # remove w in second column
     mutate(VaccinePregnantTime2 = gsub('w', '', VaccinePregnantTime2),
            VaccinePregnantTime2 = as.numeric(VaccinePregnantTime2),
            VaccinePregnantTime1 = as.numeric(VaccinePregnantTime1)) |> 
     left_join(DataCountry, by = c('CODE' = 'ISO3')) |> 
     arrange(VaccinePregnantTime1, VaccinePregnantTime2)

# periods of pregnancy
DataPeriods <- data.frame(
     Periods = c('1st Trimester', '2nd Trimester', '3rd Trimester'),
     Start = c(0, 14, 28),
     End = c(14, 28, 36)
) |> 
     mutate(Periods = factor(Periods, 
                             levels = c('1st Trimester', '2nd Trimester', '3rd Trimester'),
                             labels = c('1st Trimester', '2nd Trimester', '3rd Trimester')))

# Visualization
fill_color <- c("#C50084", "#7D5CC6", "#E37222")

p2 <- ggplot(DataTime) +
    geom_rect(data = DataPeriods,
              aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Periods),
              alpha = 0.5) +
    geom_linerange(aes(y = CODE, xmin = VaccinePregnantTime1, xmax = VaccinePregnantTime2),
                   linewidth = 2, color = "#2A6EBB") +
    coord_cartesian(xlim = c(8, 36)) +
    scale_y_discrete(limits = rev(DataTime$CODE),
                     labels = rev(DataTime$location_name_1)) +
    scale_x_continuous(breaks = seq(0, 36, 7),
                       labels = seq(0, 36, 7),
                       expand = c(0, 0)) +
    scale_fill_manual(values = fill_color,
                      name = NULL,
                      labels = c('1st Trimester\n(0-13 weeks)', 
                                 '2nd Trimester\n(14-27 weeks)', 
                                 '3rd Trimester\n(28-40 weeks)')) +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.justification.bottom = 'right',
          plot.title = element_text(hjust = -0.1),
          plot.title.position = 'panel',
          legend.text = element_text(margin = margin(l = 5, r = 20))) +
    labs(title = 'B',
         x = 'Gestational week of vaccination',
         y = NULL)

## panel C------------------------------------

plot_map_density <- function(x, fill_color) {
     # remove NA
     x <- x[!is.na(x)]
     y <- density(x, n = 2^12)
     data <- data.frame(x = y$x, y = y$y)
     
     ggplot(data = data, aes(x, y)) +
          geom_line() + 
          geom_segment(aes(xend = x, yend = 0, colour = x))+
          geom_hline(yintercept = 0, linetype = 'dashed')+
          scale_x_continuous(limits = c(0, 1),
                             breaks = seq(0, 1, by = 0.25),
                             labels = scales::percent_format(accuracy = 1),
                             expand = c(0, 0))+
          scale_y_continuous(limits = c(0, NA),
                             expand = expansion(mult = c(0, 0.2)))+
          scale_color_gradientn(colors = fill_color,
                                values = scales::rescale(c(0, 0.4, 0.6, 0.8, 0.9, 0.95, 0.97, 1)),
                                limits = c(0, 1),
                                labels = scales::percent_format(accuracy = 1))+
          theme_bw()+
          theme(panel.grid = element_blank(),
                axis.text.x = element_text(color = 'black', face = 'plain'),
                axis.title = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                panel.border = element_rect(color = 'black', fill = NA),
                plot.title = element_text(size = 11),
                plot.title.position = 'plot',
                legend.position = 'none')+
          labs(colour = NULL)+
          guides(colour = guide_colorbar(barwidth = 1))
}

fill_color_1 <- paletteer_d("MapPalettes::bruiser", direction = -1)
DataDTP <- DataCoverageLatest |> 
     rename(CoverageDTP = Coverage)

DataMapPlot1 <- DataMap |> 
     left_join(DataDTP, by = c('SOC' = 'CODE'))

p3_m <- plot_map_density(DataDTP$CoverageDTP/100, fill_color_1)+
     labs(title = 'Vaccine coverage')

p3 <- ggplot(data = DataMapPlot1) +
     geom_sf(data = DataMapBorder, color = 'grey', fill = NA) +
     geom_sf(aes(fill = CoverageDTP/100)) +
     # add x, y tick labels
     theme(axis.text.x = element_text(size = 8),
           axis.text.y = element_text(size = 8)) +
     coord_sf(xlim = c(-180, 180), ylim = c(-60, 75), expand = FALSE) +
     scale_fill_gradientn(colors = fill_color_1,
                          values = scales::rescale(c(0, 0.4, 0.6, 0.8, 0.9, 0.95, 0.97, 1)),
                          limits = c(0, 1),
                          breaks = c(0, 0.4, 0.6, 0.8, 0.9, 0.95, 0.97, 1),
                          labels = scales::percent_format(accuracy = 1),
                          na.value = "white")+
     ggthemes::theme_map()+
     theme(plot.background = element_rect(fill = 'white', color = NA),
           legend.position = 'none',
           plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
           panel.border = element_rect(fill = NA, color = 'black')) +
     labs(title = 'C', x = NULL, y = NULL)

p3 <- p3 + inset_element(p3_m, left = 0.01, bottom = 0.01, right = 0.25, top = 0.55)

# save figure data --------------------------------------------------------

design <- "
AB
CB
"

fig_1 <- p1 + free(p2) + p3 +
     plot_layout(design = design, widths = c(1.5, 1), heights = c(1, 1))

# save
ggsave('./Output/Figure 3.png',
       plot = fig_1,
       width = 12, height = 6, dpi = 300)

ggsave('./Output/Figure 3.pdf',
       plot = fig_1,
       width = 12, height = 6,
       device = cairo_pdf,
       family = "Helvetica")

write.xlsx(list('A' = DataVaccine,
                'B' = DataTime,
                'C' = DataCoverageLatest),
           './Output/Figure 3.xlsx')

# ---- Table S1: policy/timing + latest maternal coverage -----------------

to_percent <- function(x) {
     dplyr::if_else(is.na(x), NA_real_, dplyr::if_else(x <= 1, x * 100, x))
}

TableS1 <- DataVaccine |>
     select(CODE, location_name_1, VaccinePregnant, VaccinePregnantTime) |>
     full_join(
          DataCoverageLatest |>
               select(CODE, YEAR, Coverage, COVERAGE_CATEGORY_DESCRIPTION),
          by = "CODE"
     ) |>
     transmute(
          ISO3 = CODE,
          Country = location_name_1,
          MaternalPertussisRecommendation = VaccinePregnant,
          RecommendedTimingWindow = VaccinePregnantTime,
          MaternalCoverageYear = YEAR,
          MaternalCoveragePercent = to_percent(Coverage),
          MaternalCoverageType = COVERAGE_CATEGORY_DESCRIPTION
     ) |>
     arrange(Country) |> 
     filter(!(MaternalPertussisRecommendation == 'Not recommended' & is.na(MaternalCoveragePercent))) |> 
     filter(!is.na(Country))

write.csv(TableS1, "./Output/Table S1.csv", row.names = FALSE)
