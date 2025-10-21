
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

# visualization: booster vaccination program ------------------------------------

fill_color <- paletteer_d("werpals::pan")

p1 <- ggplot(DataMapPlot) +
     geom_sf(data = DataMapBorder, color = 'grey', fill = NA) +
     geom_sf(aes(fill = VaccinePregnant), color = 'white') +
     coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE) +
     scale_fill_manual(values = fill_color,
                       na.translate = FALSE,
                       na.value = 'grey') +
     labs(title = 'A',
          fill = NULL) +
     ggthemes::theme_map()+
     theme(legend.position = 'inside',
           legend.justification.inside = 'bottom',
           plot.background = element_rect(fill = 'white', color = NA),
           legend.position.inside = c(0.1, 0.1),
           plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
           panel.border = element_rect(fill = NA, color = 'black'))

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
fill_color <- paletteer_d("werpals::pan")

# split data into two parts for better visualization
n <- nrow(DataTime)
idx <- seq_len(n)
DataTime_A <- DataTime[idx <= n/2, ]
DataTime_B <- DataTime[idx > n/2, ]

p2_1 <- ggplot(DataTime_A) +
    geom_rect(data = DataPeriods,
              aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Periods),
              alpha = 0.5) +
    geom_linerange(aes(y = CODE, xmin = VaccinePregnantTime1, xmax = VaccinePregnantTime2),
                   linewidth = 2, color = "#5785C1FF") +
    coord_cartesian(xlim = c(8, 36)) +
    scale_y_discrete(limits = rev(DataTime_A$CODE),
                     labels = rev(DataTime_A$location_name_1)) +
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
          plot.margin = margin(t = 0, r = 10, b = 0, l = 40),
          plot.title = element_text(hjust = -0.1),
          plot.title.position = 'panel',
          legend.text = element_text(margin = margin(l = 5, r = 20))) +
    labs(title = 'B',
         x = 'Recommended time for maternal vaccination (weeks)',
         y = NULL)

p2_2 <- ggplot(DataTime_B) +
    geom_rect(data = DataPeriods,
              aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Periods),
              alpha = 0.5) +
    geom_linerange(aes(y = CODE, xmin = VaccinePregnantTime1, xmax = VaccinePregnantTime2),
                   linewidth = 2, color = "#5785C1FF") +
    coord_cartesian(xlim = c(8, 36)) +
    scale_y_discrete(limits = rev(DataTime_B$CODE),
                     labels = rev(DataTime_B$location_name_1),
                     position = "right") +
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
          plot.margin = margin(t = 0, r = 10, b = 0, l = 40),
          plot.title = element_text(hjust = -0.1),
          plot.title.position = 'panel',
          legend.text = element_text(margin = margin(l = 5, r = 20))) +
    labs(title = 'C',
         x = 'Recommended time for maternal vaccination (weeks)',
         y = NULL)

p2 <- p2_1 + p2_2 +
     plot_layout(guides = 'collect', axis_titles = 'collect') &
     theme(legend.position = 'bottom')

# save
ggsave('./Output/Figure 3.png',
       plot = cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.54, 1)),
       # plot = p1 / p2,
       width = 7, height = 9.3, dpi = 300)

ggsave('./Output/Figure 3.pdf',
       plot = cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.54, 1)),
       # plot = p1 / p2,
       width = 7, height = 9.3,
       device = cairo_pdf,
       family = "Helvetica")

# save figure data --------------------------------------------------------

write.xlsx(list('A' = DataVaccine,
                'B' = DataTime_A,
                'C' = DataTime_B),
           './Output/Figure 3.xlsx')
