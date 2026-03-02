#####################################
## @Description: Figure 2: Global routine pertussis vaccination program
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-21 18:40:46
## @LastEditors: Li Kangguo
## @LastEditTime: 2025-10-22 09:18:37
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

DataVaccine <- read.xlsx('./Data/VaccineSchedule.xlsx')

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

DataVaccine <- DataVaccine |> 
     mutate(TimeFirstShotG = cut(as.numeric(TimeFirstShot), 
                                 breaks = c(-Inf, 1.5, 2.5, Inf),
                                 labels = c('[1.0,1.5)m', 
                                            '[1.5,2.5)m',
                                            '[2.5,3.0]m'),
                                 include.lowest = TRUE),
            TimeLastShotG = cut(as.numeric(TimeLastShot), 
                                breaks = c(-Inf, 12, 30, 66, 150, Inf),
                                labels = c('[3,6]m', 
                                           '[1,2]y',
                                           '[3,5]y',
                                           '[6,12]y',
                                           '[13,18]y'),
                                include.lowest = TRUE),
            VaccineDose = as.factor(VaccineDose))


DataMapPlot <- DataMap |> 
     left_join(DataVaccine, by = c('SOC' = 'CODE'))

# visualization: routine vaccination program ------------------------------------

# Routine vaccine dose
Data <- DataVaccine |> 
     group_by(VaccineDose) |> 
     summarise(Count = n()) |> 
     ungroup()

breaks <- pretty(c(Data$Count, 0), n = 5)

fill_color <- c("#2A6EBB", "#F0AB00", "#C50084", "#7D5CC6")
names(fill_color) <- levels(DataVaccine$VaccineDose)

p1 <- ggplot(Data) +
     geom_bar(aes(x = VaccineDose, y = Count, fill = VaccineDose), 
              stat = 'identity', width = 0.75) +
     geom_text(aes(x = VaccineDose, y = Count, label = Count), 
               vjust = -0.3,
               size = 4) +
     scale_fill_manual(values = fill_color) +
     scale_y_continuous(breaks = breaks,
                        limits = range(breaks),
                        expand = expansion(mult = c(0, 0))) +
     theme_bw()+
     theme(plot.title.position = 'plot',
           legend.position = 'none')+
     labs(title = 'A: Number of scheduled doses',
          x = NULL,
          y = 'Number of countries')

p2 <- ggplot(DataMapPlot) +
     geom_sf(data = DataMapBorder, color = 'grey', fill = NA) +
     geom_sf(aes(fill = VaccineDose), color = 'white') +
     coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE) +
     scale_fill_manual(values = fill_color,
                       na.value = 'grey') +
     labs(title = NULL) +
     ggthemes::theme_map()+
     theme(legend.position = 'none',
           panel.border = element_rect(fill = NA, color = 'black'))

# First shot time
Data <- DataVaccine |> 
     group_by(TimeFirstShotG) |> 
     summarise(Count = n()) |> 
     ungroup()

fill_color <- c("#7D5CC6", "#E37222", "#69BE28")
names(fill_color) <- levels(DataVaccine$TimeFirstShotG)
breaks <- pretty(c(Data$Count, 0), n = 5)

p3 <- ggplot(Data) +
     geom_bar(aes(x = TimeFirstShotG, y = Count, fill = TimeFirstShotG), 
              stat = 'identity', width = 0.75) +
     geom_text(aes(x = TimeFirstShotG, y = Count, label = Count), 
               vjust = -0.3,
               size = 4) +
     scale_fill_manual(values = fill_color) +
     scale_y_continuous(breaks = breaks,
                        limits = range(breaks),
                        expand = expansion(mult = c(0, 0))) +
     theme_bw()+
     theme(plot.title.position = 'plot',
           legend.position = 'none')+
     labs(title = 'B: Timing of first scheduled dose',
          x = NULL,
          y = 'Number of countries')

p4 <- ggplot(DataMapPlot) +
     geom_sf(data = DataMapBorder, color = 'grey', fill = NA) +
     geom_sf(aes(fill = TimeFirstShotG), color = 'white') +
     coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE) +
     scale_fill_manual(values = fill_color,
                       na.value = 'grey') +
     labs(title = NULL) +
     ggthemes::theme_map()+
     theme(legend.position = 'none',
           panel.border = element_rect(fill = NA, color = 'black'))

# Last shot time
Data <- DataVaccine |> 
     group_by(TimeLastShotG) |> 
     summarise(Count = n()) |> 
     ungroup()

fill_color <- c("#E37222", "#69BE28", "#00B2A9", "#CD202C", "#747678")
names(fill_color) <- levels(DataVaccine$TimeLastShotG)
breaks <- pretty(c(Data$Count, 0), n = 5)

p5 <- ggplot(Data) +
     geom_bar(aes(x = TimeLastShotG, y = Count, fill = TimeLastShotG), 
              stat = 'identity', width = 0.75) +
     geom_text(aes(x = TimeLastShotG, y = Count, label = Count), 
               vjust = -0.3,
               size = 4) +
     scale_fill_manual(values = fill_color) +
     scale_y_continuous(breaks = breaks,
                        limits = range(breaks),
                        expand = expansion(mult = c(0, 0))) +
     theme_bw()+
     theme(plot.title.position = 'plot',
           legend.position = 'none')+
     labs(title = 'C: Timing of last scheduled childhood dose',
          x = NULL,
          y = 'Number of countries',
          fill = 'Timing of last childhood dose')

p6 <- ggplot(DataMapPlot) +
     geom_sf(data = DataMapBorder, color = 'grey', fill = NA) +
     geom_sf(aes(fill = TimeLastShotG), color = 'white') +
     coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE) +
     scale_fill_manual(values = fill_color,
                       na.value = 'grey') +
     labs(title = NULL) +
     ggthemes::theme_map()+
     theme(legend.position = 'none',
           panel.border = element_rect(fill = NA, color = 'black'))

# combine plots

fig <- p1 + p2 + p3 + p4 + p5 + p6 +
     plot_layout(ncol = 2, nrow = 3, byrow = TRUE, widths = c(0.5, 1))

# save
ggsave('./Output/Figure 2.png',
       plot = fig,
       width = 9, height = 9, dpi = 300)

ggsave('./Output/Figure 2.pdf',
       plot = fig,
       width = 9, height = 9,
       device = cairo_pdf,
       family = 'Helvetica')

# save figure data --------------------------------------------------------

write.xlsx(list('A' = DataVaccine |> select(CODE, NAME, VaccineDose),
                'B' = DataVaccine |> select(CODE, NAME, TimeFirstShotG),
                'C' = DataVaccine |> select(CODE, NAME, TimeLastShotG)),
           './Output/Figure 2.xlsx')
