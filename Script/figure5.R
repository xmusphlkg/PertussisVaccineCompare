#####################################
## @Description: Figure 5: Global vaccine coverage of DTP1 and DTP3 in 2019, 2021, and 2024
## @version: 
## @Author: Li Kangguo
## @Date: 2025-10-21 19:20:50
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-02 11:44:42
#####################################

# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(openxlsx)
library(paletteer)
library(Cairo)
library(sf)

rm(list = ls())


# function ----------------------------------------------------------------

plot_map_density <- function(x, fill_color) {
     # remove NA
     x <- x[!is.na(x)]
     y <- density(x, n = 2^12)
     data <- data.frame(x = y$x, y = y$y)
     
     ggplot(data = data, aes(x, y)) +
          geom_line() + 
          geom_segment(aes(xend = x, yend = 0, colour = x))+
          geom_hline(yintercept = 0, linetype = 'dashed')+
          scale_x_continuous(limits = c(0.2, 1),
                             breaks = seq(0.2, 1, by = 0.2),
                             labels = scales::percent_format(accuracy = 1),
                             expand = c(0, 0))+
          scale_y_continuous(limits = c(0, NA),
                             expand = expansion(mult = c(0, 0.2)))+
          scale_color_gradientn(colors = fill_color,
                                values = scales::rescale(c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.97, 1)),
                                limits = c(0.2, 1),
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

# Data --------------------------------------------------------------------

DataCoverage <- read.xlsx('./Data/Diphtheria tetanus toxoid and pertussis (DTP) vaccination coverage 2026-25-02 14-33 UTC.xlsx')

DataDTP3 <- DataCoverage |>
     filter(ANTIGEN %in% c("DTPCV1", "DTPCV3"), 
            COVERAGE_CATEGORY == 'WUENIC', 
            GROUP == 'COUNTRIES') |> 
     transmute(Location_ID = CODE,
               Year = as.integer(YEAR),
               Antigen = ANTIGEN,
               Coverage = as.numeric(COVERAGE)/100) |>
     distinct(Location_ID, Year, Antigen, .keep_all = TRUE) |>
     mutate(Antigen = recode(Antigen,
                             "DTPCV1" = "DTP1",
                             "DTPCV3" = "DTP3")) |>
     pivot_wider(names_from = Antigen, values_from = Coverage)

# Load map data
DataMap <- st_read('./Data/Map GS(2021)648 - geojson/globalmap.shp',
                   quiet = TRUE)

DataMapBorder <- st_read('./Data/Map GS(2021)648 - geojson/china_border.shp',
                         quiet = TRUE)

# find the country not in Map data
DataCoverage[!DataCoverage$Location_ID %in% DataMap$SOC, 'NAME']

# visual ------------------------------------------------------------------

plot_DTP <- function(i) {
     selection <- DataCoverage |> 
          select(ANTIGEN, YEAR) |> 
          filter(YEAR %in% c(2019, 2021, 2024)) |>
          unique() |> 
          arrange(ANTIGEN, YEAR) |> 
          mutate(ANTIGEN = recode(ANTIGEN,
                                  "DTPCV1" = "DTP1",
                                  "DTPCV3" = "DTP3"))
     antigen <- selection$ANTIGEN[i]
     year <- selection$YEAR[i]
     
     if (antigen == 'DTP1') {
          fill_color <- paletteer_d("MapPalettes::bruiser", direction = -1)
          DataDTP <- DataDTP3 |> 
               filter(Year == year) |> 
               select(Location_ID, CoverageDTP = DTP1)
     } else if (antigen == 'DTP3') {
          fill_color <- paletteer_d("MapPalettes::green_machine", direction = -1)
          
          DataDTP <- DataDTP3 |> 
               filter(Year == year) |> 
               select(Location_ID, CoverageDTP = DTP3)
     } else {
          stop('Invalid antigen. Please choose either "DTP1" or "DTP3".')
     }
     
     DataMapPlot <- DataMap |> 
          left_join(DataDTP, by = c('SOC' = 'Location_ID'))
     
     fig_1_m <- plot_map_density(DataMapPlot$CoverageDTP, fill_color)+
          labs(title = paste('Vaccine coverage,\n', antigen, ' in ', year, sep = ''))
     
     fig_1 <- ggplot(data = DataMapPlot) +
          geom_sf(data = DataMapBorder, color = 'grey', fill = NA) +
          geom_sf(aes(fill = CoverageDTP)) +
          # add x, y tick labels
          theme(axis.text.x = element_text(size = 8),
                axis.text.y = element_text(size = 8)) +
          scale_x_continuous(limits = c(-180, 180),
                             expand = c(0, 0)) + 
          scale_y_continuous(limits = c(-60, 75)) +
          scale_fill_gradientn(colors = fill_color,
                               values = scales::rescale(c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.97, 1)),
                               limits = c(0.2, 1),
                               breaks = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.97, 1),
                               labels = scales::percent_format(accuracy = 1),
                               na.value = "white")+
          ggthemes::theme_map()+
              theme(plot.background = element_rect(fill = 'white', color = NA),
                   legend.position = 'none',
                   plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
                   panel.border = element_rect(fill = NA, color = 'black')) +
          labs(title = LETTERS[i], x = NULL, y = NULL)
     
     fig_1 <- fig_1 + inset_element(fig_1_m, left = 0.01, bottom = 0.01, right = 0.25, top = 0.55)
}

fig_1 <- plot_DTP(1)

fig_2 <- plot_DTP(2)

fig_3 <- plot_DTP(3)

fig_4 <- plot_DTP(4)

fig_5 <- plot_DTP(5)

fig_6 <- plot_DTP(6)

fig <- cowplot::plot_grid(
     fig_1, fig_2, fig_3, fig_4, fig_5, fig_6,
     ncol = 2, align = 'v', rel_widths = c(1, 1),
     byrow = F
)

ggsave("./Output/Figure 5.pdf",
       fig,
       width = 12,
       height = 8.5,
       device = cairo_pdf)

ggsave("./Output/Figure 5.png",
       fig,
       width = 12,
       height = 8.5)

write.xlsx(DataDTP3 |> 
                filter(Year %in% c(2019, 2021, 2024)),
           "./Output/Figure 5.xlsx")