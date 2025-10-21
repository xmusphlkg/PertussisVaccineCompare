# packages ----------------------------------------------------------------

library(nih.joinpoint)
library(segmented)

library(tidyverse)
library(openxlsx)
library(stringr)
library(patchwork)
library(paletteer)

rm(list = ls())

source('./script/joinpoint_setting.R')

fill_color_source <- paletteer_d("werpals::pan")[1:2]
fill_color_heatmap <- paletteer_d("awtools::a_palette")

# data --------------------------------------------------------------------

## WHO --------------------------------------------------------------------
region_names <- c(
     "African Region", "Eastern Mediterranean Region", "European Region",
     "Region of the Americas", "South-East Asia Region", "Western Pacific Region",
     "Global", 'WHO'
)

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
model_incidence_who_gloabl <- joinpoint(df_incidence_who_region |> filter(Location == 'Global', Year %in% 2000:2019),
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
df_aapc_who <- rbind(get_aapc(model_incidence_who_gloabl) |> mutate(location = 'Global'),
                     get_aapc(model_incidence_who_region)) |>
     mutate(`AAPC (95%CI)` = paste0(aapc, p_value_label),
            aapc = as.numeric(aapc),
            location = factor(location, levels = region_names))

df_aapc_gbd <- rbind(get_aapc(model_incidence_gbd_global) |> mutate(location = 'Global'),
                     get_aapc(model_incidence_gbd_region)) |>
     mutate(`AAPC (95%CI)` = paste0(aapc, p_value_label),
            aapc = as.numeric(aapc),
            location = factor(location, levels = region_names))

# figure ------------------------------------------------------------------

## panel line --------------------------------------------------------------

i <- 1

panel_line_function <- function(i){
     data_region <- df_incidence_region |> 
          filter(Location == region_names[i])
     
     breaks <- pretty(c(0, data_region$Incidence, data_region$upper))
     
     ggplot(data_region, aes(x = Year, y = Incidence, color = Source)) +
          geom_line() +
          geom_point()+
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source), alpha = 0.2, color = NA) +
          scale_y_continuous(breaks = breaks, limits = range(breaks), expand = c(0,0)) +
          scale_x_continuous(breaks = seq(2000, 2025, by = 5), limits = c(2000, 2025), expand = c(0,0)) +
          scale_color_manual(values = fill_color_source) +
          scale_fill_manual(values = fill_color_source) +
          labs(title = paste(LETTERS[i], region_names[i], sep = ': '),
               x = NULL,
               y = "Incidence rate (per 100,000 population)") +
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

ggsave('./Output/Figure 1_1.pdf',
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

ggsave('./Output/Figure 1_2.pdf',
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

# save figure data --------------------------------------------------------

write.xlsx(list('A-G' = df_incidence_region,
                'H' = df_aapc_who,
                'I' = df_aapc_gbd),
           './Output/Figure 1.xlsx')

