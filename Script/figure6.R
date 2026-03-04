#####################################
## @Description: Figure 6 - Exploratory analysis of DTP1/DTP3 coverage vs pertussis incidence
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2026-03-02
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-02 18:58:00
#####################################

# packages ----------------------------------------------------------------

library(tidyverse)
library(openxlsx)
library(patchwork)
library(Cairo)

rm(list = ls())

# data --------------------------------------------------------------------

## Load incidence and vaccine program data (from Figure 4 pipeline)
load("./Output/DataAll.RData")

# DataAll contains WHO/GBD incidence (2019, 2021, 2024) and identifiers

## Load DTP1/DTP3 coverage data (same source as Figures 5 and 7) ---------

DataCoverage <- read.xlsx("./Data/Diphtheria tetanus toxoid and pertussis (DTP) vaccination coverage 2026-25-02 14-33 UTC.xlsx")

DataDTP <- DataCoverage |>
     filter(ANTIGEN %in% c("DTPCV1", "DTPCV3"),
            COVERAGE_CATEGORY == "WUENIC",
            GROUP == "COUNTRIES") |>
     transmute(
          Location_ID = CODE,
          Year        = as.integer(YEAR),
          Antigen     = ANTIGEN,
          Coverage    = as.numeric(COVERAGE) # percentage (0–100)
     ) |>
     distinct(Location_ID, Year, Antigen, .keep_all = TRUE) |>
     mutate(Antigen = recode(Antigen,
                             "DTPCV1" = "DTP1",
                             "DTPCV3" = "DTP3")) |>
     pivot_wider(names_from = Antigen, values_from = Coverage)

# long-format incidence by source/year -----------------------------------

who_long <- DataAll |>
     select(Location_ID, location_name = Location, starts_with("WHO_Inci_")) |>
     pivot_longer(
          cols      = matches("WHO_Inci_\\d{4}"),
          names_to  = "InciVar",
          values_to = "Incidence"
     ) |>
     mutate(
          Year   = as.integer(stringr::str_extract(InciVar, "\\d{4}")),
          source = "WHO"
     ) |>
     select(-InciVar)

gbd_long <- DataAll |>
     select(Location_ID, location_name = Location, starts_with("GBD_Inci_")) |>
     pivot_longer(
          cols      = matches("GBD_Inci_\\d{4}"),
          names_to  = "InciVar",
          values_to = "Incidence"
     ) |>
     mutate(
          Year   = as.integer(stringr::str_extract(InciVar, "\\d{4}")),
          source = "GBD"
     ) |>
     select(-InciVar)

DataLong <- bind_rows(who_long, gbd_long) |>
     filter(Year %in% c(2019, 2021, 2024)) |>
     left_join(DataDTP, by = c("Location_ID", "Year")) |>
     filter(!is.na(Incidence), !is.na(DTP1), !is.na(DTP3))

# sensitivity dataset: GBD all-age (crude) incidence (requires Figure 4 to have
# created GBD_AllAge_Inci_* columns in DataAll.RData)
gbd_allage_long <- DataAll |>
     select(Location_ID, location_name = Location, starts_with("GBD_AllAge_Inci_")) |>
     pivot_longer(
          cols      = matches("GBD_AllAge_Inci_\\d{4}"),
          names_to  = "InciVar",
          values_to = "Incidence"
     ) |>
     mutate(
          Year   = as.integer(stringr::str_extract(InciVar, "\\d{4}")),
          source = "GBD all-age"
     ) |>
     select(-InciVar) |>
     filter(Year %in% c(2019, 2021, 2024)) |>
     left_join(DataDTP, by = c("Location_ID", "Year")) |>
     filter(!is.na(Incidence), !is.na(DTP1), !is.na(DTP3))

# harmonise colour/legend with Figure 4 ---------------------------------

labs_vars_y <- c(
     "Incidence rate, 2019 (WHO)",
     "Incidence rate, 2021 (WHO)",
     "Incidence rate, 2024 (WHO)",
     "Incidence rate, 2019 (GBD)",
     "Incidence rate, 2021 (GBD)"
)

fill_color <- c("#2A6EBB", "#F0AB00", "#C50084", "#7D5CC6", "#E37222")
names(fill_color) <- labs_vars_y

DataLong <- DataLong |>
     mutate(
          source_label = dplyr::case_when(
               source == "WHO" & Year == 2019 ~ labs_vars_y[1],
               source == "WHO" & Year == 2021 ~ labs_vars_y[2],
               source == "WHO" & Year == 2024 ~ labs_vars_y[3],
               source == "GBD" & Year == 2019 ~ labs_vars_y[4],
               source == "GBD" & Year == 2021 ~ labs_vars_y[5],
               TRUE ~ NA_character_
          ),
          source_label = factor(source_label, levels = labs_vars_y)
     ) |>
     filter(!is.na(source_label))

# create DTP coverage categories for grouped comparisons -----------------

DataLong <- DataLong |>
     mutate(
          DTP1_cat = cut(DTP1,
                         breaks = c(-Inf, 60, 80, 90, Inf),
                         labels = c("<60%", "60–79%", "80–89%", "≥90%"),
                         right = FALSE),
          DTP3_cat = cut(DTP3,
                         breaks = c(-Inf, 60, 80, 90, Inf),
                         labels = c("<60%", "60–79%", "80–89%", "≥90%"),
                         right = FALSE)
     )

# simple correlation summary (Spearman) ----------------------------------

cor_summary <- DataLong |>
     group_by(source, Year) |>
     summarise(
          n        = n(),
          cor_DTP1 = suppressWarnings(cor(Incidence, DTP1, method = "spearman", use = "complete.obs")),
          cor_DTP3 = suppressWarnings(cor(Incidence, DTP3, method = "spearman", use = "complete.obs")),
          .groups  = "drop"
     )

print(cor_summary)

# resampling-based grouped comparisons (similar to Figure 4) -------------

bootstrap_n   <- 1000
permutation_n <- 5000

perm_test_pair <- function(vec, grp, n_perm = 5000, seed = NULL) {
     if (!is.null(seed)) set.seed(seed)
     levs <- levels(factor(grp))
     if (length(levs) != 2) return(NA_real_)
     obs <- mean(vec[grp == levs[2]], na.rm = TRUE) - mean(vec[grp == levs[1]], na.rm = TRUE)
     combined <- vec
     n <- length(vec)
     perm_diffs <- numeric(n_perm)
     for (b in seq_len(n_perm)) {
          perm_grp <- sample(grp, size = n, replace = FALSE)
          perm_diffs[b] <- mean(combined[perm_grp == levs[2]], na.rm = TRUE) -
               mean(combined[perm_grp == levs[1]], na.rm = TRUE)
     }
     mean(abs(perm_diffs) >= abs(obs), na.rm = TRUE)
}

compute_dtp_panel <- function(df, coverage_var, year, source = "WHO") {
     var_sym <- rlang::sym(coverage_var)
     dsub <- df |>
          filter(Year == year, source == !!source) |>
          select(Incidence, group = !!var_sym) |>
          drop_na()
     
     if (nrow(dsub) == 0) {
          return(list(bootstrap_summary = NULL, p_df = NULL))
     }
     
     dsub$group <- droplevels(dsub$group)
     unique_groups <- levels(dsub$group)
     
     # bootstrap mean and 95% CI per coverage category
     bootstrap_summary <- dsub |>
          group_by(group) |>
          summarise(n = n(), mean = mean(Incidence, na.rm = TRUE), .groups = "drop") |>
          rowwise() |>
          mutate(
               ci = list({
                    vals <- dsub$Incidence[dsub$group == group]
                    if (length(vals) <= 1) {
                         c(lower = mean, upper = mean)
                    } else {
                         boots <- replicate(bootstrap_n, mean(sample(vals, replace = TRUE), na.rm = TRUE))
                         stats::quantile(boots, probs = c(0.025, 0.975), na.rm = TRUE)
                    }
               })
          ) |>
          mutate(ci_lower = ci[[1]], ci_upper = ci[[2]]) |>
          select(-ci, -n) |>
          ungroup()
     
     # permutation-based pairwise tests between coverage categories
     if (length(unique_groups) > 1) {
          my_comparisons <- combn(unique_groups, 2, simplify = FALSE)
          p_list <- lapply(my_comparisons, function(pair) {
               g1 <- pair[1]; g2 <- pair[2]
               sub <- dsub |>
                    filter(group %in% c(g1, g2))
               if (nrow(sub) < 2) return(tibble(group1 = g1, group2 = g2, p = NA_real_))
               pval <- tryCatch(perm_test_pair(sub$Incidence, sub$group, n_perm = permutation_n),
                                error = function(e) NA_real_)
               tibble(group1 = g1, group2 = g2, p = pval)
          })
          p_df <- bind_rows(p_list) |>
               mutate(p.adj = p.adjust(p, method = "holm")) |>
               mutate(p.label = ifelse(is.na(p.adj), NA_character_,
                                       formatC(p.adj, format = "f", digits = 2)),
                      p.signif = dplyr::case_when(
                           is.na(p.adj) ~ "",
                           p.adj < 0.001 ~ "***",
                           p.adj < 0.01  ~ "**",
                           p.adj < 0.05  ~ "*",
                           TRUE          ~ "ns"
                      )) |>
               mutate(coverage_var = coverage_var,
                      Year = year,
                      source = source)
     } else {
          p_df <- NULL
     }
     
     bootstrap_summary <- bootstrap_summary |>
          mutate(coverage_var = coverage_var,
                 Year = year,
                 source = source)
     
     list(bootstrap_summary = bootstrap_summary, p_df = p_df)
}

years_use   <- sort(unique(DataLong$Year))
sources_use <- sort(unique(DataLong$source))
covars_use  <- c("DTP1_cat", "DTP3_cat")

dtp_boot_list <- list()
dtp_p_list    <- list()

for (cv in covars_use) {
     for (yr in years_use) {
          for (src in sources_use) {
               res <- compute_dtp_panel(DataLong, cv, yr, source = src)
               if (!is.null(res$bootstrap_summary)) {
                    dtp_boot_list[[length(dtp_boot_list) + 1]] <- res$bootstrap_summary
               }
               if (!is.null(res$p_df)) {
                    dtp_p_list[[length(dtp_p_list) + 1]] <- res$p_df
               }
          }
     }
}

dtp_boot_summary <- if (length(dtp_boot_list) > 0) bind_rows(dtp_boot_list) else NULL
dtp_p_values     <- if (length(dtp_p_list) > 0) bind_rows(dtp_p_list) else NULL

write.xlsx(list(dtp_boot_summary, dtp_p_values),
           file = "./Output/Figure 6.xlsx",
           overwrite = TRUE)

# visualization -----------------------------------------------------------

make_scatter_panels <- function(cov_var, cov_label, letter_offset = 0) {
     lapply(seq_along(years_use), function(i) {
          yr <- years_use[i]
          d <- DataLong |>
               filter(Year == yr)
          
          ggplot(d, aes(x = .data[[cov_var]], y = Incidence, colour = source_label)) +
               geom_point(alpha = 0.6, size = 1.2, show.legend = F) +
               geom_smooth(method = "loess", se = FALSE, show.legend = F) +
               scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
               scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, by = 20)) +
               labs(
                    x = paste0(cov_label, " coverage (%)"),
                    y = "Incidence rate per 100,000 population",
                    title = LETTERS[letter_offset + i],
                    colour = 'Data source'
               ) +
               scale_color_manual(values = fill_color, breaks = labs_vars_y, drop = FALSE) +
               theme_bw()
     })
}

p_dtp1_list <- make_scatter_panels("DTP1", "DTP1", letter_offset = 0)   # A, B, C
p_dtp3_list <- make_scatter_panels("DTP3", "DTP3", letter_offset = 3)   # D, E, F

## 2) Grouped comparisons (coverage categories, aligned with Figure 4) ---

dtp_boot_long <- dtp_boot_summary |>
     mutate(
          type_label = dplyr::case_when(
               source == "WHO" & Year == 2019 ~ labs_vars_y[1],
               source == "WHO" & Year == 2021 ~ labs_vars_y[2],
               source == "WHO" & Year == 2024 ~ labs_vars_y[3],
               source == "GBD" & Year == 2019 ~ labs_vars_y[4],
               source == "GBD" & Year == 2021 ~ labs_vars_y[5],
               TRUE ~ NA_character_
          ),
          type_label = factor(type_label, levels = labs_vars_y)
     ) |>
     filter(!is.na(type_label)) |>
     mutate(
          coverage_label = dplyr::case_when(
               coverage_var == "DTP1_cat" ~ "DTP1 coverage category",
               coverage_var == "DTP3_cat" ~ "DTP3 coverage category",
               TRUE ~ coverage_var
          )
     )

plot_grouped <- function(cov_var, src, letter) {
     d <- dtp_boot_long |>
          filter(coverage_var == cov_var, source == src)
     if (nrow(d) == 0) return(NULL)
     
     d |>
          ggplot(aes(x = group, y = mean, ymin = ci_lower, ymax = ci_upper,
                     colour = type_label)) +
          geom_pointrange(position = position_dodge(width = 0.6), show.legend = T) +
          geom_line(aes(group = type_label), position = position_dodge(width = 0.6), show.legend = F) +
          scale_color_manual(values = fill_color, breaks = labs_vars_y, drop = FALSE) +
          scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
          labs(
               x = unique(d$coverage_label),
               y = "Incidence rate per 100,000 population",
               title = letter,
               colour = 'Data source'
          ) +
          theme_bw()+
          guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
          theme(legend.position = "bottom",
                legend.title.position = 'top')
}

# DTP1 (G: WHO, H: GBD); DTP3 (I: WHO, J: GBD)
p_grp_dtp1_who <- plot_grouped("DTP1_cat", "WHO", LETTERS[7])
p_grp_dtp1_gbd <- plot_grouped("DTP1_cat", "GBD", LETTERS[8])
p_grp_dtp3_who <- plot_grouped("DTP3_cat", "WHO", LETTERS[9])
p_grp_dtp3_gbd <- plot_grouped("DTP3_cat", "GBD", LETTERS[10])

left_panels <- wrap_plots(c(p_dtp1_list, p_dtp3_list),
                          ncol = 2, guides = "collect", byrow = F,
                          axes = 'collect',
                          axis_titles = 'collect') &
     theme(legend.position = "bottom")

right_panels <- wrap_plots(c(p_grp_dtp1_who, p_grp_dtp1_gbd, p_grp_dtp3_who, p_grp_dtp3_gbd),
                           ncol = 2, guides = "collect", byrow = T,
                           axis_titles = 'collect') &
     theme(legend.position = "bottom")

fig6 <- cowplot::plot_grid(left_panels, right_panels, ncol = 2, rel_widths = c(1, 1.1))

ggsave("./Output/Figure 6.png",
       plot = fig6,
       width = 12, height = 8, dpi = 300)

ggsave("./Output/Figure 6.pdf",
       plot = fig6,
       width = 12, height = 12,
       device = cairo_pdf,
       family = "Helvetica")

# Sensitivity figure (Supplement Figure S2): GBD all-age (crude) only
# Keep it lightweight: DTP3 vs incidence, years 2019 & 2021.
DataLong_allage_s2 <- gbd_allage_long |>
     filter(Year %in% c(2019, 2021)) |>
     mutate(year = factor(Year, levels = c(2019, 2021)))

fill_year_s2 <- c("2019" = "#7D5CC6", "2021" = "#E37222")

p_s2 <- ggplot(DataLong_allage_s2, aes(x = DTP3, y = Incidence, colour = year)) +
     geom_point(alpha = 0.6, size = 1.2, show.legend = TRUE) +
     geom_smooth(method = "loess", se = FALSE, show.legend = FALSE) +
     scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
     scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, by = 20)) +
     scale_color_manual(values = fill_year_s2) +
     labs(
          x = "DTP3 coverage (%)",
          y = "Incidence rate per 100,000 population",
          colour = "Year"
     ) +
     theme_bw() +
     theme(panel.grid = element_blank(),
           legend.position = "top")

ggsave("./Output/Supplement_Figure_S2.png",
       plot = p_s2,
       width = 6, height = 4, dpi = 300)
