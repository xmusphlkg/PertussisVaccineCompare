##############################################################
##  GLMM Forest Plot Pipeline (Negative Binomial / ZINB)
##  Author: Nicole Phillips
##  Purpose: Stabilize RR estimation via shrinkage & pooling
##  Date: 2025-10-22
##############################################################

# ======================
# 1. Load dependencies
# ======================
pkgs <- c("tidyverse","forcats","glmmTMB","DHARMa","performance","emmeans")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(pkgs, library, character.only = TRUE))

# ======================
# 2. Preprocess data
# ======================
prep_data <- function(df,
                      ref_dose = "4",
                      ref_first = "[1.5,2.5)m",
                      ref_last  = "[3,6]m",
                      ref_preg  = "Not recommended",
                      ref_src   = "WHO") {
     df %>%
          mutate(
               VaccineDose      = fct_relevel(as.factor(VaccineDose), ref_dose),
               TimeFirstShotG   = fct_relevel(as.factor(TimeFirstShotG), ref_first),
               TimeLastShotG    = fct_relevel(as.factor(TimeLastShotG),  ref_last),
               VaccinePregnant  = fct_relevel(as.factor(VaccinePregnant), ref_preg),
               source           = fct_relevel(as.factor(source), ref_src),
               Year             = as.integer(Year),
               location_id      = as.factor(location_id),
               Cases            = round(Cases)    # ensure integer counts
          )
}

# ======================
# 3. Common-country selector
# ======================
common_countries <- function(df, years, sources, mode = c("within_year","across_years")) {
     mode <- match.arg(mode)
     if (mode == "within_year") {
          keep_ids <- df %>%
               filter(Year %in% years, source %in% sources) %>%
               distinct(Year, source, location_id) %>%
               count(Year, location_id, name = "n_src") %>%
               filter(n_src == length(sources)) %>%
               select(Year, location_id)
          df %>% inner_join(keep_ids, by = c("Year","location_id"))
     } else {
          ids <- df %>%
               filter(Year %in% years, source %in% sources) %>%
               distinct(Year, source, location_id)
          ids_all <- ids %>% count(location_id) %>%
               filter(n == length(years) * length(sources)) %>% pull(location_id)
          df %>% filter(location_id %in% ids_all, Year %in% years, source %in% sources)
     }
}

# ======================
# 4. Model fitting: NB ↔ ZINB + shrinkage control
# ======================
fit_nb_or_zinb <- function(df, year, src, zi_test = TRUE,
                           shrinkage = c("random","fixed")) {
     shrinkage <- match.arg(shrinkage)
     dsub <- df %>% filter(Year == year, source == src) %>% droplevels()
     stopifnot(nrow(dsub) > 0)
     
     if (shrinkage == "random") {
          base_form <- as.formula(
               "Cases ~ 1 + VaccinePregnant + (1|VaccineDose) + (1|TimeFirstShotG) + (1|TimeLastShotG) + (1|location_id)"
          )
     } else {
          base_form <- as.formula(
               "Cases ~ 1 + VaccineDose + TimeFirstShotG + TimeLastShotG + VaccinePregnant + (1|location_id)"
          )
     }
     
     fit_nb <- glmmTMB(
          formula   = base_form,
          data      = dsub,
          offset    = log(Population),
          family    = nbinom2,
          ziformula = ~0,
          control   = glmmTMBControl(optCtrl = list(iter.max=1e4, eval.max=1e4))
     )
     
     best_fit  <- fit_nb
     best_name <- if (shrinkage=="random") "NB (shrink)" else "NB"
     
     if (isTRUE(zi_test)) {
          res <- tryCatch(DHARMa::simulateResiduals(best_fit), error=function(e) NULL)
          if (!is.null(res)) {
               zi_p <- tryCatch(DHARMa::testZeroInflation(res)$p.value, error=function(e) NA)
               if (!is.na(zi_p) && zi_p < 0.05) {
                    fit_zinb <- tryCatch(update(fit_nb, ziformula = ~1), error=function(e) NULL)
                    if (!is.null(fit_zinb)) {
                         comp <- tryCatch(AIC(fit_nb, fit_zinb), error=function(e) data.frame(AIC=c(NA,NA)))
                         if (!any(is.na(comp$AIC)) && comp$AIC[2] + 2 < comp$AIC[1]) {
                              best_fit  <- fit_zinb
                              best_name <- if (shrinkage=="random") "ZINB (shrink)" else "ZINB"
                         }
                    }
               }
          }
     }
     
     list(fit = best_fit, model = best_name, data = dsub, shrinkage = shrinkage)
}

# ======================
# 5. RR extraction helpers
# ======================

## 5.1 fixed-effect RR (emmeans)
emm_rr <- function(fit_obj, var, ref_level = NULL) {
     rg <- emmeans::emmeans(fit_obj, specs = stats::as.formula(paste0("~ ", var)))
     levs <- levels(rg@grid[[var]])
     if (length(levs) <= 1L) {
          return(tibble::tibble(var = character(), level = character(), RR = numeric(),
                                lower = numeric(), upper = numeric()))
     }
     ref_idx <- if (is.null(ref_level)) 1L else match(ref_level, levs)
     if (is.na(ref_idx)) ref_idx <- 1L
     contr <- emmeans::contrast(rg, method = "trt.vs.ctrl", ref = ref_idx)
     sm <- summary(contr, type = "response")
     lvl <- sub("\\s*(/|-)\\s*.*$", "", as.character(sm$contrast))
     tibble::tibble(
          var   = rep(var, length(lvl)),
          level = lvl,
          RR    = sm$ratio,
          lower = sm$lower.CL,
          upper = sm$upper.CL
     )
}

## 5.2 random-effect RR (via BLUPs)
rand_rr_glmmTMB <- function(fit_obj, group, ref_level) {
     re_all <- glmmTMB::ranef(fit_obj, condVar = TRUE)
     empty <- tibble::tibble(var = character(), level = character(), RR = numeric(),
                             lower = numeric(), upper = numeric())
     if (!(group %in% names(re_all))) return(empty)
     re_grp <- re_all[[group]]
     lev <- rownames(re_grp)
     if (length(lev) <= 1L) return(empty)
     ref_idx <- match(ref_level, lev)
     if (is.na(ref_idx)) ref_idx <- 1L
     u <- re_grp[["(Intercept)"]]
     cv <- attr(re_grp, "condVar")
     vars <- tryCatch({
          if (is.null(cv)) rep(NA_real_, length(lev))
          else if (!is.null(dim(cv)) && length(dim(cv)) == 3L) sapply(seq_along(lev), function(i) cv[1,1,i])
          else if (is.list(cv)) sapply(cv, function(m) m[1,1])
          else rep(NA_real_, length(lev))
     }, error = function(e) rep(NA_real_, length(lev)))
     idx <- setdiff(seq_along(lev), ref_idx)
     out <- lapply(idx, function(i) {
          logRR <- as.numeric(u[i] - u[ref_idx])
          se    <- if (is.na(vars[i]) || is.na(vars[ref_idx])) NA_real_ else sqrt(vars[i] + vars[ref_idx])
          RR <- exp(logRR)
          lo <- if (is.na(se)) NA_real_ else exp(logRR - 1.96 * se)
          hi <- if (is.na(se)) NA_real_ else exp(logRR + 1.96 * se)
          tibble::tibble(var = group, level = lev[i], RR = RR, lower = lo, upper = hi)
     })
     dplyr::bind_rows(out)
}

## 5.3 unified wrapper
extract_all_rr <- function(fit, ref_levels, shrinkage) {
     parts <- list()
     if (shrinkage == "fixed") {
          parts$Dose   = emm_rr(fit, "VaccineDose",     ref_levels[["VaccineDose"]])
          parts$First  = emm_rr(fit, "TimeFirstShotG",  ref_levels[["TimeFirstShotG"]])
          parts$Last   = emm_rr(fit, "TimeLastShotG",   ref_levels[["TimeLastShotG"]])
          parts$Preg   = emm_rr(fit, "VaccinePregnant", ref_levels[["VaccinePregnant"]])
     } else {
          parts$Dose   = rand_rr_glmmTMB(fit, "VaccineDose",    ref_levels[["VaccineDose"]])
          parts$First  = rand_rr_glmmTMB(fit, "TimeFirstShotG", ref_levels[["TimeFirstShotG"]])
          parts$Last   = rand_rr_glmmTMB(fit, "TimeLastShotG",  ref_levels[["TimeLastShotG"]])
          parts$Preg   = emm_rr(fit, "VaccinePregnant", ref_levels[["VaccinePregnant"]])
     }
     parts <- parts[vapply(parts, nrow, integer(1L)) > 0]
     if (length(parts) == 0) return(tibble::tibble(var = character(), level = character(),
                                                   RR = numeric(), lower = numeric(), upper = numeric()))
     dplyr::bind_rows(parts)
}

# ======================
# 6. Pipeline runner
# ======================
run_pipeline <- function(df,
                         years   = c(2019, 2021, 2024),
                         sources = c("WHO","GBD"),
                         mode_common = "within_year",
                         zi_test = TRUE,
                         shrinkage = c("random","fixed")) {
     
     shrinkage <- match.arg(shrinkage)
     df_keep <- common_countries(df, years, sources, mode = mode_common)
     
     ref_levels <- list(
          VaccineDose     = levels(df_keep$VaccineDose)[1],
          TimeFirstShotG  = levels(df_keep$TimeFirstShotG)[1],
          TimeLastShotG   = levels(df_keep$TimeLastShotG)[1],
          VaccinePregnant = levels(df_keep$VaccinePregnant)[1]
     )
     
     results <- list(); fits <- list()
     
     for (yr in years) for (src in sources) {
          dsub <- df_keep %>% filter(Year == yr, source == src)
          if (nrow(dsub) == 0) next
          fitobj <- fit_nb_or_zinb(df_keep, yr, src, zi_test = zi_test, shrinkage = shrinkage)
          fits[[paste(yr, src, sep = "_")]] <- fitobj
          rr <- extract_all_rr(fitobj$fit, ref_levels, shrinkage) %>%
               mutate(Year = yr, Source = src, Model = fitobj$model)
          results[[paste(yr, src, sep = "_")]] <- rr
     }
     
     list(rr = bind_rows(results), fits = fits, data_used = df_keep)
}

# ======================
# 7. Forest plot
# ======================
plot_forest <- function(rr_df) {
     rr_df %>%
          mutate(
               var = factor(var, levels = c("VaccineDose","TimeFirstShotG","TimeLastShotG","VaccinePregnant"),
                            labels = c("Dose","First shot","Last shot","Pregnancy")),
               panel = factor(Year, levels = c(2019,2021,2024),
                              labels = c("A: 2019","B: 2021","C: 2024"))
          ) %>%
          ggplot(aes(x = RR, y = interaction(var, level, lex.order = TRUE), color = Source)) +
          geom_point(position = position_dodge(width = 0.6)) +
          geom_errorbarh(aes(xmin = pmax(lower, .Machine$double.eps), xmax = upper),
                         height = 0.22, position = position_dodge(width = 0.6)) +
          geom_vline(xintercept = 1, linetype = 2) +
          scale_x_log10() +
          labs(x = "Rate ratio (log scale)", y = NULL, color = "Source") +
          facet_wrap(~ panel, nrow = 1, scales = "free_y") +
          theme_minimal(base_size = 12)
}

# ======================
# 8. Example usage
# ======================
# Example:
DataLongCounts2 <- prep_data(DataLongCounts)
out <- run_pipeline(DataLongCounts2,
                    years = c(2019,2021,2024),
                    sources = c("WHO","GBD"),
                    mode_common = "within_year",
                    zi_test = TRUE,
                    shrinkage = "random")
rr_table <- out$rr
print(head(rr_table))
p <- plot_forest(rr_table)
print(p)

##############################################################
# End of Script
##############################################################
