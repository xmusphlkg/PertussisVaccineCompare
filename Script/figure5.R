#####################################
## @Description: Figure 5 - GLMM analysis for pertussis incidence vs vaccination program
## @version: 1.0.0
## @Author: Li Kangguo
## @Date: 2025-10-22 11:08:31
## @LastEditors: Li Kangguo
## @LastEditTime: 2025-10-22 11:40:00
#####################################

library(tidyverse)
library(glmmTMB)
library(gtsummary)
library(broom.mixed)
library(openxlsx)
library(flextable)
library(gt)
library(DHARMa)

rm(list = ls())

# Load precomputed DataAll (created by Script/figure4.R)
if(!file.exists('./Output/DataAll.RData')) stop('DataAll.RData not found. Run Script/figure4.R first to create ./Output/DataAll.RData')
load('./Output/DataAll.RData')

DataLongCounts <- DataAll |>
	select(location_id, starts_with('Cases_'), starts_with('Population_'), VaccineDose, TimeFirstShotG, TimeLastShotG, VaccinePregnant) |>
	pivot_longer(cols = starts_with('Cases_'), names_to = 'CasesVar', values_to = 'Cases') |>
	mutate(Year = as.integer(str_extract(CasesVar, '\\d{4}')))

# attach population
DataLongCounts <- DataLongCounts |>
	mutate(PopVar = paste0('Population_', Year)) |>
	rowwise() |>
	mutate(Population = as.numeric(cur_data_all()[[PopVar]])) |>
	ungroup() |>
	select(-PopVar, -CasesVar)

# prepare predictors and filter
DataLongCounts <- DataLongCounts |>
	filter(Year %in% c(2019,2021,2024)) |>
	mutate(
		VaccineDose = factor(VaccineDose),
		TimeFirstShotG = factor(TimeFirstShotG),
		TimeLastShotG = factor(TimeLastShotG),
		VaccinePregnant = factor(VaccinePregnant)
	) |>
	filter(!is.na(Cases) & !is.na(Population) & Population > 0)

# modelling function for counts (negative binomial)
fit_nb_glmm_counts <- function(df, year) {
	dsub <- df |> filter(Year == year)
	if(nrow(dsub) < 10) { message('Too few obs for ', year); return(NULL) }
	fit <- tryCatch({
		glmmTMB(Cases ~ VaccineDose + TimeFirstShotG + TimeLastShotG + VaccinePregnant + (1 | location_id),
						data = dsub,
						offset = log(Population),
						family = nbinom2,
						ziformula = ~0)
	}, error = function(e) { message('fit error ', year, ': ', e$message); return(NULL) })
	return(fit)
}

years_to_fit <- c(2019,2021,2024)
models <- list()
for(y in years_to_fit) models[[as.character(y)]] <- fit_nb_glmm_counts(DataLongCounts, y)

# tidy and exponentiate coefficients (rate ratios)
result_rows <- list()
for(y in names(models)) {
	m <- models[[y]]
	if(is.null(m)) next
	tidy_m <- broom.mixed::tidy(m, effects = 'fixed', conf.int = TRUE, exponentiate = TRUE)
	tidy_m <- tidy_m |> mutate(year = y)
	result_rows[[y]] <- tidy_m
}

results_all <- bind_rows(result_rows) |>
	select(year, term, estimate, conf.low, conf.high, p.value, std.error)

# Create gtsummary tables
tbls <- list()
for(y in names(models)) {
	if(is.null(models[[y]])) next
	tbls[[y]] <- tbl_regression(models[[y]], exponentiate = has_counts, estimate_fun = ~style_ratio(.x, digits = 2)) |> modify_caption(paste0('GLMM - ', y))
}
if(length(tbls) > 0) {
	merged_tbl <- tbl_merge(tbls = tbls)
	try({ merged_tbl |> as_flex_table() |> flextable::save_as_docx(path = './Output/Figure 5_GLMM_Results.docx') }, silent = TRUE)
	try({ merged_gt <- merged_tbl |> as_gt(); gtsave(merged_gt, filename = './Output/Figure 5_GLMM_Results.png') }, silent = TRUE)
}

## Enhanced diagnostics and plots
try({
	diag_list <- list()
	for(y in names(models)) {
		m <- models[[y]]
		if(is.null(m)) next

		# Data used
		dsub <- DataLongCounts |> filter(Year == as.integer(y))
		n_obs <- nrow(dsub)

		# Pearson dispersion statistic
		res_p <- residuals(m, type = 'pearson')
		p_len <- length(res_p)
		p_fixef <- length(fixef(m)$cond)
		df_resid <- max(1, p_len - p_fixef)
		pearson_chisq <- sum(res_p^2, na.rm = TRUE)
		dispersion <- pearson_chisq / df_resid

		# DHARMa simulated residuals and dispersion test (if available)
		dh <- tryCatch({
			sim <- simulateResiduals(m, plot = FALSE)
			testDispersion(sim)
		}, error = function(e) { NULL })
		dh_p <- if(!is.null(dh)) dh$p.value else NA_real_

		# random effect variance
		vc <- tryCatch({ as.data.frame(VarCorr(m)$cond) }, error = function(e) NULL)
		re_var <- NA_real_
		if(!is.null(vc) && 'vcov' %in% names(vc)) {
			# pick location_id group if present, else first
			idx <- which(vc$grp == 'location_id')
			if(length(idx) == 0) idx <- 1
			re_var <- vc$vcov[idx[1]]
		}

		# AIC / logLik
		aic_val <- AIC(m)
		loglik_val <- as.numeric(logLik(m))

		diag_list[[y]] <- tibble(year = y, n = n_obs, pearson_chisq = pearson_chisq, df_resid = df_resid, dispersion = dispersion, DHARMa_p = dh_p, re_var = re_var, AIC = aic_val, logLik = loglik_val)
	}

	diag_df <- bind_rows(diag_list)
	# write diagnostics to Excel
	addWorksheet(wb, 'ModelDiagnostics')
	writeData(wb, sheet = 'ModelDiagnostics', x = diag_df)
	saveWorkbook(wb, './Output/Figure 5_GLMM_Results.xlsx', overwrite = TRUE)

	# Enhance forest plot: color by significance and add RR label
	plot_df <- results_all |> rename(estimate = estimate, lower = conf.low, upper = conf.high)
	# p-value may be NA for some models
	plot_df <- plot_df |> mutate(pval = ifelse(is.na(p.value), 1, p.value), signif = case_when(pval < 0.001 ~ '***', pval < 0.01 ~ '**', pval < 0.05 ~ '*', TRUE ~ 'ns'))
	plot_df <- plot_df |> filter(!term %in% c('(Intercept)'))
	plot_df <- plot_df |> mutate(
		term_clean = term |> str_replace_all('VaccineDose', 'Dose: ') |> str_replace_all('TimeFirstShotG', 'FirstShot: ') |> str_replace_all('TimeLastShotG', 'LastShot: ') |> str_replace_all('VaccinePregnant', 'MaternalVac: '),
		term_label = ifelse(nchar(term_clean) > 50, paste0(substr(term_clean,1,47),'...'), term_clean)
	)

	p_forest <- ggplot(plot_df, aes(x = estimate, y = fct_rev(fct_inorder(term_label)), color = signif)) +
		geom_point(aes(size = 1/std.error), show.legend = TRUE) +
		geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
		facet_wrap(~ year, scales = 'free_x') +
		scale_x_log10(labels = scales::label_number(scale = 1, accuracy = 0.01)) +
		scale_color_manual(values = c('***' = '#d73027', '**' = '#fc8d59', '*' = '#fee08b', 'ns' = '#999999')) +
		guides(size = 'none') +
		labs(x = 'Rate ratio (log scale)', y = NULL, color = 'p-value', title = 'GLMM: Rate ratios (95% CI) by year') +
		theme_minimal() + theme(axis.text.y = element_text(size = 9))

	ggsave('./Output/Figure5_GLMM_Forest.png', plot = p_forest, width = 10, height = 7, dpi = 300)
	ggsave('./Output/Figure5_GLMM_Forest.pdf', plot = p_forest, width = 10, height = 7, device = cairo_pdf)
	message('Saved enhanced forest plot and diagnostics to ./Output/')

}, silent = FALSE)

