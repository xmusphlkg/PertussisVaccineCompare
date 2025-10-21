
## joinpoint setting for number

run_opt_who = run_options(model="ln",
                          max_joinpoints=5,
                          model_selection_method = 'bic',
                          ci_method = 'parametric',
                          dependent_variable_type = 'crude rate',
                          n_cores=parallel::detectCores())

run_opt_gbd = run_options(model="ln",
                          max_joinpoints=5,
                          model_selection_method = 'bic',
                          ci_method = 'parametric',
                          dependent_variable_type = 'age-adjusted rate',
                          n_cores=parallel::detectCores())

## export options
export_opt_gbd = export_options(aapc_full_range  = FALSE,
                                export_aapc = TRUE,
                                aapc_start_range1 = 2000,
                                aapc_end_range1 = 2004,
                                aapc_start_range2 = 2005,
                                aapc_end_range2 = 2009,
                                aapc_start_range3 = 2010,
                                aapc_end_range3 = 2014)

export_opt_gbd <- paste0(
     export_opt_gbd,
     "\nAAPC Start Range4=2015",
     "\nAAPC End Range4=2019"
)

export_opt_who = export_options(aapc_full_range  = FALSE,
                                export_aapc = TRUE,
                                aapc_start_range1 = 2000,
                                aapc_end_range1 = 2004,
                                aapc_start_range2 = 2005,
                                aapc_end_range2 = 2009,
                                aapc_start_range3 = 2010,
                                aapc_end_range3 = 2014)
export_opt_who <- paste0(
     export_opt_who,
     "\nAAPC Start Range4=2015",
     "\nAAPC End Range4=2019"
)


## fig axis
scientific_10 <- function(x) {
     ifelse(x == 0, 0, parse(text = gsub("[+]", "", gsub("e", "%*%10^", scales::scientific_format()(x)))))
}

## get AAPC from jp_model
get_aapc <- function(jp_model) {
     data <- jp_model$aapc |>
          mutate(across(c(aapc, aapc_c_i_low, aapc_c_i_high), ~formatC(., format = "f", digits = 2)),
                 p_value = as.numeric(p_value),
                 p_value_label = case_when(p_value < 0.001 ~ '***',
                                           p_value < 0.01 ~ '**',
                                           p_value < 0.05 ~ '*',
                                           TRUE ~ ''),
                 legend = paste0(start_obs, '~', end_obs, '\n',
                                 aapc, '(', aapc_c_i_low, '~', aapc_c_i_high, ')', p_value_label),
                 Year = paste(start_obs, end_obs, sep = '~'),
                 Value = paste0(aapc, ' (', aapc_c_i_low, '~', aapc_c_i_high, ')'))
     
     return(data)
}