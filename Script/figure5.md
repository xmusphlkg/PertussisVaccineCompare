## Statistical Analysis using Generalized Linear Mixed-Effects Models (GLMM)

To complement the non-parametric analysis and provide a more comprehensive understanding of the factors associated with pertussis incidence, a Generalized Linear Mixed-Effects Model (GLMM) was implemented. This approach allows for the simultaneous assessment of multiple vaccination program characteristics while accounting for the inherent correlation of data within the same country over time and controlling for potential confounding variables.

### 1. Model Rationale

A GLMM is particularly well-suited for this analysis for several key reasons:

-   **Handling of Clustered Data**: The dataset consists of observations from numerous countries, with multiple data points (e.g., yearly incidence rates) from each country. A GLMM addresses the non-independence of these observations by treating `country` as a **random effect**. This accounts for unobserved country-specific factors (e.g., healthcare system quality, reporting standards, population density) that might influence incidence rates.

-   **Appropriate Distribution for Count Data**: Incidence rates are derived from count data (number of cases), which often exhibit a non-normal, right-skewed distribution. The model employed a **Negative Binomial distribution**, which is an extension of the Poisson distribution that is adept at handling **overdispersion**—a common feature of infectious disease counts where the variance is greater than the mean.

-   **Modeling Incidence Rates**: The model directly analyzes incidence rates by including the logarithm of the `Population` as an **offset term**. This effectively models the number of `Cases` while adjusting for the size of the population at risk, thereby analyzing the rate of disease.

-   **Multivariable Adjustment**: Unlike pairwise comparisons, the GLMM framework assesses the association of each vaccination characteristic (e.g., `VaccineDose`) with pertussis incidence while simultaneously adjusting for all other variables in the model (`TimeFirstShotG`, `TimeLastShotG`, `VaccinePregnant`). This provides a more robust estimate of each factor's independent contribution.

### 2. Model Specification

Separate models were fitted for the years 2019, 2021, and 2024 to investigate how the associations may have varied, particularly before and after the onset of the COVID-19 pandemic.

The general structure of the model for each year was as follows:

`Cases ~ VaccineDose + TimeFirstShotG + TimeLastShotG + VaccinePregnant + (1 | location_name) + offset(log(Population))`

-   **Response Variable**: `Cases` (the number of reported pertussis cases).
-   **Fixed Effects**: These are the main explanatory variables of interest:
    -   `VaccineDose`: Number of doses in the national vaccination schedule.
    -   `TimeFirstShotG`: Grouped time of the first vaccine dose.
    -   `TimeLastShotG`: Grouped time of the last vaccine dose.
    -   `VaccinePregnant`: Recommendation status for maternal vaccination.
-   **Random Effect**: `(1 | location_name)` specifies a random intercept for each country, allowing the baseline incidence rate to vary from one country to another.
-   **Offset**: `offset(log(Population))` adjusts for the population size.
-   **Family**: `nbinom2` (Negative Binomial distribution with a quadratic variance-mean relationship).

### 3. Interpretation of Results

The model output was summarized to present the **Rate Ratio (RR)** for each category of the fixed effects. The RR indicates how the incidence rate is expected to change multiplicatively when moving from a reference category to another category, holding all other factors constant.

-   An **RR > 1** suggests an increased risk (higher incidence rate).
-   An **RR < 1** suggests a protective effect (lower incidence rate).
-   An **RR = 1** suggests no association.

The results, including 95% confidence intervals for the RRs and corresponding p-values, were compiled into a summary table for clear interpretation and reporting. This provides a powerful, model-based complement to the visual and non-parametric findings from the initial analysis.
