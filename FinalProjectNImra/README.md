
# multiPathAIC

An R package implementing multi-path AIC model selection (lower AIC)
with variable stability analysis.  
Supports both linear and logistic regression models.

## Installation

\`\`\`r \# If using GitHub \#
remotes::install_github(“yourusername/multiPathAIC”)

# Or load from local folder

devtools::load_all(“path/to/multiPathAIC”)

## Example Usage

library(multiPathAIC)

# Linear model example

set.seed(123) res_linear \<- build_paths(mtcars, response = “mpg”, K =
3, delta = 1) stab_linear \<- stability(mtcars, response = “mpg”, B =
10, K = 3) plausible_models(res_linear, stab_linear, delta_AIC = 2, tau
= 0.6)

# Logistic model example

res_logistic \<- build_paths(mtcars, response = “am”, K = 3, delta = 1)
stab_logistic \<- stability(mtcars, response = “am”, B = 10, K = 3)
plausible_models(res_logistic, stab_logistic, delta_AIC = 2, tau = 0.6)
