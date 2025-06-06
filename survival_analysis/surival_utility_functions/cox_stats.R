cox_stats <- function(model) {
  summary_table          <- summary(model)
  coefficient_table      <- summary_table$coefficients # table with β, se, z, p-value
  confidence_intervals   <- summary_table$conf.int     # table with HR & CI
  variable_levels        <- model$xlevels              # factor levels used in the fit
  
  # initialize an empty vector to store the output
  out_label <- c()
  
  # loop for the number of coefficients
  for (i in seq_len(nrow(coefficient_table))) {
    rn <- rownames(coefficient_table)[i]         
    
    ## Try to map the row name back to “level vs. ref”
    label <- rn  # fall-back if it is a numeric covariate                     
    for (var in names(variable_levels)) {
      if (startsWith(rn, var)) {
        level <- sub(paste0("^", var), "", rn)
        ref   <- variable_levels[[var]][1]
        label <- paste(level, "vs.", ref)
        break
      }
    }
    
    HR     <- round(confidence_intervals[i, "exp(coef)"], 2)
    lower  <- round(confidence_intervals[i, "lower .95"], 1)
    upper  <- round(confidence_intervals[i, "upper .95"], 1)
    pval   <- prettyNum(coefficient_table[i, "Pr(>|z|)"])
    out_label[i] <- sprintf("%s: HR (95%% CI) = %s (%s–%s); p = %s",
                      label, HR, lower, upper, pval)
  }
  
  out_label
}

# example usage -----------------
# 
# library(survival)
# library(survminer)
# library(tidyverse)
# 
# lung.1 <- lung %>%
#   mutate(age_3cl = case_when(
#     age < quantile(age, 0.333) ~ "low",
#     age > quantile(age, 0.667) ~ "high",
#     is.na(age) ~ NA,
#     TRUE ~ "medium"
#   ) %>% factor(levels = c("low", "medium", "high")))
# 
# cox.model <- coxph(Surv(time, status) ~ age_3cl, data = lung.1)
# 
# 
# stats_label <- cox_stats(cox.model) %>%
#   paste0(collapse = "\n")
