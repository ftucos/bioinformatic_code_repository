get.cox.stats <- function(model, test = NULL) {
  x_variables            <- model$xlevels                    # predictors
  summary_table          <- summary(model)
  coefficient_table      <- summary_table$coefficients       # table with β, se, z, p-value
  confidence_intervals   <- summary_table$conf.int           # table with HR & CI
  variable_levels        <- model$xlevels                    # factor levels used in the fit
  logrank_pval           <- summary_table$sctest["pvalue"]   # logrank p-value
  lrt_pval               <- summary_table$logtest["pvalue"]  # likelihood ratio test p-value
  wald_global_pval       <- summary_table$waldtest["pvalue"] # Wald test p-value
  
  # check if the model is multivariable
  multivariable = length(x_variables) > 1
  
  # define the significance test to use --------------------------------------
  # default: if user didn't specify 'test', choose Wald for multilevel, else logrank
  if (is.null(test)) {
    test <- if (nrow(summary(model)$coefficients) > 1) "wald" else "logrank"
  }
  
  # throw a warning if the user chosed logrank or lrt for a multilevel model
  if (test %in% c("logrank", "lrt", "wald_global") && nrow(coefficient_table) > 1) {
    warning(sprintf(
      "You've requested the %s test on a model with %d levels: only a single overall p-value will be returned.",
      test, nrow(coefficient_table)+1
    ))
  }
  
  # match arguments
  test <- match.arg(test, choices = c("wald", "logrank", "lrt", "wald_global"))
  
  # build labels ---------------------------------
  
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
        
        # if the model is multivariable, prepone the variable name
        if (multivariable) {
          label <- paste(var, level, "vs.", ref)
        } else {
          label <- paste(level, "vs.", ref)
        }
        break
      }
    }
    
    HR     <- round(confidence_intervals[i, "exp(coef)"], 2)
    lower  <- round(confidence_intervals[i, "lower .95"], 1)
    upper  <- round(confidence_intervals[i, "upper .95"], 1)
    pval   <- prettyNum(switch(test,
                               "wald" = coefficient_table[i, "Pr(>|z|)"],
                               "logrank" = logrank_pval,
                               "lrt" = lrt_pval,
                               "wald_global" = wald_global_pval
    ), digits = 2)
    out_label[i] <- sprintf("%s: HR (95%% CI) = %s (%s–%s); p = %s",
                            label, HR, lower, upper, pval)
  }
  
  return(out_label)
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
