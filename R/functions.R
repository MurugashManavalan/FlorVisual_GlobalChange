# ---- Function: Install and Load Required Packages ----

install_and_load_all_packages <- function() {
  cran_packages <- c(
    "readr", "lme4", "lmerTest", "ggplot2", "vegan", "dplyr", "corrplot", 
    "igraph", "qgraph", "Hmisc", "ggrepel", "tidyr", "tidyverse","plyr", "ggeffects", "emmeans",
    "car", "reshape2", "factoextra", "purrr", "here", "devtools", "forcats", "MuMIn", "performance", "glmmTMB", "piecewiseSEM"
  )
  
  missing_cran <- cran_packages[!cran_packages %in% installed.packages()[, "Package"]]
  if (length(missing_cran) > 0) {
    install.packages(missing_cran)
  } else {
    message("All CRAN packages are already installed.")
  }
  
  if (!"PCAtest" %in% installed.packages()[, "Package"]) {
    if (!"devtools" %in% installed.packages()[, "Package"]) {
      install.packages("devtools")
    }
    devtools::install_github("arleyc/PCAtest")
  } else {
    message("PCAtest is already installed.")
  }
  
  if (!"ggbiplot" %in% installed.packages()[, "Package"]) {
    if (!"devtools" %in% installed.packages()[, "Package"]) {
      install.packages("devtools")
    }
    devtools::install_github("vqv/ggbiplot")
  } else {
    message("ggbiplot is already installed.")
  }
  
  all_packages <- c(cran_packages, "PCAtest", "ggbiplot")
  invisible(lapply(all_packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
  
  message("All packages installed and loaded successfully.")
}


# ---- Function: Log and Square Root Transform Response Variables ----

log_sqrt_transform <- function(df, cols) {
  if (is.numeric(cols)) {
    cols <- names(df)[cols]
  }
  
  for (col in cols) {
    log_col_name <- paste0("log_", col)
    sqrt_col_name <- paste0("sqrt_", col)
    df[[log_col_name]] <- log(df[[col]])
    df[[sqrt_col_name]] <- sqrt(df[[col]])
  }
  
  return(df)
}

# ---- Function: Create Histogram of Response Variables ----

hist_df <- function(df, cols) {
  if (is.numeric(cols)) {
    cols <- names(df)[cols]
  }
  
  for (col in cols) {
    hist(df[[col]], main = paste0("Histogram of ", col), xlab = col)
  }
}

# ---- Function: To Create Linear Models ----

run_trait_models <- function(data, traits, base_formula) {
  # data: dataframe
  # traits: column names or indices
  # base_formula: right-hand side of the model, as a character string
  #   e.g. "nutrient * water * Species"
  
  # Convert indices to names if needed
  if (is.numeric(traits)) {
    traits <- names(data)[traits]
  }
  
  models <- lapply(traits, function(trait) {
    # Quote the trait name safely (handles spaces, symbols, etc.)
    trait_quoted <- paste0("`", trait, "`")
    
    # Build full formula string and convert to formula
    formula <- as.formula(paste(trait_quoted, "~", base_formula))
    
    # Fit the model
    lm(formula, data = data)
  })
  
  names(models) <- traits
  return(models)
}

# ---- Function: To create ANOVA dataset ----
create_anova_dataset <- function(models) {
  # models: list of lm models, named by trait
  
  library(dplyr) # for bind_rows
  
  anova_list <- lapply(names(models), function(trait_name) {
    aov_table <- anova(models[[trait_name]])
    aov_df <- as.data.frame(aov_table)
    aov_df$Term <- rownames(aov_df)
    aov_df$Trait <- trait_name
    rownames(aov_df) <- NULL
    # Rearrange columns
    aov_df <- aov_df %>% select(Trait, Term, everything())
    return(aov_df)
  })
  
  # Combine all traits into a single dataframe
  anova_dataset <- bind_rows(anova_list)
  return(anova_dataset)
}

# ---- Function: Run ANOVA and return a dataframe ----
tidy_anova <- function(model, species_name) {
  a <- Anova(model, test.statistic = "Chisq", type = "II")
  as.data.frame(a) %>%
    rownames_to_column("term") %>%
    transmute(
      Species   = species_name,
      term      = term,
      df        = Df,
      statistic = Chisq,
      p_value   = `Pr(>Chisq)`
    )
}
# ---- Extracting Variance Partitioning (Random Effects Contribution) ----
extract_varcorr <- function(model) {
  vc <- VarCorr(model)

  tibble(
    Block_var = attr(vc$Block, "stddev")^2,
    Individual_var = attr(vc$Individual_id, "stddev")^2,
    Residual_var = attr(vc, "sc")^2
  )
}

# ---- Function: Create Tukey Dataset ----
create_tukey_dataset <- function(models) {
  imap_dfr(models, function(m, trait_name) {
    tukey_res <- TukeyHSD(aov(m))
    
    imap_dfr(tukey_res, function(term_result, term_name) {
      as.data.frame(term_result) %>%
        tibble::rownames_to_column("comparison") %>%
        mutate(term = term_name)
    }) %>%
      mutate(trait = trait_name)
  }) %>%
    select(trait, term, comparison, diff, lwr, upr, `p adj`)
}