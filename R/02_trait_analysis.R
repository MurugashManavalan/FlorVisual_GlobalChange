source(here::here("R", "functions.R"))
install_and_load_all_packages()
df <- read_csv(here::here("Data", "traits_combined_dataset.csv"))
base_traits <- c(5:13)
colnames(df)
# Adding Nutrient and Water Levels
df <- df %>%
  mutate(
    nutrient = if_else(str_detect(Treatment, "^P"), "Low", "High"),
    water = if_else(str_detect(Treatment, "D$"), "Drought", "Well-watered")
  ) %>%
  mutate(
    nutrient = factor(nutrient, levels = c("Low", "High")),
    water = factor(water, levels = c("Well-watered", "Drought"))
  )

# Checking Distribution ----
df <- log_sqrt_transform(df, base_traits)
t_traits <- c(5:13, 16:33) # Transformed Traits
df_shapiro_results <- apply(df[,t_traits], 2, shapiro.test)
df_shapiro_results

# Creating Models and Runing Anova analysis ----
n_traits <- c(5,19,20,8:11,30,13)
models <- run_trait_models(df, n_traits, base_formula = "nutrient * water * Species")
anova_dataset <- create_anova_dataset(models)
write.csv(anova_dataset, "anova_results_allspecies.csv", row.names = FALSE)

# Variance Paritioning Barplot (Including Effect of Species) ----
## Creation of Models and Anova Dataset ----
vp_models <- run_trait_models(df, n_traits, base_formula = "nutrient + water + Species + nutrient:water + nutrient:Species + water:Species + nutrient:water:Species") # Creating needed Models
vp_anova_dataset <- data.frame()
for (model_name in names(vp_models)) {
  model = vp_models[[model_name]]
  anova_results <- anova(model) # Perform ANOVA function
  anova_results_df <- as.data.frame(anova_results) # Convert Anova results into a dataframe
  anova_results_df <- anova_results_df %>%
    mutate(PropVar = `Sum Sq`/ sum(`Sum Sq`),
           Factor = rownames(anova_results_df),
           Trait = model_name) %>%
    filter(`Pr(>F)` < 0.05) %>% # Keep only significant values
    filter(Factor != "(Intercept)" & Factor != "Residuals") # Remove Intercept and Residuals

  vp_anova_dataset <- rbind(vp_anova_dataset, anova_results_df)
}
traits_desired_order <- c("Plant_height_cm","sqrt_No_of_flowers","log_Floral_area_mm", "Bee_chromatic_contrast", "Fly_chromatic_contrast", "Achromatic_contrast","Total_brightness","log_Chroma","Hue")
factors_desired_order <- c("nutrient", "water", "nutrient:water", "nutrient:Species", "water:Species", "nutrient:water:Species", "Species")
vp_anova_dataset$Trait <- factor(vp_anova_dataset$Trait , levels = rev(traits_desired_order))
vp_anova_dataset$Factor <- factor(vp_anova_dataset$Factor, levels = rev(factors_desired_order))

## Creation of Barplot ----
png("Variance Partitioning Plot - 1 September 2026.png", width = 3000, height = 2000, res = 300)
ggplot(vp_anova_dataset, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = NULL, y = "Proportion of Variance Explained",
       title = "Variance Explained by Abiotic Factors") +
  coord_flip() + 
  scale_x_discrete(labels = c(
    Plant_height_cm         = "Plant height (PH)",
    sqrt_No_of_flowers      = "Flower number (FN)",
    log_Floral_area_mm      = "Floral area (FA)",
    Bee_chromatic_contrast  = "Bee chromatic contrast (Bees(CC))",
    Fly_chromatic_contrast  = "Fly chromatic contrast (Flies(CC))",
    Achromatic_contrast     = "Achromatic contrast (ACC)",
    Total_brightness        = "Total brightness (TB)",
    log_Chroma              = "Chroma",
    Hue                     = "Hue"
  )) +
  scale_fill_manual(
    values = c(
      "nutrient" = "#1b9e77",
      "water" = "#7570b3",
      "nutrient:water" = "#d95f02",
      "nutrient:Species" = "#e7298a",
      "water:Species" = "#66a61e",
      "nutrient:water:Species" = "#e6ab02",
      "Species" = "#1f78b4"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 16),
    plot.title  = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank()
  )
dev.off()

# Fitting Linear Models to Remove Species Effects
df$PH_resid <- resid(lm(Plant_height_cm ~ Species, data = df, na.action = na.exclude))
df$FN_resid <- resid(lm(sqrt_No_of_flowers ~ Species, data = df, na.action = na.exclude))
df$FA_resid <- resid(lm(log_Floral_area_mm ~ Species, data = df, na.action = na.exclude))
df$BCC_resid <- resid(lm(Bee_chromatic_contrast ~ Species, data = df, na.action = na.exclude))
df$FCC_resid <- resid(lm(Fly_chromatic_contrast ~ Species, data = df, na.action = na.exclude))
df$ACC_resid <- resid(lm(Achromatic_contrast ~ Species, data = df, na.action = na.exclude))
df$TB_resid <- resid(lm(Total_brightness ~ Species, data = df, na.action = na.exclude))
df$Chroma_resid <- resid(lm(log_Chroma ~ Species, data = df, na.action = na.exclude))
df$Hue_resid <- resid(lm(Hue ~ Species, data = df, na.action = na.exclude))
resid_traits <- c(34:42)

# Tukey Table for Effects of Abiotic factors in All Species
models <- run_trait_models(df, resid_traits, base_formula = "nutrient * water")
tukey_dataset <- create_tukey_dataset(models)
write.csv(tukey_dataset, "tukey_results_allspecies.csv", row.names = FALSE)

## Creating Models and ANOVA for Individual Species ----
species_list <- c("ARA ARE", "DIA DEL", "LEU IRC", "POT ARG", "VER CHA")

results_list <- map(species_list, function(sp) {
  df_sp <- df %>% filter(Species == sp)
  models_sp <- run_trait_models(df_sp, n_traits, base_formula = "nutrient * water")
  
  list(
    anova = create_anova_dataset(models_sp) %>% mutate(Species = sp)  %>% relocate(Species),
    tukey = create_tukey_dataset(models_sp) %>% mutate(Species = sp)  %>% relocate(Species)
  )
})
map_dfr
all_anova <- map_dfr(results_list, "anova")
all_tukey <- map_dfr(results_list, "tukey")
write.csv(all_anova, "all_species_anova.csv", row.names = FALSE)
write.csv(all_tukey, "all_species_tukey.csv", row.names = FALSE)

# Plots for Effect of Treatments on All Species ----
mean_point <- stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red")
facet_theme <- list(
  facet_wrap(~nutrient, labeller = labeller(nutrient = function(x) paste("Nutrient:", x))),
  theme(axis.title = element_text(size = 36),
        axis.text = element_text(size = 32),
        strip.text = element_text(size = 32, face = "bold"),
        legend.position = "none")
)
violin_box <- function(data, x, y, xlab, ylab, top_padding = 0.05, ...) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_violin(trim = FALSE, fill = "skyblue", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    mean_point +
    scale_y_continuous(expand = expansion(mult = c(0.05, top_padding))) +
    labs(x = xlab, y = ylab) +
    theme(axis.title = element_text(size = 36),
          axis.text  = element_text(size = 32),
          ...)
}

plots <- list(
  W_PH    = violin_box(df, "water",    "PH_resid",    "Water Level",    "Plant Height (cm)"),
  N_FN    = violin_box(df, "nutrient", "FN_resid",    "Nutrient Level", "No. of Flowers"),
  N_FA    = violin_box(df, "nutrient", "FA_resid",    "Nutrient Level", "Floral Area (mm²)"),
  N_FCC   = violin_box(df, "nutrient", "FCC_resid",   "Nutrient Level", "Chromatic Contrast (For Flies)"),
  N_ACC   = violin_box(df, "nutrient", "ACC_resid",   "Nutrient Level", "Achromatic Contrast"),
  N_Chroma = violin_box(df, "nutrient", "Chroma_resid", "Nutrient Level", "Chroma"),
  W_Chroma = violin_box(df, "water",   "Chroma_resid", "Water Level",   "Chroma"),
  WN_Chroma = violin_box(df, "water",   "Chroma_resid",  "Water Level",    "Chroma", top_padding = 0.20) + facet_theme,
  W_Hue   = violin_box(df, "water",    "Hue_resid",   "Water Level",    "Hue")
)

for (i in seq_along(plots)) {
  ggsave(filename = paste0(names(plots)[i], ".png"), 
         plot = plots[[i]], 
         width = 12, height = 12, dpi = 300)
}

# Integrating Pollinator Data ----
dfp <- read_csv(here::here("Data", "pollinator_data.csv"))
dfp <- dfp %>%
  mutate(
    nutrient = if_else(str_detect(Target, "^P"), "Low", "High"),
    water = if_else(str_detect(Target,"D$"), "Drought", "Well-watered")
  ) %>%
  mutate(
    N_nutrient = if_else(str_detect(Neighbour, "^P"), "Low", "High"),
    N_water = if_else(str_detect(Neighbour,"D$"), "Drought", "Well-watered") 
  ) %>%
  mutate(
    nutrient = factor(nutrient, levels = c("Low", "High")),
    water = factor(water, levels = c("Well-watered", "Drought")),
    N_nutrient = factor(N_nutrient, levels = c("Low", "High")),
    N_water = factor(N_water, levels = c("Well-watered", "Drought")),
    Species = factor(Species),
    Temperature = factor(Temperature)
  )


## Creating Generalized Linear Mixel Model for Pollinators ----
pv_model <- glmer.nb(data = dfp, Pollinator ~ nutrient * water * Species + Cloud + W_avg + Temperature + N_nutrient + N_water + (1|Block))
pollinator_anova <- Anova(pv_model, type = "II")
Anova(pv_model, type = "II")
pollinator_anova_df <- as.data.frame(pollinator_anova)
pollinator_anova_df <- pollinator_anova_df %>%
  mutate(
    PropVar = Chisq / sum(Chisq, na.rm = TRUE),
    Factor = rownames(pollinator_anova_df),
    Trait = "Pollinator"
  ) %>%
  filter(`Pr(>Chisq)` < 0.05)


### Integrating Pollinator Model with other models ----
vp_anova_dataset <- bind_rows(vp_anova_dataset, pollinator_anova_df)
traits_desired_order <- c("Plant_height_cm","sqrt_No_of_flowers","log_Floral_area_mm", "Bee_chromatic_contrast", "Fly_chromatic_contrast", "Achromatic_contrast","Total_brightness","log_Chroma","Hue","Pollinator")
factors_desired_order <- c("nutrient", "water", "nutrient:water", "nutrient:Species", "water:Species", "nutrient:water:Species", "Species")
vp_anova_dataset$Trait <- factor(vp_anova_dataset$Trait , levels = rev(traits_desired_order))
vp_anova_dataset$Factor <- factor(vp_anova_dataset$Factor, levels = rev(factors_desired_order))

## Plotting Combined Variance Paritioning plot ----
png("Variance Partitioning Plot (With Pollinator) - 1 September 2026.png", width = 3000, height = 2000, res = 300)
ggplot(vp_anova_dataset, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = NULL, y = "Proportion of Variance Explained", title = "Variance Explained by Abiotic Factors") +
  coord_flip() +
  scale_x_discrete(labels = c(
    Plant_height_cm         = "Plant height (PH)",
    sqrt_No_of_flowers      = "Flower number (FN)",
    log_Floral_area_mm      = "Floral area (FA)",
    Bee_chromatic_contrast  = "Bee chromatic contrast (Bees(CC))",
    Fly_chromatic_contrast  = "Fly chromatic contrast (Flies(CC))",
    Achromatic_contrast     = "Achromatic contrast (ACC)",
    Total_brightness        = "Total brightness (TB)",
    log_Chroma              = "Chroma",
    Hue                     = "Hue",
    Pollinator              = "Pollinator Visitation"
  )) +
  scale_fill_manual(
    values = c(
      "nutrient" = "#1b9e77",
      "water" = "#7570b3",
      "nutrient:water" = "#d95f02",
      "nutrient:Species" = "#e7298a",
      "water:Species" = "#66a61e",
      "nutrient:water:Species" = "#e6ab02",
      "Species" = "#1f78b4"
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())
dev.off()

## Plot for Effect of Treatments on Pollinator Visitation in all species ----
dfp$Pv_resid <- resid(lm(Pollinator ~ Species, data = dfp, na.action = na.exclude))
violin_box_pv <- function(data, x, y, title, xlab, ylab, facet = FALSE) {
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_violin(trim = FALSE, fill = "skyblue", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    mean_point +
    labs(title = title, x = xlab, y = ylab) +
    theme(
      plot.title  = element_text(size = 40, face = "italic", hjust = 0.5),
      axis.title  = element_text(size = 36),
      axis.text   = element_text(size = 32)
    )
  if (facet) p <- p + facet_theme
  p
}

p_plots <- list(
  N_Pv  = violin_box_pv(dfp, "nutrient", "Pv_resid", "All species", "Nutrient Level", "Pollinator Visitation"),
  W_Pv  = violin_box_pv(dfp, "water",    "Pv_resid", "All species", "Water Level",    "Pollinator Visitation"),
  WN_Pv = violin_box_pv(dfp, "water",    "Pv_resid", "All species", "Water Level",    "Pollinator Visitation", facet = TRUE)
)

for (i in seq_along(p_plots)) {
  ggsave(filename = paste0(names(p_plots)[i], ".png"), 
         plot = p_plots[[i]], 
         width = 12, height = 12, dpi = 300)
}
## Creating datasets and models for individual species ----
### ARA ARE ----
dfp_araare <- dfp %>% filter(Species == "ARA ARE")
araare_pv_model <- glmer.nb(data = dfp_araare, Pollinator ~ nutrient * water + N_nutrient + N_water + (1|Block))
Anova(araare_pv_model, type = "II")
### DIA DEL ----
dfp_diadel <- dfp %>% filter(Species == "DIA DEL")
diadel_pv_model <- glmer.nb(data = dfp_diadel, Pollinator ~ nutrient * water + N_nutrient + N_water + (1|Block))
Anova(diadel_pv_model, type = "II")
### LEU IRC ----
dfp_leuric <- dfp %>% filter(Species == "LEU IRC")
leuirc_pv_model <- glmer.nb(data = dfp_leuric, Pollinator ~ nutrient * water + N_nutrient + N_water + (1|Block))
Anova(leuirc_pv_model, type = "II")
### POT ARG ----
dfp_potarg <- dfp %>% filter(Species == "POT ARG")
potarg_pv_model <- glmer.nb(data = dfp_potarg, Pollinator ~ nutrient * water + N_nutrient + N_water + (1|Block))
Anova(potarg_pv_model, type = "II")
### VER CHA ----
dfp_vercha <- dfp %>% filter(Species == "VER CHA")
vercha_pv_model <- glmer.nb(data = dfp_vercha, Pollinator ~ nutrient * water + N_nutrient + N_water + (1|Block))
Anova(vercha_pv_model, type = "II")


## Corresponding plos for effect of treatments on individual species ----

pv_plots <- list(
  LEUIRC_N_Pv  = violin_box_pv(dfp_leuric, "nutrient", "Pollinator", "L. ircutianum", "Nutrient Level", "Pollinator Visitation"),
  VERCHA_N_Pv  = violin_box_pv(dfp_vercha, "nutrient", "Pollinator", "V. chamaedrys",  "Nutrient Level", "Pollinator Visitation"),
  POTARG_N_Pv  = violin_box_pv(dfp_potarg, "nutrient", "Pollinator", "P. argentea",    "Nutrient Level", "Pollinator Visitation"),
  POTARG_WN_Pv = violin_box_pv(dfp_potarg, "water",    "Pollinator", "P. argentea",    "Water Level",    "Pollinator Visitation", facet = TRUE),
  ARAARE_N_Pv  = violin_box_pv(dfp_araare, "nutrient", "Pollinator", "A. arenosa",     "Nutrient Level", "Pollinator Visitation"),
  ARAARE_W_Pv  = violin_box_pv(dfp_araare, "water",    "Pollinator", "A. arenosa",     "Water Level",    "Pollinator Visitation")
)

for (i in seq_along(pv_plots)) {
  ggsave(
    filename = paste0(names(pv_plots)[i], ".png"),
    plot = pv_plots[[i]],
    width = 12,
    height = 12,
    dpi = 300
  )
}

# Integrating RDA ----
df_clean <- df[complete.cases(df[, n_traits]), ]
rda.as <- rda(df_clean[, n_traits] ~ nutrient * water * Species, data = df_clean, scale = TRUE) 
RDA_anova <- anova(rda.as, permutations = 9999, by = "terms")

RDA_anova
rda_anova_df <- as.data.frame(RDA_anova) %>%
  rownames_to_column(var = "Factor") %>%   # adds Factor from rownames first
  mutate(
    Variance = as.numeric(Variance),        # ensure numeric
    PropVar = Variance / sum(Variance, na.rm = TRUE),
    Trait = "RDA"
  ) %>%
  filter(`Pr(>F)` < 0.05, Factor != "Residual")
rda_anova_df

## Combined RDA with other models
vp_anova_dataset <- bind_rows(vp_anova_dataset, rda_anova_df)
traits_desired_order <- c("Plant_height_cm","sqrt_No_of_flowers","log_Floral_area_mm", "Bee_chromatic_contrast", "Fly_chromatic_contrast", "Achromatic_contrast","Total_brightness","log_Chroma","Hue","Pollinator", "RDA")
factors_desired_order <- c("nutrient", "water", "nutrient:water", "nutrient:Species", "water:Species", "nutrient:water:Species", "Species")
vp_anova_dataset$Trait <- factor(vp_anova_dataset$Trait , levels = rev(traits_desired_order))
vp_anova_dataset$Factor <- factor(vp_anova_dataset$Factor, levels = rev(factors_desired_order))

## Variance Partitioning Barplot with Pollinators and RDA ----
png("Variance Partitioning Plot - 1 September 2026 - With RDA.png", width = 3000, height = 2000, res = 300)
ggplot(vp_anova_dataset, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = NULL, y = "Proportion of Variance Explained", title = "Variance Explained by Abiotic Factors") +
  coord_flip() + 
  scale_x_discrete(labels = c(
    Plant_height_cm         = "Plant height (PH)",
    sqrt_No_of_flowers      = "Flower number (FN)",
    log_Floral_area_mm      = "Floral area (FA)",
    Bee_chromatic_contrast  = "Bee chromatic contrast (Bees(CC))",
    Fly_chromatic_contrast  = "Fly chromatic contrast (Flies(CC))",
    Achromatic_contrast     = "Achromatic contrast (ACC)",
    Total_brightness        = "Total brightness (TB)",
    log_Chroma              = "Chroma",
    Hue                     = "Hue",
    Pollinator              = "Pollinator Visitation",
    RDA = "RDA"
  )) +
  scale_fill_manual(
    values = c(
      "nutrient" = "#1b9e77",
      "water" = "#7570b3",
      "nutrient:water" = "#d95f02",
      "nutrient:Species" = "#e7298a",
      "water:Species" = "#66a61e",
      "nutrient:water:Species" = "#e6ab02",
      "Species" = "#1f78b4"
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())
dev.off()

## RDA Biplot ----
dfs <- df %>%
  dplyr::rename(
    PH          = Plant_height_cm,
    FN          = sqrt_No_of_flowers,
    FA          = log_Floral_area_mm,
    "Bees(CC)"  = Bee_chromatic_contrast,
    "Flies(CC)" = Fly_chromatic_contrast,
    ACC         = Achromatic_contrast,
    TB          = Total_brightness,
    Chr         = log_Chroma,
    Hue         = Hue
  )
df_clean <- dfs[complete.cases(dfs[, n_traits]), ]

# Grouping traits based on colours
group_colours <- c(
  Morphology          = "#02B050",
  "Pollinator Vision" = "#FF2800",
  "Colour Spectra"    = "#2078b4"
)
trait_groups <- c(
  PH          = "Morphology",
  FN          = "Morphology",
  FA          = "Morphology",
  "Bees(CC)"  = "Pollinator Vision",
  "Flies(CC)" = "Pollinator Vision",
  ACC         = "Pollinator Vision",
  TB          = "Colour Spectra",
  Chr         = "Colour Spectra",
  Hue         = "Colour Spectra"
)
trait_colors <- group_colours[trait_groups]
names(trait_colors) <- names(trait_groups)

rda.asn <- rda(df_clean[, n_traits] ~ nutrient * water + Condition(Species), data = df_clean, scale = TRUE) # Full model

# Open PNG device
png("RDA - 2 September 2026.png", width = 2000, height = 2000, res = 200)

#  Extract scores
site_scores <- scores(rda.asn, display = "sites", scaling = 2)
sp_scores   <- scores(rda.asn, display = "species", scaling = 2)
bp_scores   <- scores(rda.asn, display = "bp", scaling = 2)

# Match colours to sp_scores
sp_colors <- trait_colors[rownames(sp_scores)]

#  Compute axis labels
sum_rda <- summary(rda.asn)
total_imp <- sum_rda$cont$importance # Extract total variance table
axis_names <- colnames(total_imp) # Get axis names exactly as vegan defines them
axis1_name <- axis_names[1] # Extract first two axes
axis2_name <- axis_names[2]
axis1_var <- round(total_imp["Proportion Explained", 1] * 100, 2)
axis2_var <- round(total_imp["Proportion Explained", 2] * 100, 2)
xlab_text <- paste0(axis1_name, " (", axis1_var, "%)")
ylab_text <- paste0(axis2_name, " (", axis2_var, "%)")

#  Automatic arrow scaling
max_site <- max(abs(site_scores))
max_sp   <- max(abs(sp_scores))
arrow_scale <- max_site / max_sp * 0.8  # 80% of site range

#  Determine plot limits (with margin)
all_x <- c(site_scores[,1], sp_scores[,1]*arrow_scale, bp_scores[,1])
all_y <- c(site_scores[,2], sp_scores[,2]*arrow_scale, bp_scores[,2])
xlim <- range(all_x)
ylim <- range(all_y) 

#  Base plot
plot(rda.asn, type = "n", scaling = 2, xlim = xlim, ylim = ylim, axes = FALSE,xlab = xlab_text, ylab = ylab_text)
axis(1, lwd = 2)
axis(2, lwd = 2)
box(lwd = 2)
abline(h = 0, v = 0, lty = 2, col = "grey70")

#  Plot sites
points(site_scores[,1], site_scores[,2], pch = 21, bg = "grey60", col = "black", cex = 1)

#  Plot species arrows (coloured by trait group)
for (i in seq_len(nrow(sp_scores))) {
  arrows(0, 0,
         sp_scores[i,1]*arrow_scale, sp_scores[i,2]*arrow_scale,
         col = sp_colors[i], lwd = 2, length = 0.08)
  text(sp_scores[i,1]*arrow_scale*1.1, sp_scores[i,2]*arrow_scale*1.1,
       labels = rownames(sp_scores)[i], col = sp_colors[i], cex = 1.2, font = 2)
}

#  Plot predictor arrows
arrows(0, 0, bp_scores[,1], bp_scores[,2], col = "black", lwd = 2, length = 0.08, lty = 2)
text(bp_scores[,1]*1.1, bp_scores[,2]*1.1, labels = rownames(bp_scores), col = "black", cex = 1.3, font = 2)

# Legend showing the 3 trait groups
legend("topright", legend = names(group_colours), col = group_colours, lty = 1, lwd = 2, bty = "n", cex = 1.3)

dev.off()

## RDA Biplot (For Nutrient Only) ----
rda.asn <- rda(df_clean[, n_traits] ~ nutrient + Condition(Species, water), data = df_clean, scale = TRUE)

# Open PNG device
png("RDA - Nutrient - 2 September 2026.png", width = 2000, height = 2000, res = 200)

#  Extract scores
site_scores <- scores(rda.asn, display = "sites", scaling = 2)
sp_scores   <- scores(rda.asn, display = "species", scaling = 2)
bp_scores   <- scores(rda.asn, display = "bp", scaling = 2)

# Match colours to sp_scores
sp_colors <- trait_colors[rownames(sp_scores)]

#  Compute axis labels
sum_rda <- summary(rda.asn)
total_imp <- sum_rda$cont$importance # Extract total variance table
axis_names <- colnames(total_imp) # Get axis names exactly as vegan defines them
axis1_name <- axis_names[1] # Extract first two axes
axis2_name <- axis_names[2]
axis1_var <- round(total_imp["Proportion Explained", 1] * 100, 2)
axis2_var <- round(total_imp["Proportion Explained", 2] * 100, 2)
xlab_text <- paste0(axis1_name, " (", axis1_var, "%)")
ylab_text <- paste0(axis2_name, " (", axis2_var, "%)")

#  Automatic arrow scaling
max_site <- max(abs(site_scores))
max_sp   <- max(abs(sp_scores))
arrow_scale <- max_site / max_sp * 0.8  # 80% of site range

#  Determine plot limits (with margin)
all_x <- c(site_scores[,1], sp_scores[,1]*arrow_scale, bp_scores[,1])
all_y <- c(site_scores[,2], sp_scores[,2]*arrow_scale, bp_scores[,2])
xlim <- range(all_x) * 1.1
ylim <- range(all_y) * 1.1

#  Base plot
plot(rda.asn, type = "n", scaling = 2, xlim = xlim, ylim = ylim, axes = FALSE,xlab = xlab_text, ylab = ylab_text)
axis(1, lwd = 2)
axis(2, lwd = 2)
box(lwd = 2)
abline(h = 0, v = 0, lty = 2, col = "grey70")

#  Plot sites
points(site_scores[,1], site_scores[,2], pch = 21, bg = "grey60", col = "black", cex = 1)

#  Plot species arrows (coloured by trait group)
for (i in seq_len(nrow(sp_scores))) {
  arrows(0, 0,
         sp_scores[i,1]*arrow_scale, sp_scores[i,2]*arrow_scale,
         col = sp_colors[i], lwd = 2, length = 0.08)
  text(sp_scores[i,1]*arrow_scale*1.1, sp_scores[i,2]*arrow_scale*1.1,
       labels = rownames(sp_scores)[i], col = sp_colors[i], cex = 1.2, font = 2)
}

#  Plot predictor arrows
arrows(0, 0, bp_scores[,1], bp_scores[,2], col = "black", lwd = 2, length = 0.08, lty = 2)
text(bp_scores[,1]*1.1, bp_scores[,2]*1.1, labels = rownames(bp_scores), col = "black", cex = 1.3, font = 2)

# Legend showing the 3 trait groups
legend("topright", legend = names(group_colours), col = group_colours, lty = 1, lwd = 2, bty = "n", cex = 1.3)

dev.off()

## RDA Biplot (For Water Only) ----
rda.asn <- rda(df_clean[, n_traits] ~ water + Condition(Species, nutrient), data = df_clean, scale = TRUE)

# Open PNG device
png("RDA - Water - 2 September 2026.png", width = 2000, height = 2000, res = 200)

#  Extract scores
site_scores <- scores(rda.asn, display = "sites", scaling = 2)
sp_scores   <- scores(rda.asn, display = "species", scaling = 2)
bp_scores   <- scores(rda.asn, display = "bp", scaling = 2)

# Match colours to sp_scores
sp_colors <- trait_colors[rownames(sp_scores)]

#  Compute axis labels
sum_rda <- summary(rda.asn)
total_imp <- sum_rda$cont$importance # Extract total variance table
axis_names <- colnames(total_imp) # Get axis names exactly as vegan defines them
axis1_name <- axis_names[1] # Extract first two axes
axis2_name <- axis_names[2]
axis1_var <- round(total_imp["Proportion Explained", 1] * 100, 2)
axis2_var <- round(total_imp["Proportion Explained", 2] * 100, 2)
xlab_text <- paste0(axis1_name, " (", axis1_var, "%)")
ylab_text <- paste0(axis2_name, " (", axis2_var, "%)")

#  Automatic arrow scaling
max_site <- max(abs(site_scores))
max_sp   <- max(abs(sp_scores))
arrow_scale <- max_site / max_sp * 0.8  # 80% of site range

#  Determine plot limits (with margin)
all_x <- c(site_scores[,1], sp_scores[,1]*arrow_scale, bp_scores[,1])
all_y <- c(site_scores[,2], sp_scores[,2]*arrow_scale, bp_scores[,2])
xlim <- range(all_x) * 1.1
ylim <- range(all_y) * 1.1

#  Base plot
plot(rda.asn, type = "n", scaling = 2, xlim = xlim, ylim = ylim, axes = FALSE,xlab = xlab_text, ylab = ylab_text)
axis(1, lwd = 2)
axis(2, lwd = 2)
box(lwd = 2)
abline(h = 0, v = 0, lty = 2, col = "grey70")

#  Plot sites
points(site_scores[,1], site_scores[,2], pch = 21, bg = "grey60", col = "black", cex = 1)

#  Plot species arrows (coloured by trait group)
for (i in seq_len(nrow(sp_scores))) {
  arrows(0, 0,
         sp_scores[i,1]*arrow_scale, sp_scores[i,2]*arrow_scale,
         col = sp_colors[i], lwd = 2, length = 0.08)
  text(sp_scores[i,1]*arrow_scale*1.1, sp_scores[i,2]*arrow_scale*1.1,
       labels = rownames(sp_scores)[i], col = sp_colors[i], cex = 1.2, font = 2)
}

#  Plot predictor arrows
arrows(0, 0, bp_scores[,1], bp_scores[,2],col = "black", lwd = 2, length = 0.08, lty = 2)
text(bp_scores[,1]*1.1,bp_scores[,2]*1.1,labels = rownames(bp_scores),col = "black", cex = 1.3, font = 2)

# Legend showing the 3 trait groups
legend("topright", legend = names(group_colours), col = group_colours, lty = 1, lwd = 2, bty = "n", cex = 1.3)

dev.off()

# Correlogram (Which includes p values and central diagonal is white) ----
cor_data <- df[,resid_traits]
names(cor_data)
names(cor_data)[names(cor_data) == "PH_resid"] <- "PH"
names(cor_data)[names(cor_data) == "FN_resid"] <- "FNo"
names(cor_data)[names(cor_data) == "FA_resid"] <- "FA"
names(cor_data)[names(cor_data) == "BCC_resid"] <- "BCC"
names(cor_data)[names(cor_data) == "FCC_resid"] <- "FCC"
names(cor_data)[names(cor_data) == "ACC_resid"] <- "ACC"
names(cor_data)[names(cor_data) == "TB_resid"] <- "TB"
names(cor_data)[names(cor_data) == "Chroma_resid"] <- "Chroma"
names(cor_data)[names(cor_data) == "Hue_resid"] <- "Hue"

## Calculate correlation matrix
corr_matrix <- rcorr(as.matrix(cor_data))
r_values <- corr_matrix$r
p_values <- corr_matrix$P

## Set upper triangle and diagonal elements to NA for r values
r_values[upper.tri(r_values)] <- NA
diag(r_values) <- NA

## Set lower triangle elements to NA for p values
p_values[lower.tri(p_values)] <- NA

## Melt data into long format
r_melt <- melt(r_values, na.rm = TRUE)
p_melt <- melt(p_values, na.rm = TRUE)

## Merge correlation and p_value data
merged_data <- merge(r_melt, p_melt, by = c("Var1", "Var2"), all = TRUE)
colnames(merged_data) <- c("Var1", "Var2", "correlation", "p_value")

## Create labels and fontface based on correlation and p values
merged_data$label <- with(merged_data,
                          ifelse(is.na(correlation),
                                ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value)),
                                round(correlation,2)))
merged_data$fontface <- with(merged_data,
                            ifelse(is.na(correlation) & p_value <= 0.05, "bold", "plain"))

## Plot correlogram
ggplot(data = merged_data, aes(Var1, Var2, fill = correlation)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "red", high = "blue4", mid = "#ffffff",
                            midpoint = 0, limit = c(-1,1), space = "Lab",
                            name = "Pearson\nCorrelation")+
        geom_text(aes(label = label, fontface = fontface), color = "black", size = 5) +
        theme_minimal()+
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks = element_blank()) +
        coord_fixed()
  
## Export correlogram
ggsave("correlation_heatmap - 20 May 2026.jpeg",
      plot = ggplot(data = merged_data, aes(Var1, Var2, fill = correlation)) +
                geom_tile(color = "white") +
                scale_fill_gradient2(low = "red", high = "blue4", mid = "#ffffff",
                      midpoint = 0, limit = c(-1,1), space = "Lab",
                      name = "Pearson\nCorrelation")+
                geom_text(aes(label = label, fontface = fontface), color = "black", size = 5) +
                theme_minimal()+
                theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.ticks = element_blank()) +
              coord_fixed(),
          width = 12, height = 12, dpi = 300 )

# Plots for Effect of Treatments on Individual Species ----
df_araare <- df %>% filter(Species == "ARA ARE")
df_diadel <- df %>% filter(Species == "DIA DEL")
df_leuirc <- df %>% filter(Species == "LEU IRC")
df_potarg <- df %>% filter(Species == "POT ARG")
df_vercha <- df %>% filter(Species == "VER CHA")

plots <- list(
  # ── ARAARE ──────────────────────────────────────────────────────────────────
  ARAARE_WN_PH     = violin_box(df_araare, "water",   "Plant_height_cm",        "Water Level",    "Plant Height (cm)", top_padding = 0.20) + facet_theme,
  ARAARE_N_FN      = violin_box(df_araare, "nutrient","sqrt_No_of_flowers",      "Nutrient Level", "Flower Number"),
  ARAARE_WN_FN     = violin_box(df_araare, "water",   "sqrt_No_of_flowers",      "Water Level",    "Flower Number", top_padding = 0.20) + facet_theme,
  ARAARE_N_FA      = violin_box(df_araare, "nutrient","log_Floral_area_mm",  "Nutrient Level", "Floral Area"),
  ARAARE_W_FA      = violin_box(df_araare, "water",   "log_Floral_area_mm",  "Water Level",    "Floral Area"),
  ARAARE_WN_FA     = violin_box(df_araare, "water",   "log_Floral_area_mm",  "Water Level",    "Floral Area", top_padding = 0.20) + facet_theme,
  ARAARE_N_BCC     = violin_box(df_araare, "nutrient","Bee_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Bees)"),
  ARAARE_WN_BCC    = violin_box(df_araare, "water",   "Bee_chromatic_contrast",  "Water Level",    "Chromatic Contrast (For Bees)", top_padding = 0.20) + facet_theme,
  ARAARE_N_FCC     = violin_box(df_araare, "nutrient","Fly_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Flies)"),
  ARAARE_W_FCC     = violin_box(df_araare, "water",   "Fly_chromatic_contrast",  "Water Level",    "Chromatic Contrast (For Flies)"),
  ARAARE_WN_FCC    = violin_box(df_araare, "water",   "Fly_chromatic_contrast",  "Water Level",    "Chromatic Contrast (For Flies)", top_padding = 0.20) + facet_theme,
  ARAARE_N_ACC     = violin_box(df_araare, "nutrient","Achromatic_contrast",     "Nutrient Level", "Achromatic Contrast"),
  ARAARE_W_ACC     = violin_box(df_araare, "water",   "Achromatic_contrast",     "Water Level",    "Achromatic Contrast"),
  ARAARE_N_TB      = violin_box(df_araare, "nutrient","Total_brightness",        "Nutrient Level", "Total Brightness"),
  ARAARE_W_TB      = violin_box(df_araare, "water",   "Total_brightness",        "Water Level",    "Total Brightness"),
  ARAARE_WN_TB     = violin_box(df_araare, "water",   "Total_brightness",        "Water Level",    "Total Brightness", top_padding = 0.20) + facet_theme,
  ARAARE_W_Chroma  = violin_box(df_araare, "water",   "log_Chroma",              "Water Level",    "Chroma"),
  ARAARE_WN_Chroma = violin_box(df_araare, "water",   "log_Chroma",              "Water Level",    "Chroma", top_padding = 0.20) + facet_theme,
  ARAARE_W_Hue     = violin_box(df_araare, "water",   "Hue",                     "Water Level",    "Hue"),
  ARAARE_WN_Hue    = violin_box(df_araare, "water",   "Hue",                     "Water Level",    "Hue", top_padding = 0.20) + facet_theme,

  # ── DIADEL ──────────────────────────────────────────────────────────────────
  DIADEL_N_PH      = violin_box(df_diadel, "nutrient","Plant_height_cm",         "Nutrient Level", "Plant Height (cm)"),
  DIADEL_W_PH      = violin_box(df_diadel, "water",   "Plant_height_cm",         "Water Level",    "Plant Height (cm)"),
  DIADEL_N_FN      = violin_box(df_diadel, "nutrient","sqrt_No_of_flowers",      "Nutrient Level", "Flower Number"),
  DIADEL_W_FN      = violin_box(df_diadel, "water","sqrt_No_of_flowers",         "Water Level", "Flower Number"),
  DIADEL_WN_FN     = violin_box(df_diadel, "water",   "sqrt_No_of_flowers",      "Water Level",    "Flower Number", top_padding = 0.20) + facet_theme,
  DIADEL_N_FA      = violin_box(df_diadel, "nutrient","log_Floral_area_mm",  "Nutrient Level", "Floral Area"),
  DIADEL_N_BCC     = violin_box(df_diadel, "nutrient","Bee_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Bees)"),
  DIADEL_N_FCC     = violin_box(df_diadel, "nutrient","Fly_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Flies)"),
  DIADEL_N_ACC     = violin_box(df_diadel, "nutrient","Achromatic_contrast",     "Nutrient Level", "Achromatic Contrast"),
  DIADEL_N_TB      = violin_box(df_diadel, "nutrient","Total_brightness",        "Nutrient Level", "Total Brightness"),
  DIADEL_WN_TB     = violin_box(df_diadel, "water",   "Total_brightness",        "Water Level",    "Total Brightness", top_padding = 0.20) + facet_theme,
  DIADEL_N_Chroma  = violin_box(df_diadel, "nutrient","log_Chroma",              "Nutrient Level", "Chroma"),

  # ── LEUIRC ──────────────────────────────────────────────────────────────────
  LEUIRC_WN_PH     = violin_box(df_leuirc, "water",   "Plant_height_cm",         "Water Level",    "Plant Height (cm)", top_padding = 0.20) + facet_theme,
  LEUIRC_N_FN      = violin_box(df_leuirc, "nutrient","sqrt_No_of_flowers",      "Nutrient Level", "Flower Number"),
  LEUIRC_N_BCC     = violin_box(df_leuirc, "nutrient","Bee_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Bees)"),
  LEUIRC_WN_BCC    = violin_box(df_leuirc, "water",   "Bee_chromatic_contrast",  "Water Level",    "Chromatic Contrast (For Bees)", top_padding = 0.20) + facet_theme,
  LEUIRC_N_FCC     = violin_box(df_leuirc, "nutrient","Fly_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Flies)"),
  LEUIRC_WN_FCC    = violin_box(df_leuirc, "water",   "Fly_chromatic_contrast",  "Water Level",    "Chromatic Contrast (For Flies)", top_padding = 0.20)+ facet_theme,
  LEUIRC_N_Chroma  = violin_box(df_leuirc, "nutrient","log_Chroma",              "Nutrient Level", "Chroma"),
  LEUIRC_WN_Chroma = violin_box(df_leuirc, "water",   "log_Chroma",              "Water Level",    "Chroma", top_padding = 0.20) + facet_theme,

  # ── POTARG ──────────────────────────────────────────────────────────────────
  POTARG_N_FN      = violin_box(df_potarg, "nutrient","sqrt_No_of_flowers",      "Nutrient Level", "Flower Number"),
  POTARG_W_FN      = violin_box(df_potarg, "water",   "sqrt_No_of_flowers",      "Water Level",    "Flower Number"),
  POTARG_N_FA      = violin_box(df_potarg, "nutrient","log_Floral_area_mm",  "Nutrient Level", "Floral Area"),
  POTARG_N_ACC     = violin_box(df_potarg, "nutrient","Achromatic_contrast",     "Nutrient Level", "Achromatic Contrast"),
  POTARG_N_TB      = violin_box(df_potarg, "nutrient","Total_brightness",        "Nutrient Level", "Total Brightness"),
  POTARG_N_Chroma  = violin_box(df_potarg, "nutrient","log_Chroma",              "Nutrient Level", "Chroma"),
  POTARG_W_Hue     = violin_box(df_potarg, "water",   "Hue",                     "Water Level",    "Hue"),
  POTARG_WN_Hue    = violin_box(df_potarg, "water",   "Hue",                     "Water Level",    "Hue", top_padding = 0.20) + facet_theme,

  # ── VERCHA ──────────────────────────────────────────────────────────────────
  VERCHA_N_PH      = violin_box(df_vercha, "nutrient","Plant_height_cm",         "Nutrient Level", "Plant Height (cm)"),
  VERCHA_N_FN      = violin_box(df_vercha, "nutrient","sqrt_No_of_flowers",      "Nutrient Level", "Flower Number"),
  VERCHA_N_BCC     = violin_box(df_vercha, "nutrient","Bee_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Bees)"),
  VERCHA_W_BCC     = violin_box(df_vercha, "water",   "Bee_chromatic_contrast",  "Water Level",    "Chromatic Contrast (For Bees)"),
  VERCHA_N_FCC     = violin_box(df_vercha, "nutrient","Fly_chromatic_contrast",  "Nutrient Level", "Chromatic Contrast (For Flies)"),
  VERCHA_W_ACC     = violin_box(df_vercha, "water",   "Achromatic_contrast",     "Water Level",    "Achromatic Contrast"),
  VERCHA_W_TB      = violin_box(df_vercha, "water",   "Total_brightness",        "Water Level",    "Total Brightness"),
  VERCHA_N_Chroma  = violin_box(df_vercha, "nutrient","log_Chroma",              "Nutrient Level", "Chroma"),
  VERCHA_W_Chroma  = violin_box(df_vercha, "water",   "log_Chroma",              "Water Level",    "Chroma"),
  VERCHA_WN_Chroma = violin_box(df_vercha, "water",   "log_Chroma",              "Water Level",    "Chroma", top_padding = 0.20) + facet_theme
)

for (i in seq_along(plots)) {
  ggsave(filename = paste0(names(plots)[i], ".png"), 
         plot = plots[[i]], 
         width = 12, height = 12, dpi = 300)
}

# Merging Pollinator Visits and Trait Dataset ----
dfp <- dfp %>% mutate(Individual_num = str_extract(Individual, "\\d+")) # Extract numeric individual ID
dfp$Individual_num <- as.numeric(dfp$Individual_num)
dfp_merged <- dfp %>%
  left_join(
    df,
    by = c(
      "Species"      = "Species",
      "Target"       = "Treatment",   # match Target in dfp to Treatment in df
      "Individual_num" = "Individual" # match numeric individual ID
    )
  )
dfp_merged$Individual_id <- paste(dfp_merged$Species,dfp_merged$Target,dfp_merged$Individual_num,sep = "_") # Creating unique individual id
write.csv(dfp_merged, "pollination_traits_dataset.csv", row.names = FALSE)
write.csv(df, "traits_residuals_dataset.csv", row.names = FALSE)

