library(tidyverse)
library(ggplot2)
source(here::here("R", "functions.R"))
install_and_load_all_packages()

# Effect of Treatments on Floral Reflectance (Values based on Endler and Mielke, 2005) ----
results <- read_csv(here::here("Data", "reflectance_results.csv"))
results <- results %>% filter(species %in% c("ARA ARE", "DIA DEL", "LEU IRC", "POT ARG", "VER CHA"))
with(results, table(species, treatment))

## Separate nutrient (P/R) and water (D/W)
results <- results %>%
  mutate(
    nutrient = if_else(str_detect(treatment, "^P"), "Low", "High"),
    water = if_else(str_detect(treatment, "D$"), "Drought", "Well-watered")
  ) %>%
  mutate(
    nutrient = factor(nutrient, levels = c("Low", "High")),
    water = factor(water, levels = c("Well-watered", "Drought"))
  )

# Checking Distribution ----
colnames(results)
base_traits <- c(5:10)
results <- log_sqrt_transform(results, base_traits)
t_traits <- c(5:10,13:24) # Transformed Traits
hist_df(results, t_traits)
df_shapiro_results <- apply(results[,t_traits], 2, shapiro.test)
df_shapiro_results

## Linear Models ----
b_chro_contrast <- lm(bee_chromatic_contrast ~ nutrient * water, data = results)
f_chro_contrast <- lm(fly_chromatic_contrast ~ nutrient * water, data = results)
b_achro_contrast <- lm(achromatic_contrast ~ nutrient * water, data = results)
l_tb <- lm(total_brightness ~ nutrient * water, data = results)
l_chroma <- lm(chroma ~ nutrient * water, data = results)
l_hue <- lm(hue ~ nutrient * water, data = results)

## Mixed Models ----
b_chro_contrast <- lmerTest::lmer(log_bee_chromatic_contrast ~ nutrient * water + (1 | individual), data = results)
f_chro_contrast <- lmerTest::lmer(sqrt_fly_chromatic_contrast ~ nutrient * water + (1 | individual), data = results)
b_achro_contrast <- lmerTest::lmer(achromatic_contrast ~ nutrient * water + (1 | individual), data = results)
l_tb <- lmerTest::lmer(sqrt_total_brightness ~ nutrient * water + (1 | individual), data = results)
l_chroma <- lmerTest::lmer(log_chroma ~ nutrient * water + (1 | individual), data = results)
l_hue <- lmerTest::lmer(hue ~ nutrient * water + (1 | individual), data = results)

## Testing Significance ----
anova(b_chro_contrast)
anova(f_chro_contrast)
anova(b_achro_contrast)
anova(l_tb)
anova(l_chroma)
anova(l_hue)

# Fitting Linear Models to Remove Species Effects
results$bcc_resid <- resid(lm(log_bee_chromatic_contrast ~ species, data = results))
results$ac_resid <- resid(lm(achromatic_contrast ~ species, data = results))
results$TB_resid <- resid(lm(sqrt_total_brightness ~ species, data = results))
results$Chroma_resid <- resid(lm(log_chroma ~ species, data = results))
results$Hue_resid <- resid(lm(hue ~ species, data = results))

## Visulization ----
library(ggplot2)

### 1. Nutrient → Bee chromatic contrast
ggplot(results, aes(x = nutrient, y = bee_chromatic_contrast)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(
    title = "Effect of Nutrient on Chromatic Contrast (Bees)",
    x = "Nutrient Treatment",
    y = "Chromatic Contrast (Bees)"
  ) +
  theme_classic(base_size = 14)

### 2. Nutrient x Water → Bee chromatic contrast
ggplot(results, aes(x = water, y = bee_chromatic_contrast, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(
    title = "Effect of Nutrient and Water on Chromatic Contrast (Bees)",
    x = "Water Treatment",
    y = "Chromatic Contrast (Bees)",
    fill = "Nutrient Treatment"
  )
  theme_classic(base_size = 14)

### 3. Water → Achromatic contrast
ggplot(results, aes(x = water, y = achromatic_contrast)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(
    title = "Effect of Drought on Achromatic Contrast",
    x = "Water Treatment",
    y = "Achromatic Contrast"
  ) +
  theme_classic(base_size = 14)

### 4. Nutrient x Water → Achromatic Contrast
ggplot(results, aes(x = water, y = achromatic_contrast, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(
    title = "Effect of Nutrient and Water on Achromatic contrast",
    x = "Water Treatment",
    y = "Achromatic contrast",
    fill = "Nutrient Treatment"
  )
  theme_classic(base_size = 14)

### 5. Nutrient x Water → Total Brightness
ggplot(results, aes(x = water, y = total_brightness, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(
    title = "Effect of Nutrient and Water on Total Brightness",
    x = "Water Treatment",
    y = "Total Brightness",
    fill = "Nutrient Treatment"
  )
  theme_classic(base_size = 14)

### 6. Nutrient x Water → Chroma
ggplot(results, aes(x = water, y = chroma, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(
    title = "Effect of Nutrient and Water on Chroma",
    x = "Water Treatment",
    y = "Chroma",
    fill = "Nutrient Treatment"
  )
  theme_classic(base_size = 14)

### 7. Water → Hue
ggplot(results, aes(x = water, y = hue)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(
    title = "Effect of Drought on Hue",
    x = "Water Treatment",
    y = "Hue"
  ) +
  theme_classic(base_size = 14)

## Visulization (Controlling Species Effects) ----
plots <- list(
### 1. Nutrient → Bee chromatic contrast
N_bcc = ggplot(results, aes(x = nutrient, y = bcc_resid)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(title = "Effect of Nutrient on Chromatic Contrast (Bees)", x = "Nutrient Treatment", y = "Chromatic Contrast (Bees)") +
  theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 32)),

### 2. Nutrient x Water → Bee chromatic contrast
WN_bcc = ggplot(results, aes(x = water, y = bcc_resid, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(title = "Effect of Nutrient and Water on Chromatic Contrast (Bees)", x = "Water Treatment", y = "Chromatic Contrast (Bees)", fill = "Nutrient Treatment")+
  theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 32)),

### 3. Water → Achromatic contrast
W_ac = ggplot(results, aes(x = water, y = ac_resid)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(title = "Effect of Drought on Achromatic Contrast", x = "Water Treatment", y = "Achromatic Contrast") +
  theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 32)),

### 4. Nutrient x Water → Achromatic Contrast
WN_ac = ggplot(results, aes(x = water, y = ac_resid, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(title = "Effect of Nutrient and Water on Achromatic contrast", x = "Water Treatment", y = "Achromatic contrast", fill = "Nutrient Treatment") + 
  theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 32)),

### 5. Nutrient x Water → Total Brightness
WN_TB = ggplot(results, aes(x = water, y = TB_resid, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(title = "Effect of Nutrient and Water on Total Brightness", x = "Water Treatment", y = "Total Brightness", fill = "Nutrient Treatment") +
  theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 32)),

### 6. Nutrient x Water → Chroma
WN_Chroma = ggplot(results, aes(x = water, y = Chroma_resid, fill = nutrient)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  labs(title = "Effect of Nutrient and Water on Chroma", x = "Water Treatment", y = "Chroma", fill = "Nutrient Treatment") + 
  theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 32)),

### 7. Water → Hue
W_Hue = ggplot(results, aes(x = water, y = Hue_resid)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(title = "Effect of Drought on Hue", x = "Water Treatment", y = "Hue") +
  theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 32))
)

for (i in seq_along(plots)) {
  ggsave(filename = paste0(names(plots)[i], ".png"), 
         plot = plots[[i]], 
         width = 12, height = 12, dpi = 300)
}

cor.test(x= results$bee_chromatic_contrast, y = results$total_brightness, method = "pearson")

?cor.test
