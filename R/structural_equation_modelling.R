source(here::here("R", "functions.R"))
install_and_load_all_packages()
dft <- read_csv(here::here("Data", "traits_combined_dataset.csv"))
dftr <- read_csv(here::here("Data", "traits_residuals_dataset.csv"))
dfp <- read_csv(here::here("Data", "pollinator_traits_dataset.csv"))
# Preparing Dataframe ----
# Adding Nutrient and Water Levels
dfp <- dfp %>%
  mutate(
    nutrient = if_else(str_detect(Target, "^P"), 0, 1),
    water = if_else(str_detect(Target, "W$"), 0, 1)
  )
dft <- dft %>%
  mutate(
    nutrient = if_else(str_detect(Treatment, "^P"), 0, 1),
    water = if_else(str_detect(Treatment, "W$"), 0, 1)
  )
dftr <- dftr %>%
  mutate(
    nutrient = if_else(str_detect(Treatment, "^P"), 0, 1),
    water = if_else(str_detect(Treatment, "W$"), 0, 1)
  )
# Transforming Traits 
base_traits <- c(5:13)
dft <- log_sqrt_transform(dft, base_traits)

# Removing NAs
## Traits Dataset
vars_needed <- c("nutrient", "water", "Species", "Plant_height_cm", "sqrt_No_of_flowers",
                    "log_Floral_area_mm", "Bee_chromatic_contrast", "Fly_chromatic_contrast", "Achromatic_contrast",
                    "Total_brightness", "log_Chroma", "Hue")
dft_complete <- dft[complete.cases(dft[, vars_needed]), ]

## Pollinator Dataset
vars_needed <- c( "Pollinator", "nutrient", "water", "Species", "Plant_height_cm", "sqrt_No_of_flowers",
                    "log_Floral_area_mm", "Bee_chromatic_contrast", "Fly_chromatic_contrast", "Achromatic_contrast",
                    "Total_brightness", "log_Chroma", "Hue")
dfp_complete <- dfp[complete.cases(dfp[, vars_needed]), ]

## Traits Residuals Dataset
traitsresiduals_vars_needed <- c("PH_resid", "FN_resid", "FA_resid", "BCC_resid", "FCC_resid", "ACC_resid", "TB_resid", "Chroma_resid",
                  "Hue_resid")
dftr_complete <- dftr[complete.cases(dftr[, traitsresiduals_vars_needed]), ]

# Scaling Traits
## Traits Dataset
dft_scaled <- dft_complete
cont_vars <- c("Plant_height_cm", "sqrt_No_of_flowers", "log_Floral_area_mm", "Bee_chromatic_contrast", "Fly_chromatic_contrast",
                "Achromatic_contrast", "Total_brightness", "log_Chroma", "Hue")
dft_scaled[cont_vars] <- lapply(dft_scaled[cont_vars], scale)
dft_scaled <- as.data.frame(dft_scaled)
dft_scaled[] <- lapply(dft_scaled, function(x) { if (is.matrix(x) || is.array(x)) return(as.numeric(x)) 
                      return(x)})

## Pollinator Dataset
dfp_scaled <- dfp_complete
cont_vars <- c("Plant_height_cm", "sqrt_No_of_flowers", "log_Floral_area_mm", "Bee_chromatic_contrast", "Fly_chromatic_contrast",
                "Achromatic_contrast", "Total_brightness", "log_Chroma", "Hue", "PH_resid", "FN_resid", "FA_resid", "BCC_resid", 
                "FCC_resid", "ACC_resid", "TB_resid", "Chroma_resid", "Hue_resid")
dfp_scaled[cont_vars] <- lapply(dfp_scaled[cont_vars], scale)
dfp_scaled <- as.data.frame(dfp_scaled)
dfp_scaled[] <- lapply(dfp_scaled, function(x) { if (is.matrix(x) || is.array(x)) return(as.numeric(x)) 
                      return(x)})

## Traits Residuals Dataset
dftr_scaled <- dftr_complete
dftr_scaled[traitsresiduals_vars_needed] <- lapply(dftr_scaled[traitsresiduals_vars_needed], scale)
dftr_scaled <- as.data.frame(dftr_scaled)
dftr_scaled[] <- lapply(dftr_scaled, function(x) { if (is.matrix(x) || is.array(x)) return(as.numeric(x)) 
                      return(x)})

# SEM Model Preparation (All species)----
### Upstream ----
PH_model <- lm(PH_resid ~ nutrient, data = dftr_scaled)
FN_model <- lm(FN_resid ~ nutrient + water, data = dftr_scaled) 
FA_model <- lm(FA_resid ~ nutrient, data = dftr_scaled)
ACC_model <- lm(ACC_resid ~ nutrient, data = dftr_scaled) 
TB_model<- lm(TB_resid ~ nutrient, data = dftr_scaled) 
CH_model <- lm(Chroma_resid ~ nutrient * water, data = dftr_scaled)
H_model <- lm(Hue_resid ~ water, data = dftr_scaled)

### Downstream ----
PV_model <- lmer(Pv_resid ~ nutrient * water + FA_resid + BCC_resid + FCC_resid + ACC_resid + 
                            Chroma_resid + (1|Block), data = dfp_scaled)
Anova(PV_model, type = "II")

### SEM ----
allspecies_semmod <- piecewiseSEM::psem(
  PH_model,FN_model, FA_model, ACC_model, 
  TB_model, CH_model, H_model, PV_model,
  PH_resid %~~% FN_resid,
  PH_resid %~~% FA_resid,
  PH_resid %~~% Hue_resid,
  FN_resid %~~% FA_resid,
  FN_resid %~~% BCC_resid,
  FN_resid %~~% FCC_resid,
  FN_resid %~~% TB_resid,
  FN_resid %~~% Hue_resid,
  FA_resid %~~% Hue_resid,
  BCC_resid %~~% FCC_resid,
  BCC_resid %~~% ACC_resid,
  BCC_resid %~~% TB_resid,
  BCC_resid %~~% Chroma_resid,
  BCC_resid %~~% Hue_resid,
  FCC_resid %~~% Chroma_resid,
  FCC_resid %~~% Hue_resid,
  ACC_resid %~~% TB_resid,
  ACC_resid %~~% Hue_resid,
  TB_resid %~~% Hue_resid,
  Hue_resid %~~% Chroma_resid,
  data = dfp_scaled
)
dSep(allspecies_semmod)
fisherC(allspecies_semmod)
piecewiseSEM::coefs(allspecies_semmod, standardize = "none")
allspecies_rsquared <- piecewiseSEM::rsquared(allspecies_semmod)
allspecies_coef <- piecewiseSEM::coefs(allspecies_semmod, standardize = "none")
write.csv(allspecies_coef, "allspecies_sem-coefficients.csv", row.names = FALSE)
write.csv(allspecies_rsquared, "allspecies_rsquared.csv", row.names = FALSE)

# Individual Species dataframes ----
## Trait Dataset 
dft_araare <- dft_scaled %>% filter(Species == "ARA ARE")
dft_diadel <- dft_scaled %>% filter(Species == "DIA DEL")
dft_leuirc <- dft_scaled %>% filter(Species == "LEU IRC")
dft_potarg <- dft_scaled %>% filter(Species == "POT ARG")
dft_vercha <- dft_scaled %>% filter(Species == "VER CHA")

## Pollinator Dataset
dfp_araare <- dfp_scaled %>% filter(Species == "ARA ARE")
dfp_diadel <- dfp_scaled %>% filter(Species == "DIA DEL")
dfp_leuirc <- dfp_scaled %>% filter(Species == "LEU IRC")
dfp_potarg <- dfp_scaled %>% filter(Species == "POT ARG")
dfp_vercha <- dfp_scaled %>% filter(Species == "VER CHA")

# SEM Model Preparation ----
## ARA ARE ----
### Upstream ----
PH_model <- lm(Plant_height_cm ~ nutrient * water, data = dft_araare)
FN_model <- lm(sqrt_No_of_flowers ~ nutrient * water, data = dft_araare) 
FA_model <- lm(log_Floral_area_mm ~ nutrient * water, data = dft_araare)
BCC_model <- lm(Bee_chromatic_contrast ~ nutrient * water, data = dft_araare)
FCC_model <- lm(Fly_chromatic_contrast ~ nutrient * water, data = dft_araare)
ACC_model <- lm(Achromatic_contrast ~ nutrient + water, data = dft_araare) 
TB_model<- lm(Total_brightness ~ nutrient * water, data = dft_araare) 
CH_model <- lm(log_Chroma ~ nutrient * water, data = dft_araare)
H_model <- lm(Hue ~ nutrient * water, data = dft_araare)
### Downstream ----
PV_model <- glmer.nb(Pollinator ~ nutrient + water + Plant_height_cm + log_Floral_area_mm + log_Chroma + (1 |Block), data = dfp_araare)

### SEM ----
araare_semmod <- piecewiseSEM::psem(
  PH_model,FN_model, FA_model, BCC_model, FCC_model, ACC_model, 
  TB_model, CH_model, H_model, PV_model, 
  Plant_height_cm %~~% sqrt_No_of_flowers,
  Plant_height_cm %~~% log_Floral_area_mm,
  sqrt_No_of_flowers %~~% log_Floral_area_mm,
  sqrt_No_of_flowers %~~% nutrient,
  sqrt_No_of_flowers %~~% Pollinator,
  Bee_chromatic_contrast %~~% Fly_chromatic_contrast,
  Bee_chromatic_contrast %~~% Achromatic_contrast,
  Bee_chromatic_contrast %~~% Total_brightness,
  Bee_chromatic_contrast %~~% log_Chroma,
  Bee_chromatic_contrast %~~% Hue,
  Bee_chromatic_contrast %~~% Pollinator,
  Fly_chromatic_contrast %~~% Hue,
  Fly_chromatic_contrast %~~% Pollinator,
  Achromatic_contrast %~~% Total_brightness,
  Achromatic_contrast %~~% log_Chroma,
  Achromatic_contrast %~~% Hue,
  Achromatic_contrast %~~% Pollinator,
  Achromatic_contrast %~~% water,
  Total_brightness %~~% log_Chroma,
  Total_brightness %~~% Hue,
  Total_brightness %~~% Pollinator,
  Hue %~~% Pollinator,
  data = dfp_araare
)
dSep(araare_semmod)
fisherC(araare_semmod)
piecewiseSEM::coefs(araare_semmod, standardize = "none")
araare_rsquared <- piecewiseSEM::rsquared(araare_semmod)
araare_coef <- piecewiseSEM::coefs(araare_semmod, standardize = "none")
write.csv(araare_rsquared, "araare_rsquared.csv", row.names = FALSE)
write.csv(araare_coef, "araare_sem-coefficients.csv", row.names = FALSE)


## DIA DEL SEM ----
### Upstream ----
PH_model <- lm(Plant_height_cm ~ nutrient + water, data = dft_diadel)
FN_model <- lm(sqrt_No_of_flowers ~ nutrient * water, data = dft_diadel) 
FA_model <- lm(log_Floral_area_mm ~ nutrient + water, data = dft_diadel)
BCC_model <- lm(Bee_chromatic_contrast ~ nutrient, data = dft_diadel)
FCC_model <- lm(Fly_chromatic_contrast ~ nutrient, data = dft_diadel)
ACC_model <- lm(Achromatic_contrast ~ nutrient, data = dft_diadel) 
TB_model<- lm(Total_brightness ~ nutrient * water, data = dft_diadel) 
CH_model <- lm(log_Chroma ~ nutrient, data = dft_diadel)
H_model <- lm(Hue ~ nutrient, data = dft_diadel)

### Downstream ----
PV_model <- glmer.nb(Pollinator ~ nutrient * water + Plant_height_cm + sqrt_No_of_flowers + log_Floral_area_mm + 
                            Bee_chromatic_contrast + Fly_chromatic_contrast + Achromatic_contrast + Total_brightness + log_Chroma + Hue + 
                            (1 |Block), data = dfp_diadel) # Not included in final model due to lack of significance

### SEM ----
diadel_semmod <- piecewiseSEM::psem(
  PH_model,FN_model, FA_model, BCC_model, FCC_model, 
  ACC_model, TB_model, CH_model, H_model,
  Plant_height_cm  %~~% Total_brightness,
  Bee_chromatic_contrast %~~% Fly_chromatic_contrast,
  Bee_chromatic_contrast %~~% Achromatic_contrast,
  Bee_chromatic_contrast %~~% Total_brightness,
  Bee_chromatic_contrast %~~% log_Chroma,
  Achromatic_contrast %~~% Total_brightness,
  Achromatic_contrast %~~% log_Chroma,
  Total_brightness %~~% log_Chroma,
  log_Chroma %~~% Hue,
  data = dft_diadel)

dSep(diadel_semmod)
fisherC(diadel_semmod)
piecewiseSEM::coefs(diadel_semmod, standardize = "none")
diadel_rsquared <- piecewiseSEM::rsquared(diadel_semmod)
diadel_coef <- piecewiseSEM::coefs(diadel_semmod, standardize = "none")
write.csv(diadel_rsquared, "diadel_rsquared.csv", row.names = FALSE)
write.csv(diadel_coef, "diadel_sem-coefficients.csv", row.names = FALSE)

## LEU IRC SEM ----
### Upstream ---- 
FN_model <- lm(sqrt_No_of_flowers ~ nutrient, data = dft_leuirc)
BCC_model <- lm(Bee_chromatic_contrast ~ nutrient * water, data = dft_leuirc)
FCC_model <- lm(Fly_chromatic_contrast ~ nutrient * water, data = dft_leuirc)
ACC_model <- lm(Achromatic_contrast ~ water, data = dft_leuirc)
TB_model<- lm(Total_brightness ~ water, data = dft_leuirc)
CH_model <- lm(log_Chroma ~ nutrient * water, data = dft_leuirc)
H_model <- lm(Hue ~ nutrient * water, data = dft_leuirc)

### Downstream ----
PV_model <- glmer.nb(Pollinator ~ water + log_Floral_area_mm + Bee_chromatic_contrast + Fly_chromatic_contrast + 
                                  Achromatic_contrast + Total_brightness + log_Chroma + Hue + (1 |Block), data = dfp_leuirc)

### SEM ----
leuirc_semmod <- piecewiseSEM::psem(
  FN_model, BCC_model, FCC_model, ACC_model, 
  TB_model, CH_model, H_model, PV_model, 
  Plant_height_cm  %~~% sqrt_No_of_flowers,
  Plant_height_cm %~~% log_Floral_area_mm, 
  sqrt_No_of_flowers %~~% Fly_chromatic_contrast,
  sqrt_No_of_flowers %~~% Hue,
  sqrt_No_of_flowers %~~% water,
  sqrt_No_of_flowers %~~% Pollinator,
  Bee_chromatic_contrast %~~% Fly_chromatic_contrast,
  Bee_chromatic_contrast %~~% Achromatic_contrast,
  Bee_chromatic_contrast %~~% Total_brightness,
  Bee_chromatic_contrast %~~% log_Chroma,
  Bee_chromatic_contrast %~~% Hue,
  Bee_chromatic_contrast %~~% Pollinator,
  Fly_chromatic_contrast %~~% Total_brightness,
  Fly_chromatic_contrast %~~% log_Chroma,
  Fly_chromatic_contrast %~~% Hue,
  Achromatic_contrast %~~% Total_brightness,
  Achromatic_contrast %~~% log_Chroma,
  Achromatic_contrast %~~% Hue,
  Total_brightness %~~% Hue,
  Pollinator %~~% nutrient,
  data = dfp_leuirc)

dSep(leuirc_semmod)
fisherC(leuirc_semmod)
piecewiseSEM::coefs(leuirc_semmod, standardize = "none")
leuirc_rsquared <- piecewiseSEM::rsquared(leuirc_semmod)
leuirc_coef <- piecewiseSEM::coefs(leuirc_semmod, standardize = "none")
write.csv(leuirc_rsquared, "leuirc_rsquared.csv", row.names = FALSE)
write.csv(leuirc_coef, "leuirc_sem-coefficients.csv", row.names = FALSE)

## POT ARG SEM ----
### Upstream ----
PH_model <- lm(Plant_height_cm ~ nutrient * water, data = dft_potarg)
FN_model <- lm(sqrt_No_of_flowers ~ nutrient, data = dft_potarg) 
FA_model <- lm(log_Floral_area_mm ~ nutrient * water, data = dft_potarg)
BCC_model <- lm(Bee_chromatic_contrast ~ nutrient, data = dft_potarg)
FCC_model <- lm(Fly_chromatic_contrast ~ nutrient, data = dft_potarg)
ACC_model <- lm(Achromatic_contrast ~ nutrient, data = dft_potarg) 
TB_model <- lm(Total_brightness ~ nutrient, data = dft_potarg) 
H_model <- lm(Hue ~ nutrient * water, data = dft_potarg)

### Downstream ----
PV_model <- glmer.nb(Pollinator ~ nutrient * water + Plant_height_cm + sqrt_No_of_flowers + log_Floral_area_mm + 
                                  Bee_chromatic_contrast + Fly_chromatic_contrast + Achromatic_contrast + 
                                  Total_brightness + log_Chroma + Hue + (1 |Block), data = dfp_potarg) # Not included in final model 

### SEM ----
potarg_semmod <- piecewiseSEM::psem(
  PH_model,FN_model, FA_model, BCC_model, FCC_model, 
  ACC_model, TB_model, H_model,
  Plant_height_cm  %~~% sqrt_No_of_flowers,
  Plant_height_cm %~~% log_Floral_area_mm, 
  sqrt_No_of_flowers %~~% log_Floral_area_mm,
  Bee_chromatic_contrast %~~% log_Chroma,
  Achromatic_contrast %~~% Total_brightness,
  data = dft_potarg)

dSep(potarg_semmod)
fisherC(potarg_semmod)
piecewiseSEM::coefs(potarg_semmod, standardize = "none")
potarg_rsquared <- piecewiseSEM::rsquared(potarg_semmod)
potarg_coef <- piecewiseSEM::coefs(potarg_semmod, standardize = "none")
write.csv(potarg_rsquared, "potarg_rsquared.csv", row.names = FALSE)
write.csv(potarg_coef, "potarg_sem-coefficients.csv", row.names = FALSE)

## VER CHA SEM ----
### Upstream ----
PH_model <- lm(Plant_height_cm ~ nutrient, data = dft_vercha)
FN_model <- lm(sqrt_No_of_flowers ~ nutrient, data = dft_vercha) 
FA_model <- lm(log_Floral_area_mm ~ nutrient, data = dft_vercha)
BCC_model <- lm(Bee_chromatic_contrast ~ nutrient + water, data = dft_vercha)
FCC_model <- lm(Fly_chromatic_contrast ~ nutrient * water, data = dft_vercha)
ACC_model <- lm(Achromatic_contrast ~ water, data = dft_vercha) 
TB_model<- lm(Total_brightness ~  water, data = dft_vercha) 
CH_model <- lm(log_Chroma ~ nutrient * water, data = dft_vercha)
H_model <- lm(Hue ~ water, data = dft_vercha)

### Downstream ----
PV_model <- glmer.nb(Pollinator ~ water + sqrt_No_of_flowers + (1 |Block), data = dfp_vercha)

### SEM ----
vercha_semmod <- piecewiseSEM::psem(
  PH_model,FN_model, FA_model, BCC_model, FCC_model, ACC_model, 
  TB_model, CH_model, PV_model,
  sqrt_No_of_flowers %~~% water,
  sqrt_No_of_flowers %~~% log_Floral_area_mm,
  Bee_chromatic_contrast %~~% Fly_chromatic_contrast,
  Bee_chromatic_contrast %~~% Achromatic_contrast,
  Bee_chromatic_contrast %~~% Total_brightness,
  Bee_chromatic_contrast %~~% log_Chroma,
  Fly_chromatic_contrast %~~% Achromatic_contrast,
  Fly_chromatic_contrast %~~% Total_brightness,
  Fly_chromatic_contrast %~~% log_Chroma,
  Achromatic_contrast %~~% Total_brightness,
  Achromatic_contrast %~~% log_Chroma,
  Total_brightness %~~% Hue,
  data = dfp_vercha)

dSep(vercha_semmod)
fisherC(vercha_semmod)
piecewiseSEM::coefs(vercha_semmod, standardize = "none")
vercha_rsquared <- piecewiseSEM::rsquared(vercha_semmod)
vercha_coef <- piecewiseSEM::coefs(vercha_semmod, standardize = "none")
write.csv(vercha_rsquared, "vercha_rsquared.csv", row.names = FALSE)
write.csv(vercha_coef, "vercha_sem-coefficients.csv", row.names = FALSE)
