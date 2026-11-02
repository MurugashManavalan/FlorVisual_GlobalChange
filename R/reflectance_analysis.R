# Reflectance Analysis of all files (Using PAVO) ----
## Load required packages 
install.packages("pavo")
library(pavo)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(here)

## Importing spectra
specs <- getspec(where = here::here("Data", "Reflectance Data"), ext = "TXT", lim = c(300,700),sep = ";")

## Average repeated measurements
spectral_names <- names(specs)[names(specs) != "wl"]
species <- sub("_.*", "", spectral_names)
treatment <- sub("^[^_]+_([A-Z]{2}).*$", "\\1", spectral_names)
individual <- sub("^.*_([0-9]+)\\..*$", "\\1", spectral_names)
obs <- sub("^.*\\.([0-9]+)$", "\\1", spectral_names)

mspec <- aggspec(specs, by = individual, FUN = mean)
sppspec <- aggspec(specs, by = species, FUN = mean)

## Spectra Smoothing
plotsmooth(mspec, minsmooth = 0.05, maxsmooth = 0.5, curves = 4, ask = FALSE)
specs.sm <- procspec(specs, opt = "smooth", span = 0.2)

## Converting Negative Values to Zero
specs.smn <- procspec(specs.sm, fixneg = "zero")

## Creation of Colour Space
vis.bee <- vismodel(specs.smn, visual = "apis", qcatch = "Ei", relative = FALSE, vonkries = TRUE, achromatic = "l", illum = "D65", bkg = "green") # For bees
vis.fly <- vismodel(specs.smn, visual = "musca", qcatch = "Qi", relative = TRUE, vonkries = TRUE, illum = "D65", bkg = "green") # For flies

## Coordinates for Colour Space
hex.flowers <- colspace(vis.bee, space = "hexagon")
cat.flowers <- colspace(vis.fly, space = "categorical")

## Calculating Chromatic Contrast
hex_contrast <- sqrt(hex.flowers$x^2 + hex.flowers$y^2) #Euclidean Distance between colour loci of flower and achromatic centre (0,0)
cat_contrast <- sqrt(cat.flowers$x^2 + cat.flowers$y^2)

## Calculating Achromatic Contrast
vis.bee_green <- vis.bee$l
achro_contrast <- abs(vis.bee_green - 0.5)
View(specs.smn)

## Dataframes for Contast
bee_contrast <- data.frame(
  flower = colnames(specs.smn)[colnames(specs.smn) != "wl"],
  chromatic_contrast = hex_contrast,
  achromatic_contrast = achro_contrast
)

fly_contrast <- data.frame(
  flower = colnames(specs.smn)[colnames(specs.smn) != "wl"],
  chromatic_contrast = cat_contrast
)

## Other colourimetric variables
spec_summary <- summary(specs.smn)
B1 <- spec_summary$B1
S8 <- spec_summary$S8
H1 <- spec_summary$H1

## Save as CSV File
final_df <- data.frame(
  species = species,
  treatment = treatment,
  individual = individual,
  observation = obs,
  bee_chromatic_contrast = hex_contrast,
  fly_chromatic_contrast = cat_contrast,
  achromatic_contrast = achro_contrast,
  total_brightness = B1,
  chroma = S8,
  hue = H1
)
write_csv(final_df, here("Data", "reflectance_results.csv"))

