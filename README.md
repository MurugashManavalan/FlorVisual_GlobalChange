# Interactive effects of Nutrient Enrichment and Water Scarcity on Floral signaling and visitation

## Manuscript details

### Authors

Murugash Manavalan, Dinesh Thakur, Bidhan Panthi, Zuzana Münzbergová

#### Corresponding Author

Murugash Manavalan, Department of Botany, Faculty of Science, Charles University, Prague, Czech Republic ; Institute of Botany, Czech Academy of Sciences, Průhonice, Czech Republic. Email: murugash.manavalan@gmail.com ; manavalm@natur.cuni.cz

### Summary 
- Floral traits remain an underexplored frontier in global change biology compared to their vegetative counterparts. As anthropogenic change involves simultaneous shifts in resource availability, focusing on isolated stressors may overlook complex biological responses. Characterizing the interactive effects of these stressors is therefore crucial to move beyond additive models and accurately predict the stability of pollination services in a transforming world. 
- The individual and interactive effects of nutrient enrichment and water scarcity were investigated on the floral traits and pollinator visitation rates of five grassland species.
- While species identity emerged as the dominant driver of overall trait variation, across-species responses were also present. Among the abiotic drivers, nutrient enrichment emerged as the primary influence, significantly increasing flower number and achromatic contrast while decreasing floral area, chromatic contrast and chroma. Conversely, water scarcity primarily reduced plant height and increased floral chroma and hue. The two factors also often interacted in their effects. The changes in visual traits scaled up to changes in pollinator visitation, but the specific effects were largely species-specific. 
- Ultimately, our findings demonstrate that anthropogenic stressors act in divergent ways to reshape floral phenotypes and pollinator visitation, advancing our understanding of plant-pollinator interactions under compounding global pressures.

## Project Details

| File                          | Description                                                                                                           |
|-------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| `Data/Loggers Data`           | Data from TOMST dataloggers allowing for validation of effects of experimental treatments on plant individuals        |
| `Data/Reflectance Data`       | Reflectance measurements of flower colour samples acquired using spectrometer                                         |
| `Data/Input`                  | Data files for pollinator visitation and traits                                                                       |
| `Data/Processed`              | Datasets for reflectance traits, pollinator visitation and traits (combined), traits correcting for effect of species |
| `R/functions.R`               | Functions required for the project alongside needed packages                                                          |
| `R/logger_analysis.R`         | Script for analysis of datalogger data                                                                                |
| `R/01_reflectance_analysis.R` | Script for processing of reflectance measurements                                                                     |
| `R/02_trait_analysis`         | Script for analysis of effect of abiotic factors on floral triats and pollinator visitation                           |
| `R/03_sem.R`                  | Script for Structural Equation Modelling (SEM)                                                                        |

### Getting the code

The project is accessible in two ways:
  
  1. If a user has a [GitHub account](https://github.com/), the easiest way is to [clone](https://happygitwithr.com/clone.html) this GitHub repo.
  
  2. A user can download the latest *Release* of the project as a zip file from the [Release page](https://github.com/MurugashManavalan/FlorVisual_GlobalChange/releases).

### Usage

To reproduce the analysis used in the publication, the user can run scripts in `R/` in suggested numerical order or run whichever scripts they are interested in.

## Citation

If you use the code or data, please cite the archived repository: 10.5281/zenodo.23883868