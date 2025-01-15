# phylogeography_efp_repository
Scripts and data for Science submission Williams et al Jan 2025
(C) 2025 Geoffrey M Williams

This project was funded by USDA

## Citation

Williams GM, Augustinus B, Brockerhoff EG, Christen S, Coyle DR, Gougherty A, Lynch SC, Sakalidis ML, Roy K, Stallman JK, Franic I, Ghelardini L, Santini A, Sniezko RA. 2025. Phylogeography drives forest pathogen establishment. In review.

## Description of contents of respository

### Covariates_input
The main folder contains derived data used as covariates for imports, introduced trees and shrubs, and climate. Also contains intermediate derived data produced in the process.

#### Map_KG-Global
From a previous publication: F. Rubel, K. Brugger, K. Haslinger, I. Auer, The climate of the European Alps: Shift of very high resolution Köppen-Geiger climate zones 1800–2100. Meteorologische Zeitschrift 26, 115-125 (2017)
<https://koeppen-geiger.vu-wien.ac.at>

#### New_Trade_Data
Raw data from WITS <https://wits.worldbank.org> and US Census <https://usatrade.census.gov> used to calculate trade data.

### Figures_Sept2024 & Final_figures
Figure output

### Final_output
Contains the data analyzed in this study, along with some of the results

#### Data analyzed or used in study
1. Final_analysis_data_pathogens-Sept-2024.csv
2. pest_list_taxonomy.xlsx and .csv

#### Results of GLM and ANCOVA
1. linearmodels.RData
2. results.ancova.offset.stats.Sept_2024_2.csv
3. results.glm.centered.stats.Sept_2024_pois_2.csv
4. zinf_model_final.RData

#### Results of partial least square analysis
1. final_pls_R2.csv
2. PLS-models-9-10-23.RData
3. linearmodels.RData

#### Results of structural equation modeling
1. SEM_grouped_results.csv

### Host_list
Contains derived data on country tree lists, exotic tree lists, global tree phylogenies, distance matrices, and ordinations (in folder "Output") and raw data used to derive them (other folders)
1. Other_DB: includes data collected from BGCI and Tropicos (the latter used in initial approaches to separate Hawaii from BGCI entry for US)
2. Invasive_Trees_States and USDA_DB: data downloaded from USDA Plants database
3. Input: data from the GloNAF database

### Pathogens
Data on pathogen diversity, both raw data (Input) and derived (Output). "Expanded_data_EU_AU_NA_HI_Sept2024.csv" is the dataset of origins and invaded ranges of pathogens analyzed.

### Scripts
Scripts used in the analyes. Within individual folders, any order in which scripts should be run to produce results are denoted by labels, letter or number, prefixed to the file names. The input and output of the folders are in the repository above. In most cases, scripts are meant to be run from the parent directory.
1. preliminary_a_Global_Tree_Phylogeny: scripts used to make a plant name synonym table (folder Synonym_table) and global tree phylogeny for calculation of unifrac distances among countries (folder Tree_of_trees)
2. preliminary_b_Covariate_Scripts: scripts used to derive climate distances (derive_covariate_data.R), import volume (trade_data.R), and lists of established non native forest trees and shrubs.
3. c_Analysis_Scripts and c_Figures_Scripts: scripts used to put the final data for analysis together and run  analyses, and produce figures and tables of main results
4. d_Supplemental_Scripts: scripts used to produce supplemental material

## NOTE: Files not in the repository due to size limitations:
1. Figures_Sept2-24/climate.maps2.RData
2. Figures_Sept2024/Fig3_composite_ABCDEF.RData
3. Final_figures/Fig3_composite_ABCDEF.emf
4. Host_list/Output/host_translation_table_all_new_wfo.csv * (available upon request)