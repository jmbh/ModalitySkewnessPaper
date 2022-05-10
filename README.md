# README
This is the reproducibility archive for the paper "Multimodality and Skewness in Emotion Time Series" (preprint: [https://psyarxiv.com/qudr6](https://psyarxiv.com/qudr6)). It allows one to reproduce all analyses and results in the paper.

### Folders
- `/Data` contains the preprocessed data of the seven studies we re-analyzed in our paper. The scripts used to preprocess the data based on the data files provided by the authors can be found at [http://github.com/jmbh/EmotionTimeSeries](http://github.com/jmbh/EmotionTimeSeries).

- `/Files` contains outputs created by the below scripts.

- `/Figures` contains the figures created by the below scripts.

### Scripts
- `aux_Functions.R` contains auxiliary functions for plotting, estimating modes with various methods we compare in the appendix, and ways to generate data which is used for the validation studies in the appendix.

- `EstimateModalitySkew.R` loads the seven datasets, estimates the modality, and computes the skewness of each individual distribution, and outputs the files with file names `EMS_...` (see below).

- `Evaluation_Figures.R` takes the files `A_DS_info.RDS`, `A_Modality.RDS`, `A_Densities.RDS` and `A_Skewness.RDS` (see below) as input, and plots all the Figures in the main text.

- `Illustration_Figure1.R` plots Figure 1 with the 8 example distributions of the emotion *sad*.

- `Additional_Analyses.R` contains the code to produce Appendix F (looking into the relationship between skewness and the location of the distribution) and Appendix G (looking into the absolute model fit of the VAR model).

- `Method_Validation.R` contains the code to reproduce the validation of our Mode-estimation method in Appendices A.2 and A.3.

- `Evaluation_ML_Skew.R` loads estimates of skewness, data-set characteristics and between-person variables. Performs all multilevel analyses on skewness described in the main text Section 2.3.2 and Appendix E (Tables 3 and 4, Figure 13).

- `Evaluation_ML_MModality.R` loads estimates of modality, data-set characteristics and between-person variables. Performs the study-specific multilevel analyses on modality described in the main text Section 2.3.1; prepares data for 3-level analysis done in Julia (see file `julia_analysis.jl`).

- `julia_analysis.jl` contains the Julia code to perform the 3-level multilevel model on modality in Section 2.3.1 and Appendix D.


### Files

- `EMS_DS_info.RDS` contains a list with various information on the dataset level (which items are included, which scale has been used, etc.).

- `EMS_Modality.RDS` contains an array with the modality estimates for all individual distributions; contains an `NA` if the distribution was excluded.

- `EMS_Densities.RDS` contains a list with the density estimates for all individual distributions, on which the modality estimates are based on; empty, if the distribution was excluded.

- `EMS_Skewness.RDS` contains an array with the skewness for each individual distribution.

- `EMS_Exclusion.RDS` contains a matrix indicating the counts for each exclusions separately for exclusion type and dataset.

- `EMS_TSlength.RDS` contains a list which contains matrices for each dataset, indicating the time series length for each subject.

- `EMS_TSlength_IDs.RDS` contains IDs for each individual in each study, including distributions that have been excluded.

- `SimResults_Validation.RDS` contain the results of the validation simulation in Appendix A.2.

- `SimResults_Validation_Skew.RDS` contain the results of the validation simulation in Appendix A.3.

- `SimResults_Validation_Skew.RDS` contain the results of the validation simulation in Appendix A.3.

- `ModalityML_Data.csv` dataset created by `Evaluation_ML_MModality.R` to run the 3-level multilevel model on modality in the file `julia_analysis.jl`.

- `BetweenData.RDS` contains between-person data for all datasets. Obtained from: [http://github.com/jmbh/EmotionTimeSeries](http://github.com/jmbh/EmotionTimeSeries).
