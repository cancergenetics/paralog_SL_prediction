## Comprehensive prediction of synthetic lethality between paralog pairs in cancer cell lines


#### Data analysis notebooks overview
| Notebook                                | Figures                 | Brief description                                        |
| --------------------------------------- |:------------------------| ---------------------------------------------------------|
| 01_identify_paralog_SLs_w_depmap.ipynb  | Fig. 1                  | Derive bronze standard dataset of (non)-SL paralog pairs           |
| 02_indiv_feature_analysis.ipynb         | Fig. 2A, S1             | Look at the indv. predictive power of 22 features of paralog pairs |
| 03_train_classifier.ipynb               | Fig. 2B                 | Train a random forest classifier w/ bronze standard data |
| 04_make_predictions.ipynb               | Fig. 3, S2, S4          | Make predictions for ind. combinatorial screens + all paralog pairs |
| 05_explain_predictions.ipynb            | Fig. 4, 5A, 5B, 5D, S5  | Generate SHAP profiles for indv. predictions |
| 06_compare_slant.ipynb                  | Fig. S3                 | Compare performance of RF predictions to that of SLant predictions |
| 07_experimental_validation.ipynb        | Fig. 5C, 5E             | Graph siRNA screen results for ASF1A/B and COPS7A/B |
| 08_human_vs_cerevisiae_paralogs.ipynb   | NA                      | Compare the avg. family size for paralogs in humans vs. cerevisiae yeast |



#### To run Jupyter notebooks:
* Obtain 3rd party data - sources are listed in data_sources
    * Note: this is primarily needed for running notebooks in the 1_data_processing folder
* Set path to 3rd party data directory in environment.yml (3rd_party_dir)
  * Note: this can also be set in the activate environment with: `conda env config vars set 3rd_party_data=my_dir`
* Create [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environment from the environment.yml file: `conda env create -f environment.yml`
* Activate conda environment: `conda activate paralogSL`
* Start notebooks: `jupyter notebook`




