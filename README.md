## Comprehensive prediction of synthetic lethality between paralog pairs in cancer cell lines


### Data analysis notebooks overview:
| Notebook                               | Figures                 | Brief description                                        |
|:---------------------------------------|:------------------------|:---------------------------------------------------------|
| 01_identify_paralog_SLs_w_depmap | Fig. 1                  | Derive bronze standard dataset of (non-) SL paralog pairs           |
| 02_indiv_feature_analysis        | Fig. 2A, S1             | Investigate the indv. predictive power of 22 features of paralog pairs |
| 03_train_classifier              | Fig. 2B                 | Train a random forest classifier w/ bronze standard data |
| 04_make_predictions              | Fig. 3, S2, S4          | Make predictions for ind. combinatorial screens + all paralog pairs |
| 05_explain_predictions           | Fig. 4, 5A, 5B, 5D, S5  | Generate SHAP profiles for indv. predictions |
| 06_compare_slant                 | Fig. S3                 | Compare performance of RF predictions to that of SLant predictions |
| 07_experimental_validation       | Fig. 5C, 5E             | Graph siRNA screen results for *ASF1A/ASF1B* and *COPS7A/COPS7B* |
| 08_human_vs_cerevisiae_paralogs  | NA                      | Compare the avg. family size for paralogs in humans vs. cerevisiae yeast |

### Data processing notebooks overview:
| Notebook                               | Brief description                                        |
|:---------------------------------------|:---------------------------------------------------------|
| 01_process_hgnc_id_map                 | Create map with HGNC symbol, Entrez ID, Ensembl ID and locus type (e.g. protein-coding) for each gene |
| 02_process_ensembl_paralogs            | Generate list of protein-coding paralog pairs with min. 20% sequence identity + some feat. annotations |
| 03_process_copy_number                 | Identify homozygous deletions from copy number and gene expression data |
| paralog_features/01_process_protein_complexes | Compute protein complex membership and essentiality features for paralog pairs |
| paralog_features/02_process_ppi        | Compute protein-protein interaction features for paralog pairs  |
| paralog_features/03_process_orthologs  | Identify (essential) orthologs for paralog pairs and compute conservation score |
| paralog_features/04_process_gtex_expression   | Compute expression features for paralog pairs |
| paralog_features/05_merge_all_features | Generate full list of all paralog pairs annotated with all features - some features, incl. age, are computed in this notebook |
| GI_screens/process_thompson_pairs      | Process gene pairs screened by Thompson et al., filter to overlap with our paralog pairs list  |
| GI_screens/process_dede_pairs          | Process gene pairs screened by Dede et al., filter to overlap with our paralog pairs list  |



### To run Jupyter notebooks:
* Obtain 3rd party data - sources are listed in data_sources
    * Note: this is primarily needed for running notebooks in the 1_data_processing folder
* Set path to 3rd party data directory in environment.yml (3rd_party_dir)
  * Note: this can also be set in the activate environment with: `conda env config vars set 3rd_party_data=my_dir`
* Create [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environment from the environment.yml file: `conda env create -f environment.yml`
* Activate conda environment: `conda activate paralogSL`
* Start notebooks: `jupyter notebook`




