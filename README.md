# Sensomic profiling of coffee leaf infusions: chemical–sensory markers for quality assessment

This repository contains the scripts required to reproduce the data analysis and results of the manuscript: **Sensomic profiling of coffee leaf infusions: chemical–sensory markers for quality assessment ([DOI Link])**

**Data Availability and Confidentiality** 

Due to **confidentiality agreements regarding proprietary corporate data**, the original datasets cannot be made publicly available in this repository.

- **Original data**: The full CSV files are deposited on Zenodo under restricted access [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18266459.svg)](https://doi.org/10.5281/zenodo.18266459). Researchers wishing to access the data for validation purposes may request it via Zenodo, subject to approval. 
- **Mock data**: We have provided synthetic `mock` datasets within this repository. These files mimic the structure and format of the original data, allowing users to test the pipeline and verify the functionality of the scripts.

## Requirements

- [Miniconda/Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- Python 3.11.0

## Installation and setup
1. Clone this GitHub repository by running the following command in your terminal:

~~~
git clone https://github.com/daviderovelli/CoffeeLeafInfuses
~~~

2. Create a new conda environment and install the required packages and dependencies from `requirements.yaml` by running:

~~~
conda env create -f requirements.yaml
conda activate coffeeleafinfuses
~~~

Alternatively, if you are using macOS or Linux, you can run the setup script instead:

~~~
source activate.sh
~~~

This script checks if the environment already exists, creates it if needed, installs dependencies, and exports the repository path to `PYTHONPATH`. 

## Repository structure
- The `config` folder contains `yaml` configuration files. File names and paths are listed in `config.yaml`. Edit this file directly to change file names and paths.

- The `data` folder contains two subfolders: 

    - `input/mock`: contains the input tables from the sensory assessment (`paneldata_infuses_full.csv`) and the volatile feature table (`input_feature_table.csv`).

    - `ouput/mock`: contains the outupts generated from the python scripts (see [Usage/Scipts](#scripts)) used as prerequisite for the whole analysis carried outusing the Jupiter Notebooks (see [Usage/Notebooks](#notebooks)).

- The `scripts` and `notebooks` folders contain Python scripts and Jupyter Notebooks to reproduce and run the data analysis (see [Usage](#usage)).

- The `results` folder contains output files generated from data analysis such as summary tables and figures.

## Usage
The `paneldata_infuses_full.csv` contains the data collected during the sensory assesment of the samples under the analysis. The `input_feature_table.csv` contains all the annotations retrived in the study using a HS-SPME-GC-MS paltform (the method is described in the original publication).

### Scripts: 
- `00_query_classification.py`: query the InChiKey in the feature table and fetch the cannonical SMILES form PubChem. For every SMILES, it classifies the compounds using NPClassifier API.

- `01_gcms_data_preprocessing.py`: the script loads a feature table and a metadata table, checks consistency, splits samples by temperature (Hot and Cold), filters and cleans feature and autoscale the data. In addition, it creates CSV files formatted for MetaboAnalyst. 

- `02_sensory_data_preprocessing.py`: the script loads a sensory table and a metadata table. It transforms scores from a 0–5 scale to a 1–6 scale and separates samples into hot and cold groups. The processed subsets are then saved as individual CSV files for downstream analysis.

### Notebooks: 
- `01_sensory.ipynb`: performs multivariate analysis of processed sensory data, including PCA and MCA. It loads hot and cold sensory datasets, explores variance structures, and visualizes sample relationships. Custom plotting functions generate score plots and radar charts for clear comparative interpretation on the selected samples. 

- `02_volatilome_stats.ipynb`: performs multivariate analysis (PCA) on volatile profiles of hot and cold infusions. Validates group separation and variance homogeneity using PERMANOVA and PERMDISP. It generates hierarchical-clustering heatmaps to visualize metabolite abundance patterns and sample correlations. An in-depth analysis for the keyodorants compounds is also provided. 
