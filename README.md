# Sensomic profiling of coffee leaf infuses: unlocking the potential of a novel food

This repository contains all the scripts needed to reproduce the data analysis and results of the manuscript "Sensomic profiling of coffee leaf infuses: unlocking the potential of a novel food" (DOI here).

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
- The `cfg` folder contains `yaml` configuration files. File names and paths are listed in `config.yaml`. **Edit this file directly to change file names and paths.**
- The `data` folder contains [...] subfolders: 
    - [complete]   
- The `scripts` and `notebooks` folders contain Python scripts and Jupyter Notebooks to reproduce and run the data analysis (see [Usage](#usage)).
- The `results` folder contains output files generated from data analysis, including processed datasets, summary tables.

## Usage
The `paneldata_infuses_full.csv` contains the data collected during the sensory assesment of the samples under the analysis. The `feature_table.csv` contains all the annotations retrived in the study using a HS-SPME-GC-MS paltform (the method is described in the original publication).

Scripts: 
- `00_query_classification.py`: query PubChem to collect the SMILES for each annotation using the InChiKey values. Using this information, NPClassfier was used to retrive the chemical classes for each comopounds. Results are combined into the original `feature_table.csv` file.
- `01_data_preprocessing.py`: using the `feature_table.csv`

This script prepares GC–MS feature data for downstream multivariate analysis. It loads a features CSV and sample metadata TSV, ensures samples match, splits the dataset into hot vs cold extraction groups, imputes missing values (fills NaN with 0), filters out rare features (by prevalence), standardizes (autoscale) each metabolite across samples, and writes autoscaled hot/cold tables to outputs/


Notebooks: 
- `01_sensory.ipynb`: 

- `02_volatilome_stats.ipynb`:

## I/O files
