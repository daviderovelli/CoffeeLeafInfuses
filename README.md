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

Before, running any scripts, set the `PYTHONPATH` to the project root so that internal imports work correctly:
~~~
$env:PYTHONPATH = "$(Get-Location)"
~~~

Alternatively, if you are using macOS or Linux, you can run the setup script instead:

~~~
source activate.sh
~~~

This script checks if the environment already exists, creates it if needed, installs dependencies, and exports the repository path to `PYTHONPATH`. 

## Repository structure
- The `config` folder contains `yaml` configuration files. File names and paths are listed in `config.yaml`. Edit this file directly to change file names and paths.

- The `data` folder contains two subfolders: 

    - `input`: contains the input tables from the sensory assessment (`paneldata_infuses_full.csv`) and the volatile feature table (`input_feature_table.csv`).

    - `processed`: contains the outupt generated from the python scripts (see [Usage/Scipts](#scripts)) used as prerequisite for the whole analysis. 

- The `scripts` and `notebooks` folders contain Python scripts and Jupyter Notebooks to reproduce and run the data analysis (see [Usage](#usage)).

- The `results` folder contains output files generated from data analysis, summary tables and figures.

## Usage
The `paneldata_infuses_full.csv` contains the data collected during the sensory assesment of the samples under the analysis. The `feature_table.csv` contains all the annotations retrived in the study using a HS-SPME-GC-MS paltform (the method is described in the original publication).

### Scripts: 
- `00_query_classification.py`: query the InChiKey in the feature table and fetch the cannonical SMILES form PubChem. For every SMILES, it classifies the compounds using NPClassifier API.

- `01_gcms_data_preprocessing.py`: the script loads a feature table and a metadata table, checks consistency, splits samples by extraction temperature (HOT vs COLD), filters and cleans features, autoscale-standardizes the data per metabolite. 

- `02_sensory_data_preprocessing.py`: the script loads a sensory table and a metadata table [TODO]

### Notebooks: 
- `01_sensory.ipynb`: 

- `02_volatilome_stats.ipynb`:

## I/O files
