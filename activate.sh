#!/bin/bash

#check if conda environment CoffeeLeafInfuses already exists
if conda env list | grep -q 'CoffeeLeafInfuses'; then

    #if env exists, print message and activate environment
    echo "Conda environment 'CoffeeLeafInfuses' already exists. Activating existing environment..."
    conda activate CoffeeLeafInfuses

else
    #if env doesn't exist, create it and install packages in requirements.txt
    echo "Creating conda environment CoffeeLeafInfuses..."
    conda env create -f requirements.yaml
    conda activate CoffeeLeafInfuses
fi

#export cwd to PYTHONPATH
echo "Exporting current working directory to PYTHONPATH..."
export PYTHONPATH=$(pwd):$PYTHONPATH
echo "CoffeeLeafInfuses ready!"