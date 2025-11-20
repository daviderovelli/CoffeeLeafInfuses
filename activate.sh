#!/bin/bash

#check if conda environment CoffeeLeafInfuses already exists
if conda env list | grep -q 'coffeeleafinfuses'; then

    #if env exists, print message and activate environment
    echo "Conda environment 'coffeeleafinfuses' already exists. Activating existing environment..."
    conda activate coffeeleafinfuses

else
    #if env doesn't exist, create it and install packages in requirements.txt
    echo "Creating conda environment coffeeleafinfuses..."
    conda env create -f requirements.yaml
    conda activate coffeeleafinfuses
fi

#export cwd to PYTHONPATH
echo "Exporting current working directory to PYTHONPATH..."
export PYTHONPATH=$(pwd):$PYTHONPATH
echo "CoffeeLeafInfuses ready!"