import pandas as pd
import numpy as np
import requests
import time
from concurrent.futures import ThreadPoolExecutor
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from src.utils import load_config # noqa: E402

#load config file
config = load_config(path='config/config.yaml', section='gc-ms')

# Path to the CSV file
print("Loading feature table...")
csv_path = config['ftable_input']
df = pd.read_csv(csv_path)

# Function to get SMILES from InChIKey via PubChem
print("Retrieving SMILES structures from PubChem...")
def get_smiles_from_inchikey(inchikey):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/CSV"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            lines = r.text.strip().split("\n")
            if len(lines) > 1:
                return lines[1].split(",")[1]
    except Exception:
        pass
    return None

# Add a new SMILES column if it doesn't already exist
if "SMILES" not in df.columns:
    smiles_list = []
    if "InChIKey" not in df.columns:
        raise ValueError("The 'InChIKey' column is not present in the CSV file.")
    for inchikey in df["InChIKey"]:
        if pd.isna(inchikey) or str(inchikey).strip() == "":
            smiles_list.append(None)
            continue
        smiles = get_smiles_from_inchikey(str(inchikey).strip())
        smiles_list.append(smiles)
        time.sleep(0.2)
    df["SMILES"] = [s.replace('"', '') if isinstance(s, str) else s for s in smiles_list]

# Cache for SMILES classification results
smiles_cache = {}

# Function to classify a SMILES structure
print("Classifying SMILES structures using NPClassifier...")
def classify_smiles(smiles):
    if pd.isna(smiles) or smiles is None:
        return None
    if smiles in smiles_cache:
        return smiles_cache[smiles]
    url = f"https://npclassifier.gnps2.org/classify?smiles={smiles}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            result = response.json()
            smiles_cache[smiles] = result
            return result
        else:
            print(f"Error API for {smiles}: {response.status_code}")
            return None
    except requests.RequestException as e:
        print(f"Connection error for {smiles}: {e}")
        return None

# Function to clean classified column values
def clean_value(value):
    """Normalize classification cell values into printable strings.
    """
    # Preserve None/NaN
    if value is None:
        return None
    try:
        if pd.isna(value):
            return value
    except Exception:
        # pd.isna may raise for some types; ignore and continue
        pass

    # Handle numpy arrays and pandas Series
    if isinstance(value, (np.ndarray, pd.Series)):
        if getattr(value, 'size', 0) == 0:
            return None
        # Convert to list for consistent handling
        value = value.tolist()

    # Join lists into comma-separated strings
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)

    # Clean simple strings
    if isinstance(value, str):
        return value.replace('[', '').replace(']', '').replace("'", '')

    # Fallback: stringify other types
    return str(value)

# Classify all SMILES structures using ThreadPoolExecutor
if "SMILES" not in df.columns:
    raise ValueError("The 'SMILES' column is not present in the DataFrame.")
with ThreadPoolExecutor(max_workers=5) as executor:
    classifications = list(executor.map(classify_smiles, df['SMILES']))
with ThreadPoolExecutor(max_workers=5) as executor:
    classifications = list(executor.map(classify_smiles, df['SMILES']))

# Add classification results to the DataFrame
df['class_results'] = [x.get('class_results') if x else None for x in classifications]
df['superclass_results'] = [x.get('superclass_results') if x else None for x in classifications]
df['pathway_results'] = [x.get('pathway_results') if x else None for x in classifications]
df['isglycoside'] = [x.get('isglycoside') if x else None for x in classifications]

# Clean classified columns
for column in ['class_results', 'superclass_results', 'pathway_results', 'isglycoside']:
    if column in df.columns:
        # Apply clean_value directly; clean_value internally handles NA/None and
        # array-like inputs safely.
        df[column] = df[column].apply(clean_value)

# Save the result
df.to_csv(config['ftable'], index=False)
print("File saved as feature_table.csv")