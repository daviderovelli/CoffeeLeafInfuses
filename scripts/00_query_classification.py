import pandas as pd
import requests
import time
from concurrent.futures import ThreadPoolExecutor

# Percorso al file CSV
csv_path = r"C:\Users\david\OneDrive - Università degli Studi di Parma\PhD\Projects\Coffee_leaf_infuses\data\processed\feature_table.csv"

# Leggi il CSV
df = pd.read_csv(csv_path)

# Funzione per ottenere SMILES da InChIKey tramite PubChem
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

# Appendi una nuova colonna SMILES se non esiste già
if "SMILES" not in df.columns:
    smiles_list = []
    if "InChIKey" not in df.columns:
        raise ValueError("La colonna 'InChIKey' non è presente nel file CSV.")
    for inchikey in df["InChIKey"]:
        if pd.isna(inchikey) or str(inchikey).strip() == "":
            smiles_list.append(None)
            continue
        smiles = get_smiles_from_inchikey(str(inchikey).strip())
        smiles_list.append(smiles)
        time.sleep(0.2)  # Rispetta i limiti di PubChem
    df["SMILES"] = [s.replace('"', '') if isinstance(s, str) else s for s in smiles_list]

# Cache per i risultati della classificazione SMILES
smiles_cache = {}

# Funzione per classificare una struttura SMILES
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
            print(f"Errore API per {smiles}: {response.status_code}")
            return None
    except requests.RequestException as e:
        print(f"Errore di connessione per {smiles}: {e}")
        return None

# Funzione per pulire i valori delle colonne classificate
def clean_value(value):
    if isinstance(value, list):
        # Unisci gli elementi in una stringa separata da virgole
        return ", ".join(str(v) for v in value)
    if isinstance(value, str):
        return value.replace('[', '').replace(']', '').replace("'", '')
# Classifica tutte le strutture SMILES usando ThreadPoolExecutor
if "SMILES" not in df.columns:
    raise ValueError("La colonna 'SMILES' non è presente nel DataFrame.")
with ThreadPoolExecutor(max_workers=5) as executor:
    classifications = list(executor.map(classify_smiles, df['SMILES']))
with ThreadPoolExecutor(max_workers=5) as executor:
    classifications = list(executor.map(classify_smiles, df['SMILES']))

# Aggiungi i risultati della classificazione al DataFrame
df['class_results'] = [x.get('class_results') if x else None for x in classifications]
df['superclass_results'] = [x.get('superclass_results') if x else None for x in classifications]
df['pathway_results'] = [x.get('pathway_results') if x else None for x in classifications]
df['isglycoside'] = [x.get('isglycoside') if x else None for x in classifications]

# Pulizia delle colonne di classificazione
for column in ['class_results', 'superclass_results', 'pathway_results', 'isglycoside']:
    if column in df.columns:
        df[column] = df[column].apply(lambda x: clean_value(x) if pd.notna(x) else x)

# Salva il risultato
df.to_csv(r"C:\Users\david\OneDrive - Università degli Studi di Parma\PhD\Projects\Coffee_leaf_infuses\data\processed\feature_table_with_smiles.csv", index=False)
print("File salvato come feature_table_with_smiles.csv")