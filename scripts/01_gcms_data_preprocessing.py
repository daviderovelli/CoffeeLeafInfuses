import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from src.utils import load_config # noqa: E402

#load config file
config = load_config(path='config/config.yaml', section='gc-ms')

FEATURES = Path(config['ftable'])
METADATA = Path(config['metadata'])

# === Metadata ========================================================
sample_col = 'ATTRIBUTE_Sample'
temp_col = 'ATTRIBUTE_Extraction'
group_col = 'ATTRIBUTE_Group'

# === Import data ========================================================
print("DATA IMPORT")

feat_df = pd.read_csv(FEATURES, index_col=0)
print(f"\n Feature table loaded: {feat_df.shape[0]} metabolites × {feat_df.shape[1]} samples")

meta_raw = pd.read_csv(METADATA, sep='\t')
print(f"Metadata loaded: {meta_raw.shape[0]} entries")


# === Extract Key Odorant Compounds 
ko_col = 'Key Odorant'
name_col = 'Compound'

valid_samples = [c for c in feat_df.columns if c in meta_raw[sample_col].values]

if ko_col in feat_df.columns:
    print(f"\n--- KEY ODORANT EXTRACTION ---")
    
    # 1. Filtra le righe con 'X'
    is_ko = feat_df[ko_col].astype(str).str.strip().str.upper() == 'X'
    
    if is_ko.sum() > 0:
        print(f" Found {is_ko.sum()} Key Odorant compounds.")
        
        ko_subset = feat_df.loc[is_ko].copy()
        if name_col in ko_subset.columns:
            ko_subset = ko_subset.set_index(name_col)
            print(f" Using '{name_col}' column as molecule names.")
        ko_clean = ko_subset[valid_samples]
        ko_export = ko_clean.T
        ko_export.index.name = 'Samples'
        output_ko_file = 'keyodorants_extracted.csv'
        ko_export.to_csv(output_ko_file)
        print(f" Extracted data saved to: {output_ko_file}")
        print(f" Format: {ko_export.shape[0]} samples (rows) x {ko_export.shape[1]} molecules (columns)")
    else:
        print(f" No compounds marked with 'X' found.")

    # 7. PULIZIA MAIN DATASET per le analisi successive
    # Manteniamo nel feat_df solo le colonne numeriche dei campioni
    # Rimuovendo colonne di testo come 'Key Odorant', 'Compound', 'RT' che causano errori
    print(f"\nCleaning feature table for processing...")
    feat_df = feat_df[valid_samples]
    print(f" Feature table cleaned: {feat_df.shape[0]} metabolites × {feat_df.shape[1]} samples")

else:
    print(f"\n WARNING: Column '{ko_col}' not found. Skipping extraction.")
    # Fallback pulizia: teniamo comunque solo i campioni validi
    feat_df = feat_df[valid_samples]

meta_raw = pd.read_csv(METADATA, sep='\t')
print(f"Metadata loaded: {meta_raw.shape[0]} entries")

# === Verify metadata-features matching =======================================
print("METADATA-FEATURES CONSISTENCY CHECK")

meta_raw = meta_raw[meta_raw[sample_col].isin(feat_df.columns)].copy()
samples_in_features = set(feat_df.columns)
samples_in_metadata = set(meta_raw[sample_col])
missing_in_metadata = samples_in_features - samples_in_metadata
missing_in_features = samples_in_metadata - samples_in_features

if missing_in_metadata:
    print(f"WARNING: {len(missing_in_metadata)} samples in features but NOT in metadata:")
    print(f"   {missing_in_metadata}")
else:
    print("All feature samples have metadata")

if missing_in_features:
    print(f"WARNING: {len(missing_in_features)} samples in metadata but NOT in features:")
    print(f"   {missing_in_features}")
else:
    print("All metadata samples are in features")

print(f"\nMetadata samples matched with features: {len(samples_in_metadata)}")

# === Split by extraction temperature ==========================================
print("Dataset split in HOT and COLD")

hot_cols = meta_raw.loc[meta_raw[temp_col].str.contains('hot', case=False, na=False), sample_col].tolist()
cold_cols = meta_raw.loc[meta_raw[temp_col].str.contains('cold', case=False, na=False), sample_col].tolist()

print(f"\nHOT samples: {len(hot_cols)}")
print(f"COLD samples: {len(cold_cols)}")

hot_df = feat_df[hot_cols].copy()
cold_df = feat_df[cold_cols].copy()

grp_map = meta_raw.set_index(sample_col)[group_col].rename('Group')

# === Check for NaN values ====================================================
print("NaN check")

def check_nan_values(df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    """
    Check for NaN values in the dataset.
    """
    nan_count = df.isna().sum().sum()
    nan_percentage = (nan_count / (df.shape[0] * df.shape[1])) * 100
    
    if nan_count > 0:
        print(f"\n {dataset_name}: {nan_count} NaN values ({nan_percentage:.2f}%)")
        nan_per_feature = df.isna().sum(axis=1)
        nan_features = nan_per_feature[nan_per_feature > 0]
        print(f"    Features with NaN: {len(nan_features)}")
        print(f"    Max NaN per feature: {nan_features.max()}")
        
        # Fill NaN with 0 (common for missing abundances)
        print("Filling NaN with 0")
        df = df.fillna(0)
    else:
        print(f"{dataset_name}: No NaN values detected")
    
    return df

hot_df = check_nan_values(hot_df, "HOT dataset")
cold_df = check_nan_values(cold_df, "COLD dataset")

# === FEATURE FILTERING ===============================================
print("Feature filtering: removing empty features")

def filter_empty_features(df: pd.DataFrame, 
                         min_prevalence: float = 0.2,
                         min_abundance: float = 0,
                         dataset_name: str = "Dataset") -> pd.DataFrame:
    """
    Filter features (rows) based on prevalence and abundance thresholds.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Feature matrix (features × samples)
    min_prevalence : float
        Minimum fraction of samples where feature must be detected.
        Default: 0.2 (feature must be in ≥20% of samples)
        Range: [0, 1]
    min_abundance : float
        Minimum value to consider a feature "detected"
        Default: 0 (any value > 0)
    dataset_name : str
        Name of dataset for logging
        
    """
    n_features_original = df.shape[0]
    n_samples = df.shape[1]
    
    # Calculate prevalence: fraction of samples where feature is detected
    detected = (df > min_abundance).sum(axis=1)
    prevalence = detected / n_samples
    
    # Apply filter
    mask = prevalence >= min_prevalence
    df_filtered = df[mask].copy()
    
    n_features_filtered = df_filtered.shape[0]
    n_removed = n_features_original - n_features_filtered
    
    print(f"\n{dataset_name} Filtering Results:")
    print(f"   Original features: {n_features_original}")
    print(f"   Removed: {n_removed} ({100*n_removed/n_features_original:.1f}%)")
    print(f"   Retained: {n_features_filtered} ({100*n_features_filtered/n_features_original:.1f}%)")
    print(f"   Samples: {n_samples}")
    print(f"   Prevalence threshold: {min_prevalence*100:.1f}% (feature in ≥{int(np.ceil(min_prevalence*n_samples))} samples)")
    
    return df_filtered

# Apply filters
hot_df = filter_empty_features(hot_df, min_prevalence=0.2, dataset_name="HOT")
cold_df = filter_empty_features(cold_df, min_prevalence=0.2, dataset_name="COLD")

# === Summary of filtering ====================================================
print("Feature filtering summary")

summary_data = {
    'Dataset': ['HOT', 'COLD'],
    'Features Retained': [hot_df.shape[0], cold_df.shape[0]],
    'Samples': [hot_df.shape[1], cold_df.shape[1]]
}

summary_df = pd.DataFrame(summary_data)
print("\n" + summary_df.to_string(index=False))

# === Data preprocessing ===================================================
print("Data autoscaling")

def process_data(df: pd.DataFrame, dataset_name: str = "Dataset") -> pd.DataFrame:
    """
    Apply autoscaling
    """
    print(f"\n{dataset_name}:")
    print(f"  Input shape: {df.shape[0]} metabolites × {df.shape[1]} samples")
    
    # Step 1: Transpose to (samples × metabolites)
    df_transposed = df.T
    print(f"  Transposed: {df_transposed.shape[0]} samples × {df_transposed.shape[1]} metabolites")
    
    # Step 2: Apply StandardScaler
    scaler = StandardScaler()
    scaled_array = scaler.fit_transform(df_transposed)
    
    # Step 3: Create DataFrame with proper index/columns
    scaled_df = pd.DataFrame(
        scaled_array,
        index=df.columns,      # Sample names (from original columns)
        columns=df.index       # Metabolite names (from original index)
    )
    
    # === VALIDATION: Per-metabolite statistics ===============================
    means_per_metabolite = scaled_df.mean(axis=0)
    stds_per_metabolite = scaled_df.std(axis=0)
    
    print("\n Autoscaling Validation (per metabolite):")
    print(f"     Means: min={means_per_metabolite.min():.2e}, max={means_per_metabolite.max():.2e} (expect ≈0)")
    print(f"     Stds:  min={stds_per_metabolite.min():.4f}, max={stds_per_metabolite.max():.4f} (expect ≈1)")
    
    # Strict validation
    if not np.allclose(means_per_metabolite, 0, atol=1e-10):
        print(" WARNING: Some metabolites have non-zero mean!")
    
    if not np.allclose(stds_per_metabolite, 1, atol=1e-10):
        print(" WARNING: Some metabolites have non-unit std!")
    
    # Global statistics
    global_mean = scaled_df.values.mean()
    global_std = scaled_df.values.std()
    print("\n  Global statistics:")
    print(f"     Mean: {global_mean:.6f}")
    print(f"     Std: {global_std:.6f}")
    print(f"     Output shape: {scaled_df.shape[0]} samples × {scaled_df.shape[1]} metabolites")
    
    return scaled_df

# Apply autoscaling
hot_scaled = process_data(hot_df, dataset_name="HOT Dataset")
cold_scaled = process_data(cold_df, dataset_name="COLD Dataset")

# === Final data summary =====================================================
print("Preprocessed data summary")

print("\n HOT Dataset (autoscaled):")
print(f"   Shape: {hot_scaled.shape[0]} samples × {hot_scaled.shape[1]} metabolites")
print(f"   Sample names: {list(hot_scaled.index[:5])}... (showing first 5)")
print(f"   Metabolite names: {list(hot_scaled.columns[:5])}... (showing first 5)")

print("\n COLD Dataset (autoscaled):")
print(f"   Shape: {cold_scaled.shape[0]} samples × {cold_scaled.shape[1]} metabolites")
print(f"   Sample names: {list(cold_scaled.index[:5])}... (showing first 5)")
print(f"   Metabolite names: {list(cold_scaled.columns[:5])}... (showing first 5)")

# === Save preprocessed data ========================================

hot_scaled.to_csv(config['ftable_hot'])
cold_scaled.to_csv(config['ftable_cold'])
print("File saved as feature_table_hot.csv and feature_table_cold.csv")

# === MetaboAnalyst export ===================================================
def export_for_metaboanalyst(df: pd.DataFrame, groups: pd.Series, output_path: str):
    """
    Export DataFrame for MetaboAnalyst 
    """
    # Create a copy to avoid modifying the original dataframe
    ma_df = df.copy()
    
    # Map the groups using the Sample ID (index)
    # We insert it at position 0 (first column)
    ma_df.insert(0, 'Class', ma_df.index.map(groups))
    
    # Give the index a name (MetaboAnalyst prefers 'Sample' or empty)
    ma_df.index.name = 'Sample'
    
    # Check for missing labels
    if ma_df['Class'].isna().any():
        n_missing = ma_df['Class'].isna().sum()
        print(f"WARNING: {n_missing} samples in {output_path} are missing a Group label!")
    
    # Save
    ma_df.to_csv(output_path)
    print(f" > Saved: {output_path} ({ma_df.shape[0]} samples x {ma_df.shape[1]} cols)")
    
    # Preview format
    print(f"   Format preview:\n{ma_df.iloc[:2, :3].to_string()}\n")

# Export HOT dataset
export_for_metaboanalyst(
    hot_scaled, 
    grp_map, 
    config['ma_ftable_hot']
)

# Export COLD dataset
export_for_metaboanalyst(
    cold_scaled, 
    grp_map, 
    config['ma_ftable_cold']
)

print("MetaboAnalyst export complete")
print("Volatile preprocessing script finished successfully.")