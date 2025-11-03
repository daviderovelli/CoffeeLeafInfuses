# === Imports ==============================================================
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler


# === Path ================================================================
FEATURES = Path(r"C:\Users\david\OneDrive - Università degli Studi di Parma\PhD\Projects\Coffee_leaf_infuses\data\processed\feature_table.csv")
METADATA = Path(r"C:\Users\david\OneDrive - Università degli Studi di Parma\PhD\Projects\Coffee_leaf_infuses\data\gcms\volatile_infuses_metadata.tsv")

output_dir = Path(".") / "outputs"
output_dir.mkdir(exist_ok=True)

# === Metadata ========================================================
sample_col = 'ATTRIBUTE_Sample'
temp_col = 'ATTRIBUTE_Extraction'
group_col = 'ATTRIBUTE_Group'

# === Import data ========================================================
print("\n" + "="*70)
print("DATA IMPORT")
print("="*70)

feat_df = pd.read_csv(FEATURES, index_col=0)
print(f"\n Feature table loaded: {feat_df.shape[0]} metabolites × {feat_df.shape[1]} samples")

meta_raw = pd.read_csv(METADATA, sep='\t')
print(f"Metadata loaded: {meta_raw.shape[0]} entries")

# === Verify metadata-features matching =======================================
print("\n" + "="*70)
print("METADATA-FEATURES CONSISTENCY CHECK")
print("="*70)

meta_raw = meta_raw[meta_raw[sample_col].isin(feat_df.columns)].copy()
samples_in_features = set(feat_df.columns)
samples_in_metadata = set(meta_raw[sample_col])
missing_in_metadata = samples_in_features - samples_in_metadata
missing_in_features = samples_in_metadata - samples_in_features

if missing_in_metadata:
    print(f"WARNING: {len(missing_in_metadata)} samples in features but NOT in metadata:")
    print(f"   {missing_in_metadata}")
else:
    print(f"All feature samples have metadata")

if missing_in_features:
    print(f"WARNING: {len(missing_in_features)} samples in metadata but NOT in features:")
    print(f"   {missing_in_features}")
else:
    print(f"All metadata samples are in features")

print(f"\nMetadata samples matched with features: {len(samples_in_metadata)}")

# === Split by extraction temperature ==========================================
print("\n" + "="*70)
print("DATASET SPLIT: HOT vs COLD")
print("="*70)

hot_cols = meta_raw.loc[meta_raw[temp_col].str.contains('hot', case=False, na=False), sample_col].tolist()
cold_cols = meta_raw.loc[meta_raw[temp_col].str.contains('cold', case=False, na=False), sample_col].tolist()

print(f"\nHOT samples: {len(hot_cols)}")
print(f"COLD samples: {len(cold_cols)}")

hot_df = feat_df[hot_cols].copy()
cold_df = feat_df[cold_cols].copy()

grp_map = meta_raw.set_index(sample_col)[group_col].rename('Group')

# === Check for NaN values ====================================================
print("\n" + "="*70)
print("NaN VALUE CHECK")
print("="*70)

def check_nan_values(df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    """
    Check for NaN values in the dataset.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Feature matrix to check
    dataset_name : str
        Name of the dataset (for printing)
        
    Returns:
    --------
    pd.DataFrame
        Input dataframe unchanged (for chaining)
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
        print(f"    → Filling NaN with 0...")
        df = df.fillna(0)
    else:
        print(f"✓ {dataset_name}: No NaN values detected")
    
    return df

hot_df = check_nan_values(hot_df, "HOT dataset")
cold_df = check_nan_values(cold_df, "COLD dataset")

# === FEATURE FILTERING ===============================================
print("\n" + "="*70)
print("FEATURE FILTERING: Removing empty features")
print("="*70)

def filter_empty_features(df: pd.DataFrame, 
                         min_prevalence: float = 0.2,
                         min_abundance: float = 0,
                         dataset_name: str = "Dataset") -> pd.DataFrame:
    """
    Filter features (rows) based on prevalence and abundance thresholds.
    
    A feature is retained if it appears (with value > min_abundance) in at least
    min_prevalence fraction of samples. This is standard in metabolomics to remove
    rare or missing features that may introduce noise.
    
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
        
    Returns:
    --------
    pd.DataFrame
        Filtered feature matrix (only features meeting thresholds)
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
    
    print(f"\n📊 {dataset_name} Filtering Results:")
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
print("\n" + "="*70)
print("FEATURE FILTERING SUMMARY")
print("="*70)

summary_data = {
    'Dataset': ['HOT', 'COLD'],
    'Features Retained': [hot_df.shape[0], cold_df.shape[0]],
    'Samples': [hot_df.shape[1], cold_df.shape[1]]
}

summary_df = pd.DataFrame(summary_data)
print("\n" + summary_df.to_string(index=False))

# === Data preprocessing ===================================================
print("\n" + "="*70)
print("DATA AUTOSCALING (STANDARDIZATION)")
print("="*70)

def process_data(df: pd.DataFrame, dataset_name: str = "Dataset") -> pd.DataFrame:
    """
    Apply autoscaling (standardization) to the feature data.
    
    Autoscaling ensures each metabolite (feature) has mean=0 and std=1 across all
    samples. This is crucial for multivariate analysis methods like PCA and PLS
    that are sensitive to variable scaling.
    
    Process:
    --------
    1. Transpose matrix from (metabolites × samples) to (samples × metabolites)
    2. Apply StandardScaler to each column (metabolite) independently
    3. Each metabolite gets: (value - mean) / std
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input feature matrix (metabolites × samples)
    dataset_name : str
        Name of dataset for logging
        
    Returns:
    --------
    pd.DataFrame
        Autoscaled feature matrix (samples × metabolites)
        Each metabolite: mean≈0, std≈1
    """
    print(f"\n{dataset_name}:")
    print(f"  Input shape: {df.shape[0]} metabolites × {df.shape[1]} samples")
    
    # Step 1: Transpose to (samples × metabolites)
    df_transposed = df.T
    print(f"  Transposed: {df_transposed.shape[0]} samples × {df_transposed.shape[1]} metabolites")
    
    # Step 2: Apply StandardScaler
    # StandardScaler().fit_transform() scales COLUMNS independently
    # So each metabolite (column) gets mean=0, std=1 across samples
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
    
    print(f"\nAutoscaling Validation (per metabolite):")
    print(f"     Means: min={means_per_metabolite.min():.2e}, max={means_per_metabolite.max():.2e} (expect ≈0)")
    print(f"     Stds:  min={stds_per_metabolite.min():.4f}, max={stds_per_metabolite.max():.4f} (expect ≈1)")
    
    # Strict validation
    if not np.allclose(means_per_metabolite, 0, atol=1e-10):
        print(f"WARNING: Some metabolites have non-zero mean!")
    
    if not np.allclose(stds_per_metabolite, 1, atol=1e-10):
        print(f"WARNING: Some metabolites have non-unit std!")
    
    # Global statistics
    global_mean = scaled_df.values.mean()
    global_std = scaled_df.values.std()
    print(f"\n  Global statistics:")
    print(f"     Mean: {global_mean:.6f}")
    print(f"     Std: {global_std:.6f}")
    print(f"     Output shape: {scaled_df.shape[0]} samples × {scaled_df.shape[1]} metabolites")
    
    return scaled_df

# Apply autoscaling
hot_scaled = process_data(hot_df, dataset_name="HOT Dataset")
cold_scaled = process_data(cold_df, dataset_name="COLD Dataset")

# === Final data summary =====================================================
print("\n" + "="*70)
print("PREPROCESSING COMPLETE - FINAL DATA SUMMARY")
print("="*70)

print("\n HOT Dataset (autoscaled):")
print(f"   Shape: {hot_scaled.shape[0]} samples × {hot_scaled.shape[1]} metabolites")
print(f"   Sample names: {list(hot_scaled.index[:5])}... (showing first 5)")
print(f"   Metabolite names: {list(hot_scaled.columns[:5])}... (showing first 5)")

print("\n COLD Dataset (autoscaled):")
print(f"   Shape: {cold_scaled.shape[0]} samples × {cold_scaled.shape[1]} metabolites")
print(f"   Sample names: {list(cold_scaled.index[:5])}... (showing first 5)")
print(f"   Metabolite names: {list(cold_scaled.columns[:5])}... (showing first 5)")

# === Optional: Save preprocessed data ========================================
print("\n" + "="*70)
print("OPTIONAL: SAVING PREPROCESSED DATA")
print("="*70)

output_dir = Path(".") / "outputs"
output_dir.mkdir(exist_ok=True)

try:
    hot_scaled.to_csv(output_dir / "hot_scaled.csv")
    cold_scaled.to_csv(output_dir / "cold_scaled.csv")
    print(f"\n✓ Preprocessed data saved to {output_dir}/")
    print(f"  - hot_scaled.csv")
    print(f"  - cold_scaled.csv")
except Exception as e:
    print(f"\n Error saving files: {e}")

# === Display sample data ====================================================
print("\n" + "="*70)
print("SAMPLE DATA")
print("="*70)

print("\nHOT Dataset:")
print(hot_scaled.iloc[:5, :5])

print("\nCOLD Dataset:")
print(cold_scaled.iloc[:5, :5])

print("\n" + "="*70)
print("PREPROCESSING COMPLETE")
print("="*70 + "\n")

