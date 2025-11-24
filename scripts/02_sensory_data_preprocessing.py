import pandas as pd
from pathlib import Path
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from src.utils import load_config # noqa: E402

#load config file
config = load_config(path='config/config.yaml', section='sensory')
SENSORY = Path(config['sensory_table_input'])

# Data import
panel_df = pd.read_csv(SENSORY, index_col=0)

print("Original dataset:")
print(panel_df.head())
print(f"Shape: {panel_df.shape}")

# Trasformation from 0-5 to 1-6 scale
panel_df = panel_df + 1
print("\nAfter trasformation (0-5 → 1-6):")
print(panel_df.describe())

# Split hot and cold samples
sensory_table_hot = panel_df[panel_df.index.str.startswith("H-")]
sensory_table_cold = panel_df[panel_df.index.str.startswith("C-")]

print(f"\nHot samples: {len(sensory_table_hot)}")
print(sensory_table_hot.index.tolist())
print(f"\nCold samples: {len(sensory_table_cold)}")
print(sensory_table_cold.index.tolist())

# Save processed data
sensory_table_hot.to_csv(config['stable_hot'])
sensory_table_cold.to_csv(config['stable_cold'])
print("File saved as stable_hot.csv and stable_cold.csv")