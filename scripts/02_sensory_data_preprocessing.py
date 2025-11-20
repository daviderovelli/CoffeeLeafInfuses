import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import prince
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from src.utils import load_config # noqa: E402

# === Data import ===
panel_path = r"c:\Users\david\OneDrive - Università degli Studi di Parma\PhD\Projects\Coffee_leaf_infuses\data\processed\paneldata_infuses_full.csv"
panel_df = pd.read_csv(panel_path, index_col=0)

print("Original dataset:")
print(panel_df.head())
print(f"Shape: {panel_df.shape}")

# === Trasformation 0-5 to 1-6 ===
panel_df = panel_df + 1
print("\nAfter trasformation (0-5 → 1-6):")
print(panel_df.describe())

# === Split hot and cold infusions ===
hot_panel = panel_df[panel_df.index.str.startswith("H-")]
cold_panel = panel_df[panel_df.index.str.startswith("C-")]

print(f"\nHot samples: {len(hot_panel)}")
print(hot_panel.index.tolist())
print(f"\nCold samples: {len(cold_panel)}")
print(cold_panel.index.tolist())

# === Function to calculate variable contributions ===
def calculate_variable_contributions(mca, df_cat, n_components=2):
    col_coords = mca.column_coordinates(df_cat)
    variable_contributions = {}
    
    for var in df_cat.columns:
        var_modalities = [col for col in col_coords.index if col.startswith(f"{var}_")]
        for comp in range(n_components):
            if comp not in variable_contributions:
                variable_contributions[comp] = {}
            contrib = sum([col_coords.loc[mod, comp]**2 for mod in var_modalities if mod in col_coords.index])
            variable_contributions[comp][var] = contrib
    
    return variable_contributions

# === Aggregate coordinates of variables ===
def get_variable_coordinates(mca, df_cat):
    """Calculate the aggregate coordinates of variables as the mean of the coordinates of their modalities"""
    col_coords = mca.column_coordinates(df_cat)
    variable_coords = {}
    for var in df_cat.columns:
        var_modalities = [col for col in col_coords.index if col.startswith(f"{var}_")]
        if var_modalities:
            coords = col_coords.loc[var_modalities].mean(axis=0)
            variable_coords[var] = coords
    return pd.DataFrame(variable_coords).T

# === MCA (all samples) ===
def run_mca(df, title, save_prefix):
    print(f"\n{'='*50}")
    print(f"MCA analysis: {title}")
    print(f"{'='*50}")
    
    df_cat = df.astype(str)
    mca = prince.MCA(n_components=5, n_iter=10, random_state=42).fit(df_cat)

    # Explained variance
    eigenvalues = mca.eigenvalues_
    explained_variance = (eigenvalues / eigenvalues.sum() * 100)
    print("Explained variance per component:")
    for i, var in enumerate(explained_variance[:5]):
        print(f"  Component {i+1}: {var:.2f}%")
    print(f"Total first 2 components: {explained_variance[:2].sum():.2f}%")

    row_coords = mca.row_coordinates(df_cat)
    var_coords = get_variable_coordinates(mca, df_cat)
    
    # Plot
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Sample points
    sns.scatterplot(x=row_coords.iloc[:, 0], y=row_coords.iloc[:, 1],
                   s=100, edgecolor="black", linewidth=0.8, alpha=0.9,
                   ax=ax, color="tab:blue", label="Sample")
    for i, sample in enumerate(df.index):
        ax.text(row_coords.iloc[i, 0], row_coords.iloc[i, 1], sample,
               fontsize=9, ha='right', va='bottom', weight='bold')

    # Variable arrows
    for var, (x, y) in var_coords.iloc[:, [0, 1]].iterrows():
        ax.arrow(0, 0, x*0.8, y*0.8, color='black', alpha=0.7,
                 head_width=0.05, length_includes_head=True, linewidth=1.5)
        ax.text(x*0.9, y*0.9, var, color='black', fontsize=9,
                ha='center', va='center', weight='bold')
    
    ax.set_title(f"{title} – MCA Biplot", fontsize=14, weight='bold')
    ax.set_xlabel(f"Comp 1 ({explained_variance[0]:.1f}%)", fontsize=12)
    ax.set_ylabel(f"Comp 2 ({explained_variance[1]:.1f}%)", fontsize=12)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.axhline(0, color="grey", lw=1)
    ax.axvline(0, color="grey", lw=1)
    
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:blue',
               markeredgecolor='black', markersize=8, label='Sample'),
        Line2D([0], [0], color='black', lw=2, label='Variable')
    ]
    ax.legend(handles=handles, loc='best', fontsize=10, frameon=True)
    plt.tight_layout()
    
    os.makedirs("outputs/MCA", exist_ok=True)
    plt.savefig(f"outputs/MCA/{save_prefix}_panel_mca_biplot.jpg", dpi=600, bbox_inches='tight')
    plt.savefig(f"outputs/MCA/{save_prefix}_panel_mca_biplot.svg", bbox_inches='tight')
    plt.show()
    
    return mca, df_cat

# === MCA (on selected samples) ===
selected_samples = [
    "C-MA","C-M-EC","C-M-T5","C-N-GE-WS","C-WI",
    "H-MA","H-M-EC","H-M-T5","H-N-GE-WS","H-WI"
]

def run_mca_selected(df, sample_prefix, title, save_prefix):
    print(f"\n{'='*50}")
    print(f"MCA analysis: {title}")
    print(f"{'='*50}")
    
    filtered = df.loc[[s for s in df.index if s in selected_samples and s.startswith(sample_prefix)]]
    print(f"Sample points analyzed: {filtered.index.tolist()}")
    
    df_cat = filtered.astype(str)
    mca = prince.MCA(n_components=5, n_iter=10, random_state=42).fit(df_cat)
    
    eigenvalues = mca.eigenvalues_
    explained_variance = (eigenvalues / eigenvalues.sum() * 100)
    print("Explained variance per component:")
    for i, var in enumerate(explained_variance[:5]):
        print(f"  Component {i+1}: {var:.2f}%")
    print(f"Total first 2 components: {explained_variance[:2].sum():.2f}%")

    row_coords = mca.row_coordinates(df_cat)
    var_coords = get_variable_coordinates(mca, df_cat)
    
    # Plot
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(10, 8))
    
    color = "red" if sample_prefix == "H" else "blue"
    sns.scatterplot(x=row_coords.iloc[:, 0], y=row_coords.iloc[:, 1],
                   s=120, edgecolor="black", linewidth=0.8, alpha=0.9,
                   ax=ax, color=color, label="Sample")
    for i, sample in enumerate(filtered.index):
        ax.text(row_coords.iloc[i, 0], row_coords.iloc[i, 1], sample,
               fontsize=10, ha='right', va='bottom', weight='bold')
    
    for var, (x, y) in var_coords.iloc[:, [0, 1]].iterrows():
        ax.arrow(0, 0, x*0.8, y*0.8, color='darkgreen', alpha=0.7,
                 head_width=0.05, length_includes_head=True, linewidth=1.5)
        ax.text(x*0.9, y*0.9, var, color='darkgreen', fontsize=9,
                ha='center', va='center', weight='bold')
    
    ax.set_title(f"{title} – MCA Biplot", fontsize=14, weight='bold')
    ax.set_xlabel(f"Comp 1 ({explained_variance[0]:.1f}%)", fontsize=12)
    ax.set_ylabel(f"Comp 2 ({explained_variance[1]:.1f}%)", fontsize=12)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.axhline(0, color="grey", lw=1)
    ax.axvline(0, color="grey", lw=1)
    
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
               markeredgecolor='black', markersize=10, label='Campione'),
        Line2D([0], [0], color='darkgreen', lw=2, label='Variabile')
    ]
    ax.legend(handles=handles, loc='best', fontsize=10, frameon=True)
    plt.tight_layout()
    
    plt.savefig(f"outputs/MCA/{save_prefix}_panel_mca_biplot.jpg", dpi=600, bbox_inches='tight')
    plt.savefig(f"outputs/MCA/{save_prefix}_panel_mca_biplot.svg", bbox_inches='tight')
    plt.show()
    
    # Variable contributions
    print(f"\n{'='*40}")
    print("VARIABLE CONTRIBUTIONS")
    print(f"{'='*40}")
    
    var_contributions = calculate_variable_contributions(mca, df_cat)
    
    comp1_contrib = sorted(var_contributions[0].items(), key=lambda x: x[1], reverse=True)
    comp2_contrib = sorted(var_contributions[1].items(), key=lambda x: x[1], reverse=True)

    print("\nTop 10 variabiles Component 1:")
    for var, contrib in comp1_contrib[:10]:
        print(f"  {var}: {contrib:.4f}")

    print("\nTop 10 variabiles Component 2:")
    for var, contrib in comp2_contrib[:10]:
        print(f"  {var}: {contrib:.4f}")
    
    total_contrib = {var: var_contributions[0][var] + var_contributions[1][var] for var in var_contributions[0]}
    total_sorted = sorted(total_contrib.items(), key=lambda x: x[1], reverse=True)
    
    print("\nTop 10 variabiles Total contribution (Comp1 + Comp2):")
    for var, contrib in total_sorted[:10]:
        print(f"  {var}: {contrib:.4f}")

    # === Esporta dataset CSV con campioni e variabili top per ciascuna componente ===
    def export_top_variables_csv(filtered_df, top_vars, filename):
        df_export = filtered_df[top_vars].copy()
        df_export.insert(0, "Sample", df_export.index)
        df_export.to_csv(f"outputs/MCA/{filename}", index=False)
        print(f"Salvato: outputs/MCA/{filename}")

    # Top 10 variabili per ciascuna componente
    top10_comp1 = [var for var, _ in comp1_contrib[:10]]
    top10_comp2 = [var for var, _ in comp2_contrib[:10]]

    # Esporta per componente 1
    export_top_variables_csv(filtered, top10_comp1, f"{save_prefix}_top_comp1.csv")
    # Esporta per componente 2
    export_top_variables_csv(filtered, top10_comp2, f"{save_prefix}_top_comp2.csv")
    
    # Graphs of contributions
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    vars1, contribs1 = zip(*comp1_contrib[:10])
    axes[0].barh(range(len(vars1)), contribs1, color='skyblue')
    axes[0].set_yticks(range(len(vars1)))
    axes[0].set_yticklabels(vars1)
    axes[0].set_title(f'{title}\nContribution - Comp 1')
    axes[0].set_xlabel('Contribution')
    axes[0].invert_yaxis()
    
    vars2, contribs2 = zip(*comp2_contrib[:10])
    axes[1].barh(range(len(vars2)), contribs2, color='lightcoral')
    axes[1].set_yticks(range(len(vars2)))
    axes[1].set_yticklabels(vars2)
    axes[1].set_title(f'{title}\nContribution - Comp 2')
    axes[1].set_xlabel('Contribution')
    axes[1].invert_yaxis()
    
    vars_tot, contribs_tot = zip(*total_sorted[:10])
    axes[2].barh(range(len(vars_tot)), contribs_tot, color='lightgreen')
    axes[2].set_yticks(range(len(vars_tot)))
    axes[2].set_yticklabels(vars_tot)
    axes[2].set_title(f'{title}\nTotal Contribution (Comp1 + Comp2)')
    axes[2].set_xlabel('Contribution')
    axes[2].invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(f"outputs/MCA/{save_prefix}_variable_contributions.jpg", dpi=600, bbox_inches='tight')
    plt.savefig(f"outputs/MCA/{save_prefix}_variable_contributions.svg", bbox_inches='tight')
    plt.show()
    
    return mca, df_cat

# === Execute ===
print("Creazione cartella outputs...")
os.makedirs("outputs/MCA", exist_ok=True)

mca_hot_all, _ = run_mca(hot_panel, "Hot infuses", "hot_all")
mca_cold_all, _ = run_mca(cold_panel, "Cold infuses", "cold_all")

mca_hot_sel, _ = run_mca_selected(panel_df, "H", "Hot infuses selected", "hot_selected")
mca_cold_sel, _ = run_mca_selected(panel_df, "C", "Cold infuses selected", "cold_selected")