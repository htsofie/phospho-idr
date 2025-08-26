# %% [markdown]
    # # AF3 Output Analysis

# %%
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np

# %%
# Analyze plDDT scores for all prediction folders

def extract_atom_residue_mapping(cif_file):
    """Extract atom-to-residue mapping from CIF file"""
    atom_mapping = []
    
    with open(cif_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                parts = line.split()
                if len(parts) >= 17:
                    atom_id = int(parts[1])
                    residue_name = parts[5]
                    chain_id = parts[6]
                    residue_number = int(parts[8])  # Column 9 (index 8)
                    atom_mapping.append((atom_id, residue_name, chain_id, residue_number))
    
    return atom_mapping

# Define prediction folders and their corresponding data files
prediction_folders = {}

# Q13625 Ankyrin Repeats
prediction_folders['Q13625_Ankyrin_Repeats'] = {
    'folder': 'P53_Q13625_ankyrin_repeats_prediction',
    'json_file': 'fold_p04637_full_q13625_ankyrin_repeats_full_data_0.json',
    'cif_file': 'fold_p04637_full_q13625_ankyrin_repeats_model_0.cif',
    'summary_json_file': 'fold_p04637_full_q13625_ankyrin_repeats_summary_confidences_0.json',
    'color': 'blue'
}

# Q8WUF5 Ankyrin Repeats
prediction_folders['Q8WUF5_Ankyrin_Repeats'] = {
    'folder': 'P53_Q8WUF5_Ankyrin_repeats_prediction',
    'json_file': 'fold_p04637_full_q8wuf5_ankyrin_repeats_full_data_0.json',
    'cif_file': 'fold_p04637_full_q8wuf5_ankyrin_repeats_model_0.cif',
    'summary_json_file': 'fold_p04637_full_q8wuf5_ankyrin_repeats_summary_confidences_0.json',
    'color': 'green'
}

# Q8WUF5 SH3 Domain
prediction_folders['Q8WUF5_SH3_Domain'] = {
    'folder': 'P53_Q8WUF5_Variant_SH3_domain_predictions',
    'json_file': 'fold_p04637_full_q8wuf5_variant_sh3_domain_full_data_0.json',
    'cif_file': 'fold_p04637_full_q8wuf5_variant_sh3_domain_model_0.cif',
    'summary_json_file': 'fold_p04637_full_q8wuf5_variant_sh3_domain_summary_confidences_0.json',
    'color': 'red'
}

# Rep FactorA NTD
prediction_folders['Rep_FactorA_NTD'] = {
    'folder': 'P53_Rep_FactorA_NTD_Prediction',
    'json_file': 'fold_p04637_full_p27694_replication_factor_a_protein_1_n_terminal_domain_full_data_0.json',
    'cif_file': 'fold_p04637_full_p27694_replication_factor_a_protein_1_n_terminal_domain_model_0.cif',
    'summary_json_file': 'fold_p04637_full_p27694_replication_factor_a_protein_1_n_terminal_domain_summary_confidences_0.json',
    'color': 'orange'
}

# SH3 Domain
prediction_folders['SH3_Domain'] = {
    'folder': 'P53_SH3_domain_prediction',
    'json_file': 'fold_p04637_full_q13625_sh3_domain_full_data_0.json',
    'cif_file': 'fold_p04637_full_q13625_sh3_domain_model_0.cif',
    'summary_json_file': 'fold_p04637_full_q13625_sh3_domain_summary_confidences_0.json',
    'color': 'purple'
}

# TFIIH P62 NTD
prediction_folders['TFIIH_P62_NTD'] = {
    'folder': 'P53_TFIIH_P62_NTD_Prediction',
    'json_file': 'fold_p04637_full_p32780_tfiih_p62_subunit_n_terminal_domain_full_data_0.json',
    'cif_file': 'fold_p04637_full_p32780_tfiih_p62_subunit_n_terminal_domain_model_0.cif',
    'summary_json_file': 'fold_p04637_full_p32780_tfiih_p62_subunit_n_terminal_domain_summary_confidences_0.json',
    'color': 'brown'
}

# Dictionary to store results for all predictions
all_results = {}

# Analyze each prediction folder
for prediction_name, config in prediction_folders.items():
    try:
        # Load JSON data
        json_path = f"sequence_files/{config['folder']}/{config['json_file']}"
        with open(json_path, 'r') as file:
            data = json.load(file)
        
        # Extract atom-to-residue mapping
        cif_path = f"sequence_files/{config['folder']}/{config['cif_file']}"
        atom_mapping = extract_atom_residue_mapping(cif_path)
        
        # Calculate chain A and B averages and residue-level plDDT scores
        chain_a_scores = []
        chain_b_scores = []
        
        # Create a dictionary to group atoms by residue for chain A
        residue_atoms_A = {}
        
        for i, (atom_id, residue_name, chain_id, residue_number) in enumerate(atom_mapping):
            if chain_id == 'A':
                chain_a_scores.append(data['atom_plddts'][i])
                # Group atoms by residue for chain A
                if residue_number not in residue_atoms_A:
                    residue_atoms_A[residue_number] = []
                residue_atoms_A[residue_number].append(data['atom_plddts'][i])
            elif chain_id == 'B':
                chain_b_scores.append(data['atom_plddts'][i])
        
        # Calculate average plDDT for each residue in chain A
        residue_numbers_A = sorted(residue_atoms_A.keys())
        residue_plddt_A = [np.mean(residue_atoms_A[res_num]) for res_num in residue_numbers_A]
        
        # Store results
        all_results[prediction_name] = {
            'chain_a_mean': np.mean(chain_a_scores) if chain_a_scores else 0,
            'chain_b_mean': np.mean(chain_b_scores) if chain_b_scores else 0,
            'residue_numbers': residue_numbers_A,
            'residue_plddt': residue_plddt_A,
            'color': config['color']
        }
        
    except FileNotFoundError as e:
        print(f"Warning: Could not find files for {prediction_name}")
    except Exception as e:
        print(f"Error analyzing {prediction_name}: {e}")

# %%
# View summary values for each prediction
print(f"{'Prediction':<30} {'pTM':<10} {'ipTM':<10} {'PAE_min':<10} {'Clash':<8} {'Disordered':<12}")

for prediction_name, config in prediction_folders.items():
    try:
        summary_json_path = f"sequence_files/{config['folder']}/{config['summary_json_file']}"
        with open(summary_json_path, 'r') as file:
            summary_data = json.load(file)
        
        # Extract key metrics (with fallbacks if keys don't exist)
        ptm = summary_data.get('ptm', 'N/A')
        iptm = summary_data.get('iptm', 'N/A')
        # PAE is stored as chain_pair_pae_min, take the minimum value
        pae_data = summary_data.get('chain_pair_pae_min', [])
        if pae_data and isinstance(pae_data, list) and len(pae_data) > 0:
            # Flatten the nested list and find the minimum PAE value
            flat_pae = [item for sublist in pae_data for item in sublist]
            pae_min = min(flat_pae) if flat_pae else 'N/A'
        else:
            pae_min = 'N/A'
        
        clash = summary_data.get('has_clash', 'N/A')
        disordered = summary_data.get('fraction_disordered', 'N/A')
        
        print(f"{prediction_name:<30} {ptm:<10} {iptm:<10} {pae_min:<10} {clash:<8} {disordered:<12}")
        
    except Exception as e:
        print(f"{prediction_name:<30} {'Error':<10} {'Error':<10} {'Error':<10} {'Error':<8} {'Error':<12}")
        print(f"  Error details: {e}")

# %% [markdown]
# Summary Table of pLDDT scores for each chain

print(f"{'Prediction':<30} {'Chain A Mean':<15} {'Chain B Mean':<15}")

for prediction_name, results in all_results.items():
    print(f"{prediction_name:<30} {results['chain_a_mean']:<15.3f} {results['chain_b_mean']:<15.3f}")

# %% [markdown]
# Plot of pLDDT scores for p53 across all residues

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), 
                                gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.1})

# Main plot (top)
for prediction_name, results in all_results.items():
    if 'residue_numbers' in results and 'residue_plddt' in results:
        ax1.plot(results['residue_numbers'], results['residue_plddt'], 
                 color=results['color'], linewidth=1.5, alpha=0.8, 
                 label=prediction_name)

ax1.set_ylabel('plDDT Score')
ax1.set_title('Combined plDDT Scores for Chain A (p53) Across All Domain Predictions')
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_xticks([])  # Remove x-axis labels from main plot

# Set consistent x-axis limits for both plots
ax1.set_xlim(1, 393)

# Domain diagram (bottom)
# Define p53 domains with residue ranges
domains = [
    {'name': 'Transactivation', 'start': 1, 'end': 42, 'color': 'purple', 'alpha': 0.7},
    {'name': 'Proline-rich', 'start': 63, 'end': 97, 'color': 'orange', 'alpha': 0.7},
    {'name': 'DNA-binding', 'start': 98, 'end': 292, 'color': 'lightgreen', 'alpha': 0.7},
    {'name': 'Tetramerization', 'start': 324, 'end': 355, 'color': 'lightblue', 'alpha': 0.7},
    {'name': 'C-terminal', 'start': 363, 'end': 393, 'color': 'steelblue', 'alpha': 0.7}
]

# Plot domain rectangles
for domain in domains:
    ax2.add_patch(plt.Rectangle((domain['start'], 0), 
                                domain['end'] - domain['start'], 1, 
                                facecolor=domain['color'], alpha=domain['alpha'],
                                edgecolor='black', linewidth=0.5))
    
    # Add domain name
    center_x = (domain['start'] + domain['end']) / 2
    ax2.text(center_x, 0.5, domain['name'], ha='center', va='center', 
             fontsize=8, fontweight='bold', rotation=45)

# Set up domain plot
ax2.set_xlim(1, 393)
ax2.set_ylim(0, 1)
ax2.set_xlabel('Residue Number (Chain A - p53)')
ax2.set_yticks([])
ax2.set_ylabel('Domains')
ax2.set_title('p53 Protein Domains', fontsize=10, pad=10)

# Add residue number ticks at domain boundaries
tick_positions = [1, 42, 63, 97, 98, 292, 324, 355, 363, 393]
tick_labels = ['1', '42', '63', '97', '98', '292', '324', '355', '363', '393']
ax2.set_xticks(tick_positions)
ax2.set_xticklabels(tick_labels, rotation=45, ha='right')

# Ensure both plots have identical x-axis setup
ax1.set_xlim(1, 393)
ax2.set_xlim(1, 393)

plt.show()

# %%
