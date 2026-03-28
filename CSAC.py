import os
import pandas as pd
import requests
from Bio.PDB import PDBParser, Selection, NeighborSearch, is_aa
import numpy as np

# Prompt user for input file path
input_file = input("Please enter the path to your CSAC_input.csv file: ")

# Read the existing PDB ID and metadata with header
pdb_list_df = pd.read_csv(input_file)

# Prompt user for metals to search for
metals_input = input(
    "Enter metals to search for (2-letter element codes in capitals, separated by commas. Example: FE, NI, MN): "
)
metals = [metal.strip().upper() for metal in metals_input.split(',') if metal.strip()]

# Analysis unit is fixed to ligand-cluster centroid for v2 consistency.
analysis_unit = "ligand_centroid"

# Initialize an empty DataFrame to store the analysis data
additional_columns = [
    'Analysis Unit',
    'Center Definition',
    'Center X',
    'Center Y',
    'Center Z',
    'Target Metal Count in Ligand',
    'Chain',
    'Metal',
    'ligand type',
    'ligand radius',
    'Nearby Atoms',
    'Atom Types',
    'Average Electronegativity',
    'Electronegativity SD',
    'Final Nearby Residues',
    'Average Hydropathy',
    'Hydropathy SD',
    'Hydropathy Scores',
]
columns = list(pdb_list_df.columns) + additional_columns
dataframe = pd.DataFrame(columns=columns)

# Keep track of failed downloads and missing metals
failed_downloads = []
missing_metals = []


# Function to calculate the center of a cluster
def calculate_center(atoms):
    coords = np.array([atom.get_coord() for atom in atoms])
    center = np.mean(coords, axis=0)
    return center


# Function to calculate ligand radius from center
def calculate_ligand_radius(atoms, center):
    if not atoms:
        return 0.0
    distances = [np.linalg.norm(atom.get_coord() - center) for atom in atoms]
    return float(np.max(distances)) if distances else 0.0


# Dictionary to map three-letter codes to one-letter codes
aa3to1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Kyte-Doolittle hydropathy scale
hydropathy = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Electronegativity values for selected elements
electronegativity = {
    'FE': 1.83, 'S': 2.58, 'ZN': 1.65, 'P': 2.19, 'H': 2.20,
    'O': 3.44, 'MO': 2.16, 'C': 2.55, 'NI': 1.91, 'N': 3.04
}


# Function to find nearby residues given a model, center point, and distance threshold
def find_any_nearby_residues(model, center, threshold):
    atoms = Selection.unfold_entities(model, 'A')
    ns = NeighborSearch(atoms)
    nearby_atoms = ns.search(center, threshold, 'A')
    residues = set()
    for atom in nearby_atoms:
        residue = atom.get_parent()
        if is_aa(residue):
            residues.add(residue)
    return residues


# Function to find nearby atoms given a model, center point, and distance threshold
def find_nearby_atoms(model, center, threshold=4):
    atoms = Selection.unfold_entities(model, 'A')
    ns = NeighborSearch(atoms)
    return ns.search(center, threshold, 'A')


# Process one analysis unit and append row
def process_center(
    dataframe,
    model,
    pdb_id,
    metadata,
    chain_id,
    residue,
    target_metal,
    center,
    ligand_radius,
    center_definition,
    target_metal_count,
    central_atom=None,
):
    nearby_atoms = find_nearby_atoms(model, center)

    # Central atom exclusion is only relevant to per-atom mode (not used in v2).
    if central_atom is not None:
        nearby_atoms = [a for a in nearby_atoms if a != central_atom]

    nearby_atoms_identifiers = [
        f"{a.get_parent().get_resname()}{a.get_serial_number()}" for a in nearby_atoms
    ]
    nearby_atom_types = [a.element for a in nearby_atoms]

    electroneg_vals = [
        electronegativity[elem] for elem in nearby_atom_types if elem in electronegativity
    ]
    avg_electronegativity = float(np.mean(electroneg_vals)) if electroneg_vals else np.nan
    electronegativity_sd = float(np.std(electroneg_vals)) if electroneg_vals else np.nan

    nearby_residues_5a = find_any_nearby_residues(model, center, 5.0)
    nearby_residues_5a_one_letter = [
        aa3to1[res.resname] for res in nearby_residues_5a if res.resname in aa3to1
    ]

    nearby_residues_outer_shell = find_any_nearby_residues(model, center, ligand_radius + 2.0)
    nearby_residues_outer_shell_one_letter = [
        aa3to1[res.resname] for res in nearby_residues_outer_shell if res.resname in aa3to1
    ]

    # Ion uses 5 A sphere; ligand uses radius+2 A sphere.
    use_ion_shell = np.isclose(ligand_radius, 0.0, atol=1e-8)
    final_nearby_residues = (
        nearby_residues_5a_one_letter if use_ion_shell else nearby_residues_outer_shell_one_letter
    )

    scores = [hydropathy[res] for res in final_nearby_residues]
    avg_score = float(np.mean(scores)) if scores else np.nan
    std_dev = float(np.std(scores)) if scores else np.nan

    new_row = {
        'PDB ID': pdb_id,
        'Analysis Unit': analysis_unit,
        'Center Definition': center_definition,
        'Center X': float(center[0]),
        'Center Y': float(center[1]),
        'Center Z': float(center[2]),
        'Target Metal Count in Ligand': int(target_metal_count),
        'Chain': chain_id,
        'Metal': target_metal,
        'ligand type': residue.resname,
        'ligand radius': float(ligand_radius),
        'Nearby Atoms': ', '.join(nearby_atoms_identifiers),
        'Atom Types': ', '.join(nearby_atom_types),
        'Average Electronegativity': avg_electronegativity,
        'Electronegativity SD': electronegativity_sd,
        'Final Nearby Residues': ', '.join(final_nearby_residues),
        'Average Hydropathy': avg_score,
        'Hydropathy SD': std_dev,
        'Hydropathy Scores': scores,
    }

    new_row.update(metadata)
    return pd.concat([dataframe, pd.DataFrame([new_row])], ignore_index=True)


# Iterate over each PDB ID
for _, row in pdb_list_df.iterrows():
    pdb_id = str(row['PDB ID']).strip()
    print(pdb_id)

    metadata = {col: row[col] for col in pdb_list_df.columns if col != 'PDB ID'}
    for key, value in metadata.items():
        print(f"{key}: {value}")

    found_metals = set()

    pdb_file = f"{pdb_id}.pdb"
    if not os.path.exists(pdb_file):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                with open(pdb_file, "w") as file:
                    file.write(response.text)
            else:
                print(f"Failed to download the file for PDB ID: {pdb_id}")
                failed_downloads.append(pdb_id)
                continue
        except requests.RequestException:
            print(f"Failed to download the file for PDB ID: {pdb_id}")
            failed_downloads.append(pdb_id)
            continue

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)

    # Use only the first model per structure to avoid ensemble-level pseudoreplication
    # (especially relevant for NMR entries with multiple MODEL blocks).
    for model_idx, model in enumerate(structure):
        if model_idx > 0:
            break
        for chain in model:
            for residue in chain:
                residue_atoms = list(residue.get_atoms())
                target_atoms = [a for a in residue_atoms if a.element in metals]
                if not target_atoms:
                    continue

                center = calculate_center(residue_atoms)
                ligand_radius = calculate_ligand_radius(residue_atoms, center)
                metals_in_residue = sorted({a.element for a in target_atoms})
                for target_metal in metals_in_residue:
                    found_metals.add(target_metal)
                    target_metal_count = sum(1 for a in residue_atoms if a.element == target_metal)
                    dataframe = process_center(
                        dataframe=dataframe,
                        model=model,
                        pdb_id=pdb_id,
                        metadata=metadata,
                        chain_id=chain.id,
                        residue=residue,
                        target_metal=target_metal,
                        center=center,
                        ligand_radius=ligand_radius,
                        center_definition='ligand_cluster_centroid',
                        target_metal_count=target_metal_count,
                        central_atom=None,
                    )

    missing = set(metals) - found_metals
    if missing:
        for metal in missing:
            missing_entry = {'PDB ID': pdb_id, 'Missing Metal': metal}
            missing_entry.update(metadata)
            missing_metals.append(missing_entry)

# Save the dataframe to a CSV file
input_filename = os.path.splitext(os.path.basename(input_file))[0]
output_filename = input_filename + '_' + '_'.join(metals) + '_DATABASE.csv'
dataframe.to_csv(output_filename, index=False)

# Print and save failed downloads
if failed_downloads:
    print("\nFailed downloads:")
    for pdb_id in failed_downloads:
        print(pdb_id)
    with open('failed_downloads.txt', 'w') as f:
        f.write('\n'.join(failed_downloads))

# Print structures missing metals
if missing_metals:
    print("\nStructures missing metals:")
    for entry in missing_metals:
        print(', '.join([f"{k}: {v}" for k, v in entry.items()]))
    missing_metals_df = pd.DataFrame(missing_metals)
    missing_metals_df.to_csv('missing_metals.csv', index=False)
