import os
import pandas as pd
import requests  # Added requests library import
from Bio.PDB import PDBParser, Selection, NeighborSearch, is_aa  # Added is_aa import
import numpy as np

# Prompt user for input file path
input_file = input("Please enter the path to your CSAC_input.csv file: ")

# Read the existing PDB ID and Metabolism data with header
pdb_list_df = pd.read_csv(input_file)

# Prompt user for metals to search for
metals_input = input("Enter metals to search for (2-letter element codes in capitals, separated by commas. Example: FE, NI, MN): ")
metals = [metal.strip() for metal in metals_input.split(',')]
nife_data = []
# Initialize an empty DataFrame to store the analysis data
# Initialize a new DataFrame to store the NiFe database
additional_columns = ['Chain', 'Metal', 'ligand type', 'ligand radius', 'Nearby Atoms', 'Atom Types', 'Average Electronegativity', 'Electronegativity SD', 'Final Nearby Residues', 'Average Hydropathy', 'Hydropathy SD']
columns = list(pdb_list_df.columns) + additional_columns
dataframe = pd.DataFrame(columns=columns)

# Keep track of failed downloads and missing metals
failed_downloads = []
missing_metals = []

# Function to calculate the center of a cluster
def calculate_center(atoms):
    coords = [atom.get_coord() for atom in atoms]
    center = np.sum(coords, axis=0) / len(coords)
    return center

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

# Function to find nearby residues given a structure, a center point, and a distance threshold
# This version counts residues if any of their atoms are within the threshold
def find_any_nearby_residues(structure, center, threshold):
    atoms = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atoms)
    nearby_atoms = ns.search(center, threshold, 'A')
    residues = set()
    for atom in nearby_atoms:
        residue = atom.get_parent()
        if is_aa(residue):  # Fixed the NameError by using is_aa from Bio.PDB
            residues.add(residue)
    return residues

# Function to find nearby atoms given a structure, a center point, and a distance threshold
def find_nearby_atoms(structure, center, threshold=4):
    atoms = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atoms)
    nearby_atoms = ns.search(center, threshold, 'A')
    return nearby_atoms

# Iterate over each PDB ID
for index, row in pdb_list_df.iterrows():
    pdb_id = row['PDB ID']
    print(pdb_id)
    
    # Get all columns except PDB ID and store as metadata
    metadata = {col: row[col] for col in pdb_list_df.columns if col != 'PDB ID'}
    for key, value in metadata.items():
        print(f"{key}: {value}")
    
    found_metals = set()  # Track which metals were found in this structure
    
    # Check if the PDB file has already been downloaded
    if os.path.exists(f"{pdb_id}.pdb"):
        # Parse the existing structure file
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
    else:
        # Define the URL for the protein structure
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        # Send a GET request to the URL
        response = requests.get(url)
        # Check if the request was successful
        if response.status_code == 200:
            # Write the content of the response to a file
            with open(f"{pdb_id}.pdb", "w") as file:
                file.write(response.text)
            # Parse the downloaded structure
            parser = PDBParser()
            structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
        else:
            print(f"Failed to download the file for PDB ID: {pdb_id}")
            failed_downloads.append(pdb_id)
            continue

    # Iterate over each atom in the structure to find occurrences of metals
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element in metals:
                        found_metals.add(atom.element)
                        print(atom)
                        print(chain)
                        print(residue)
                        cluster_atoms = [atom]
                        center = calculate_center(cluster_atoms)
                        nearby_atoms = find_nearby_atoms(structure, center)
                        # Remove the metal atom itself from the nearby atoms list
                        nearby_atoms = [a for a in nearby_atoms if a != atom]
                        nearby_atoms_identifiers = [f"{a.get_parent().get_resname()}{a.get_serial_number()}" for a in nearby_atoms]
                        nearby_atom_types = [a.element for a in nearby_atoms]
                        # Calculate average electronegativity
                        avg_electronegativity = np.mean([electronegativity[elem] for elem in nearby_atom_types if elem in electronegativity])
                        electronegativity_sd = np.std([electronegativity[elem] for elem in nearby_atom_types if elem in electronegativity])  # Calculate standard deviation of electronegativity
                        ligand_radius = np.max([np.linalg.norm(atom.get_coord() - center) for atom in residue])+0.000001
                        nearby_residues = find_any_nearby_residues(structure, center, 5)
                        nearby_residues_one_letter = [aa3to1[res.resname] for res in nearby_residues if res.resname in aa3to1]
                        nearby_ligres = find_any_nearby_residues(structure, center, ligand_radius+2)
                        nearby_ligres_one_letter = [aa3to1[res.resname] for res in nearby_ligres if res.resname in aa3to1]
                        # Determine the final nearby residues based on ligand radius
                        final_nearby_residues = nearby_residues_one_letter if ligand_radius == 0.000001 else nearby_ligres_one_letter
                        scores = [hydropathy[res] for res in final_nearby_residues]
                        avg_score = np.mean(scores)
                        std_dev = np.std(scores)
                        # Add the data to the DataFrame
                        new_row = {
                            'PDB ID': pdb_id,
                            'Chain': chain.id,
                            'Metal': atom.element,
                            'ligand type': residue.resname,
                            'ligand radius': ligand_radius,
                            'Nearby Atoms': ', '.join(nearby_atoms_identifiers),
                            'Atom Types': ', '.join(nearby_atom_types),
                            'Average Electronegativity': avg_electronegativity,
                            'Electronegativity SD': electronegativity_sd,  # Added column for Electronegativity SD
                            'Final Nearby Residues': ', '.join(final_nearby_residues),
                            'Average Hydropathy': avg_score,
                            'Hydropathy SD': std_dev,
                            'Hydropathy Scores': scores  # Added column for individual hydropathy scores
                        }
                        # Add metadata from input CSV
                        new_row.update(metadata)
                        dataframe = pd.concat([dataframe, pd.DataFrame([new_row])], ignore_index=True)
    
    # Check which metals from the search list were not found
    missing = set(metals) - found_metals
    if missing:
        for metal in missing:
            missing_entry = {'PDB ID': pdb_id, 'Missing Metal': metal}
            missing_entry.update(metadata)  # Add metadata from input CSV
            missing_metals.append(missing_entry)

    # Save the dataframe to a CSV file
    input_filename = os.path.splitext(os.path.basename(input_file))[0]  # Get base name of input file without extension
    output_filename = input_filename + '_' + '_'.join(metals) + '_DATABASE.csv'
    dataframe.to_csv(output_filename, index=False)

# Print and save failed downloads
if failed_downloads:
    print("\nFailed downloads:")
    for pdb_id in failed_downloads:
        print(pdb_id)
    
    # Save failed downloads to a file
    with open('failed_downloads.txt', 'w') as f:
        f.write('\n'.join(failed_downloads))

# Print structures missing metals
if missing_metals:
    print("\nStructures missing metals:")
    for entry in missing_metals:
        print(', '.join([f"{k}: {v}" for k, v in entry.items()]))
    
    # Save missing metals to a file
    missing_metals_df = pd.DataFrame(missing_metals)
    missing_metals_df.to_csv('missing_metals.csv', index=False)
