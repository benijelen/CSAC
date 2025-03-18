import pandas as pd
from collections import Counter
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

# Prompt user for input file path
input_file = input("Please enter the path to your DATABASE.csv file: ")

# Load data from CSV file
data = pd.read_csv(input_file)

# Extract necessary columns, checking if they exist
required_cols = ['Final Nearby Residues', 'Metal', 'ligand type', 'Chain']
for col in required_cols:
    if col not in data.columns:
        raise ValueError(f"Required column '{col}' not found in CSV file")

# Get metadata columns (all columns before 'Chain')
chain_idx = data.columns.get_loc('Chain')
metadata_cols = list(data.columns[:chain_idx])

# Extract data
lists = ["".join(filter(str.isupper, seq)) for seq in data['Final Nearby Residues']]
metals = data['Metal'].tolist()
ligand_types = data['ligand type'].tolist()
metadata = data[metadata_cols].to_dict('records')

# Create frequency vectors
freq_vectors = [Counter(seq) for seq in lists]

# Convert to DataFrame
df = pd.DataFrame(freq_vectors).fillna(0)

# Standardize the numerical data in the vectors by rows
df = df.div(df.sum(axis=1), axis=0).fillna(0)

# Add metadata and categories to the DataFrame
for col in metadata_cols:
    df[col] = data[col]
df['Metal'] = metals
df['Ligand Type'] = ligand_types

# Save frequency vectors to a new CSV file
output_file = input_file.replace('DATABASE.csv', 'FREQ_VECTOR.csv')
df.to_csv(output_file, index=False)
