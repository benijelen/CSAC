import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Prompt user for the database file path
database_path = input("Please enter the path to your DATABASE.csv file: ")

# Read the CSV file
data = pd.read_csv(database_path)

# Create a figure and axis
plt.figure(figsize=(12, 6))

# Create a violin plot
sns.violinplot(x='OxygenTOL', y='ligand radius', data=data)

# Customize the plot
plt.title('Distribution of Ligand Sizes for Aerobic and Anaerobic Proteins', fontsize=16)
plt.xlabel('Oxygen Tolerance', fontsize=16)
plt.ylabel('Ligand Radius (Angstroms)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

# Enhance the style
sns.set_style("whitegrid")
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Save the figure as SVG
plt.savefig('ligand_size_distribution.svg', format='svg', bbox_inches='tight')

# Show the plot
plt.show()

# Find ligand types with 0 ligand radius
zero_radius_ligands = data[data['ligand radius'] == 0]['ligand type'].unique()

print("Ligand types with 0 ligand radius:")
for ligand in zero_radius_ligands:
    print(f"- {ligand}")

# Get unique ligands and their radii
ligand_radii = data[['ligand type', 'ligand radius']].groupby('ligand type')['ligand radius'].first().reset_index()
ligand_radii = ligand_radii.sort_values('ligand type')

print("\nAll unique ligand types and their radii:")
for _, row in ligand_radii.iterrows():
    print(f"- {row['ligand type']}: {row['ligand radius']:.2f} Å")
