import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Prompt user for the database file path
database_path = input("Please enter the path to your DATABASE.csv file: ")

# Read the CSV file
data = pd.read_csv(database_path)

# Count the occurrences of each ligand type
ligand_counts = data['ligand type'].value_counts()

# Calculate percentages
ligand_percentages = ligand_counts / ligand_counts.sum() * 100

# Create a figure and axis
plt.figure(figsize=(12, 6))

# Create a bar plot
sns.barplot(x=ligand_percentages.index, y=ligand_percentages.values)

# Customize the plot
plt.title('Distribution of Ligand Types', fontsize=24)
plt.xlabel('Ligand Type', fontsize=16)
plt.ylabel('Percentage', fontsize=16)
plt.xticks(rotation=45, ha='right', fontsize=16)
plt.yticks(fontsize=18)
plt.tight_layout()

# Enhance the style
sns.set_style("whitegrid")
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Add percentage labels on top of each bar
for i, v in enumerate(ligand_percentages.values):
    plt.text(i, v, f'{v:.1f}%', ha='center', va='bottom', fontsize=16)

# Save the figure as SVG
plt.savefig('ligand_type_distribution.svg', format='svg', bbox_inches='tight')

# Show the plot
plt.show()

