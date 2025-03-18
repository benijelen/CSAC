import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Prompt user for the database file path
database_path = input("Please enter the path to your DATABASE.csv file: ")

# Load CSV file
data = pd.read_csv(database_path)

# Remove specific rows from the database
# rows_to_remove = [188, 189, 190, 191, 231, 237, 241, 247, 368]
# data = data.drop(rows_to_remove)

# Extract necessary columns
average_hydropathy = data['Average Hydropathy']
oxygen_tol = data['OxygenTOL']

# Create a scatter plot with jitter
plt.figure(figsize=(10, 6))
sns.stripplot(x=oxygen_tol, y=average_hydropathy, hue=oxygen_tol, palette='viridis', 
             jitter=0.1, size=8, alpha=0.6, legend=False)

# Adjust x-axis to reduce whitespace
plt.margins(x=0.1)

# Set plot labels and title
plt.xlabel('Oxygen Tolerance', fontsize=16)
plt.ylabel('Average Hydropathy (Kyte-Doolittle)', fontsize=16)
plt.title('Scatter Plot of Average Hydropathy by Oxygen Tolerance', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# Adjust layout to prevent legend cutoff
plt.tight_layout()

# Show plot
plt.show()

# Perform pairwise t-tests for each pair of Oxygen Tolerance groups
groups = data.groupby('OxygenTOL')
oxygen_tol_values = data['OxygenTOL'].unique()

for i, oxygen_tol_1 in enumerate(oxygen_tol_values):
    for oxygen_tol_2 in oxygen_tol_values[i+1:]:
        group1 = groups.get_group(oxygen_tol_1)['Average Hydropathy']
        group2 = groups.get_group(oxygen_tol_2)['Average Hydropathy']
        t_statistic, p_value = stats.ttest_ind(group1, group2)
        print(f"{oxygen_tol_1} vs {oxygen_tol_2} - T-Statistic: {t_statistic}, P-Value: {p_value}")

# Print all rows with Average Hydropathy less than -2.0
low_hydropathy_rows = data[data['Average Hydropathy'] < -2.0]
print("Rows with Average Hydropathy less than -2.0:")
print(low_hydropathy_rows)

# Create a boxplot of Average Hydropathy by Oxygen Tolerance
plt.figure(figsize=(10, 6))
sns.boxplot(x=oxygen_tol, y=average_hydropathy, palette='viridis')

# Set plot labels and title
plt.xlabel('Oxygen Tolerance')
plt.ylabel('Average Hydropathy')
plt.title('Box Plot of Average Hydropathy by Oxygen Tolerance')

# Show plot
plt.show()

# Print PDB IDs for aerobic proteins with high hydropathy
high_hydropathy_aerobic = data[(data['OxygenTOL'] == 'aerobic') & (data['Average Hydropathy'] > 2.1)]
if not high_hydropathy_aerobic.empty:
    print("\nAerobic PDB IDs with Average Hydropathy > 2.1:")
    print(high_hydropathy_aerobic['PDB ID'].tolist())
else:
    print("\nNo aerobic PDB IDs found with Average Hydropathy > 2.1")
