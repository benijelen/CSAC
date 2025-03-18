import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Prompt user for the frequency vector file path
file_path = input("Please enter the path to your FREQ_VECTOR.csv file: ")

# Load the data
data = pd.read_csv(file_path)

# Find index of 'PDB ID' column to separate features and metadata
pdb_id_idx = data.columns.get_loc('PDB ID')

# Separate features (columns before PDB ID) and metadata
features = data.columns[:pdb_id_idx]  # Columns before PDB ID are amino acid frequencies
X = data[features]
y = data.iloc[:, pdb_id_idx:]  # PDB ID and all columns after are metadata

# Perform PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X)
principal_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
final_df = pd.concat([principal_df, y], axis=1)

# Use OxygenTOL for grouping
final_df['Group'] = final_df['OxygenTOL']
# Plot the PCA results
plt.figure(figsize=(10, 6))
sns.scatterplot(x='PC1', y='PC2', hue='Group', data=final_df, 
                palette=['blue', 'orange', 'purple', 'yellow'],
                style='Group', markers=['o', 's', '^', 'D'])
plt.title('PCA of Freq_Vector_FE', fontsize=16)
plt.xlabel('Principal Component 1', fontsize=16)
plt.ylabel('Principal Component 2', fontsize=16)
plt.legend(title='Group', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()

# Perform PCA with all components
pca_all = PCA()
pca_all.fit(X)

# Print the loadings of all principal components
loadings_all = pd.DataFrame(pca_all.components_.T, columns=[f'PC{i+1}' for i in range(len(pca_all.components_))], index=features)
print("Loadings for all Principal Components:")
print(loadings_all)

# Plot a scree plot for all principal components
explained_variance_all = pca_all.explained_variance_ratio_
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(explained_variance_all) + 1), explained_variance_all, marker='o', linestyle='--')
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance Ratio')
plt.xticks(range(1, len(explained_variance_all) + 1))  # Ensure all principal components are shown on the x-axis
plt.show()

# Perform PCA to get PC3 and PC4
pca_4 = PCA(n_components=4)
principal_components_4 = pca_4.fit_transform(X)
principal_df_4 = pd.DataFrame(data=principal_components_4, columns=['PC1', 'PC2', 'PC3', 'PC4'])
final_df_4 = pd.concat([principal_df_4, y], axis=1)

# Combine 'OxygenTOL' and 'Metal' into a single column for better color differentiation
final_df_4['Group'] = final_df_4['OxygenTOL'] + ' - ' + final_df_4['Metal']

# Plot PC3 vs. PC4
plt.figure(figsize=(10, 6))
sns.scatterplot(x='PC3', y='PC4', hue='Group', data=final_df_4, palette=['red', 'blue', 'green', 'yellow'])
plt.title('PCA of Freq_Vector_FE (PC3 vs. PC4)')
plt.xlabel('Principal Component 3')
plt.ylabel('Principal Component 4')
plt.legend(title='Group')
plt.show()

# Scatter plot the "C" values from the freq_vector, split points into columns and colors based on "OxygenTOL"
plt.figure(figsize=(10, 6))
sns.boxplot(x=data['OxygenTOL'], y=data['C'], palette='viridis')
plt.title('Box Plot of "C" Values from Freq_Vector_FE')
plt.xlabel('OxygenTOL')
plt.ylabel('"C" Value')
plt.show()

# Determine whether the "C" values for the two groups are statistically significant
group1 = data[data['OxygenTOL'] == 'anaerobic']['C']
group2 = data[data['OxygenTOL'] == 'aerobic']['C']

t_statistic, p_value = stats.ttest_ind(group1, group2)
print(f"T-Statistic: {t_statistic}, P-Value: {p_value}")

if p_value < 0.05:
    print("The difference between the two groups is statistically significant.")
else:
    print("The difference between the two groups is not statistically significant.")
