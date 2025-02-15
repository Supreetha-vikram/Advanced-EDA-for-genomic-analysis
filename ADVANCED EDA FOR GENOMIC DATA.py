#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# In[2]:


# Define column names (First 6 columns are metadata, rest are genotypic data)
col_names = ["Family_ID", "Individual_ID", "Paternal_ID", "Maternal_ID", "Sex", "Phenotype"]

# Loads only first 6 metadata columns
df_metadata = pd.read_csv("donors.ped", delim_whitespace=True, names=col_names, usecols=range(6))

# Display first few rows
print(df_metadata.head())


# In[3]:


df_full = pd.read_csv("donors.ped", delim_whitespace=True, header=None)

# Assign column names (first 6 are metadata, rest are SNPs)
metadata_cols = ["Family_ID", "Individual_ID", "Paternal_ID", "Maternal_ID", "Sex", "Phenotype"]
snp_cols = [f"SNP_{i}" for i in range(1, len(df_full.columns) - 5)]  # SNPs start from column 7
df_full.columns = metadata_cols + snp_cols

# Display first few rows
print(df_full.head())


# In[4]:


print(df_full.info())  # Check column types
print(df_full.describe())  # Get summary statistics
print(df_full.isnull().sum())  # Check for missing values


# In[8]:


# Select only SNP columns
snp_data = df_full.iloc[:, 6:]  # Excluding metadata columns

# Count unique values (alleles) for each SNP
allele_counts = snp_data.apply(lambda col: col.value_counts())

# Display allele frequencies for first 5 SNPs
print(allele_counts.head())


# In[9]:


# Count missing values in the dataset
missing_values = df_full.isnull().sum()

# Print columns with missing values
print(missing_values[missing_values > 0])


# In[10]:


import matplotlib.pyplot as plt

# Get the first SNP column (example)
snp_column = snp_data.columns[0]

# Count allele occurrences
allele_counts = snp_data[snp_column].value_counts()

# Plot
plt.figure(figsize=(6, 4))
allele_counts.plot(kind='bar', color=['blue', 'orange'])
plt.xlabel("Allele Type")
plt.ylabel("Frequency")
plt.title(f"Allele Frequency Distribution for {snp_column}")
plt.show()


# In[11]:


import seaborn as sns

# Count the number of unique alleles per SNP
num_unique_alleles = snp_data.nunique()

# Plot distribution
plt.figure(figsize=(8, 5))
sns.histplot(num_unique_alleles, bins=20, kde=True, color="purple")
plt.xlabel("Number of Unique Alleles per SNP")
plt.ylabel("Count")
plt.title("Distribution of Unique SNP Variants")
plt.show()


# In[12]:


df_full.to_csv("processed_genetic_data.csv", index=False)
print("Processed genetic data saved as 'processed_genetic_data.csv'.")


# In[13]:


from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder, StandardScaler
import numpy as np

# Select only SNP columns (excluding metadata)
snp_data = df_full.iloc[:, 6:].copy()

# Convert allele pairs into numeric values (A/G → 0, G/G → 1, A/A → 2)
def encode_snp(col):
    unique_vals = col.unique()
    mapping = {val: i for i, val in enumerate(unique_vals)}
    return col.map(mapping)

# Apply encoding to all SNP columns
snp_encoded = snp_data.apply(encode_snp)

# Standardize data before PCA
scaler = StandardScaler()
snp_scaled = scaler.fit_transform(snp_encoded)

print("SNP data successfully encoded and standardized!")


# In[14]:


# Apply PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(snp_scaled)

# Convert to DataFrame
df_pca = pd.DataFrame(pca_result, columns=["PC1", "PC2"])

# Add metadata for visualization
df_pca["Phenotype"] = df_full["Phenotype"]

print(df_pca.head())  # Show first few PCA-transformed rows


# In[15]:


import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(8, 6))
sns.scatterplot(x="PC1", y="PC2", hue=df_pca["Phenotype"], palette="coolwarm", alpha=0.7, data=df_pca)

plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.title("PCA Visualization of SNP Data")
plt.legend(title="Phenotype")
plt.show()


# In[16]:


explained_variance = pca.explained_variance_ratio_
print(f"PC1 explains {explained_variance[0]*100:.2f}% of the variance")
print(f"PC2 explains {explained_variance[1]*100:.2f}% of the variance")


# In[17]:


from sklearn.decomposition import PCA

# Apply PCA with 3 components
pca_3d = PCA(n_components=3)
pca_result_3d = pca_3d.fit_transform(snp_scaled)

# Convert to DataFrame
df_pca_3d = pd.DataFrame(pca_result_3d, columns=["PC1", "PC2", "PC3"])

# Add metadata for visualization
df_pca_3d["Phenotype"] = df_full["Phenotype"]

print(df_pca_3d.head())  # Show first few PCA-transformed rows


# In[18]:


import plotly.express as px

# Create a 3D scatter plot
fig = px.scatter_3d(df_pca_3d, x="PC1", y="PC2", z="PC3", 
                     color=df_pca_3d["Phenotype"].astype(str),  # Convert phenotype to string for color coding
                     title="3D PCA Visualization of SNP Data",
                     labels={"Phenotype": "Phenotype"},
                     opacity=0.7)

fig.show()


# In[19]:


explained_variance = pca_3d.explained_variance_ratio_
print(f"PC1 explains {explained_variance[0]*100:.2f}% of the variance")
print(f"PC2 explains {explained_variance[1]*100:.2f}% of the variance")
print(f"PC3 explains {explained_variance[2]*100:.2f}% of the variance")


# In[20]:


# Standardize the data before PCA
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df_full.select_dtypes(include=[np.number]))

# Apply PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data)

# Convert to DataFrame for visualization
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])

# Scatter plot of PCA results
plt.figure(figsize=(8, 6))
sns.scatterplot(x=pca_df['PC1'], y=pca_df['PC2'])
plt.title("PCA of Genomic Data")
plt.show()


# In[ ]:




