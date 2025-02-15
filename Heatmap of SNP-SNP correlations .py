#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np

# Example SNP dataset (Replace this with actual data)
df_snp = pd.DataFrame({
    'SNP1': ['A/A', 'A/G', 'G/G', 'A/G'],
    'SNP2': ['C/C', 'C/T', 'T/T', 'C/T'],
    'SNP3': ['G/G', 'A/G', 'A/A', 'A/G']
})

# Function to encode SNPs numerically
def encode_snp(genotype):
    mapping = {'A/A': 0, 'A/G': 1, 'G/G': 2, 
               'C/C': 0, 'C/T': 1, 'T/T': 2}
    return mapping.get(genotype, np.nan)

# Apply encoding
snp_encoded = df_snp.applymap(encode_snp)

print(snp_encoded.head())  # Check encoded SNP data


# In[4]:


# Compute correlation
snp_corr = snp_encoded.corr()

print(snp_corr)  # Check correlation values


# In[5]:


import seaborn as sns
import matplotlib.pyplot as plt

# Plot heatmap
plt.figure(figsize=(6, 5))
sns.heatmap(snp_corr, annot=True, cmap="coolwarm", fmt=".2f", linewidths=0.5)
plt.title("SNP Correlation Heatmap")
plt.show()


# In[ ]:




