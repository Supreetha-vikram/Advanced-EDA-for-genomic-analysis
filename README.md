Advanced EDA for Genomic Data Analysis: Identifying Genetic Variations Through Visualization

📌 Project Overview
This project explores genomic data analysis using Exploratory Data Analysis (EDA) and Principal Component Analysis (PCA) to identify genetic variations in SNP (Single Nucleotide Polymorphism) datasets. Clustering methods such as K-Means and DBSCAN are applied to find patterns in genomic data.

📂 Dataset
We used a genomic dataset (PED format) containing SNP data. The dataset includes:
Genotypes (A/G, C/T, etc.) for multiple SNPs
Metadata (Phenotypes, IDs, etc.)

🔹 Data Preprocessing
Loaded PED file and extracted SNP genotypes.
Converted nucleotide bases to numerical format.
Handled missing values.
Scaled SNP data for PCA & clustering.

📊 Exploratory Data Analysis (EDA)
EDA helps understand genetic distributions and variation. We performed:

✅ 1️⃣ Allele Frequency Analysis
Bar plots to visualize allele distributions.
Identified common vs. rare SNPs.

✅ 2️⃣ Correlation Analysis
Heatmap of SNP-SNP correlations to find related genetic markers.

✅ 3️⃣ PCA for Genetic Variation
Applied Principal Component Analysis (PCA) to reduce SNP dimensions.
Visualized PCA clusters (2D & 3D plots).
Checked variance explained by top PCs.

🖥️ PCA Visualization:
fig = px.scatter_3d(df_pca_3d, x="PC1", y="PC2", z="PC3", color="Phenotype")
fig.show()

📌 Results & Insights

PCA revealed distinct genetic variation patterns.

SNP correlation heatmaps showed co-inherited genetic markers.

🔧 Installation & Usage

📌 1️⃣ Install Dependencies
pip install pandas numpy scikit-learn seaborn plotly matplotlib

📌 2️⃣ Run the Notebook
jupyter notebook

Open the notebook and follow the steps.


