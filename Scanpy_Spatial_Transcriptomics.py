#!/usr/bin/env python
# coding: utf-8

# In[1]:


# https://scanpy.readthedocs.io/en/stable/tutorials/spatial/basic-analysis.html
# https://squidpy.readthedocs.io/en/stable/
# used data : https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Human_Lymph_Node


# In[2]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
np.random.seed(42)


# In[3]:


sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# In[ ]:





# In[4]:


print("Reading the data:")


# In[5]:


adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


# In[6]:


str(adata)


# In[ ]:





# In[7]:


print('QC and preprocessing:')


# In[8]:


fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)


# In[9]:


sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)


# In[ ]:





# In[10]:


print('Normalization :')


# In[11]:


sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)


# In[ ]:





# In[12]:


print("Manifold embedding and clustering based on transcriptional similarity :")


# In[13]:


sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.tl.leiden(
    adata,
    key_added="clusters",
    directed=False,
    n_iterations=2
)


# In[14]:


# Visualize the UMAP plot
# sc.pl.umap(adata)


# In[15]:


plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)


# In[16]:


print("Visualization in spatial coordinates :")


# In[17]:


# using sc.pl.spatial :
# img_key: key where the img is stored in the adata.uns element
# crop_coord: coordinates to use for cropping (left, right, top, bottom)
# alpha_img: alpha value for the transcparency of the image
# bw: flag to convert the image into gray scale
# size parameter : it becomes a scaling factor for the spot sizes

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])


# In[18]:


plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)


# In[19]:


# By changing the alpha values of the spots, we can visualize better the underlying tissue morphology from the H&E image.

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,
    size=1.3,
)


# In[20]:


print("Cluster marker genes :")

# Let us further inspect cluster 5, which occurs in small groups of spots across the image.
# Compute marker genes and plot a heatmap with expression levels of its top 10 marker genes across clusters.

plt.rcParams["figure.figsize"] = (4, 4)
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="9", n_genes=10, groupby="clusters")


# In[21]:


# We see that CR2 recapitulates the spatial structure.

sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])
sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)


# In[ ]:





# In[22]:


print("A Merfish Example :")
import os

# Get the current working directory
current_folder = os.getcwd()
print("Current folder path:", current_folder)


# In[23]:


# coordinates = pd.read_excel("./merfish/pnas.1912459116.sd15.xlsx", index_col=0)
# counts = sc.read_csv("./merfish/pnas.1912459116.sd12.csv").transpose()

import pandas as pd

# Read the coordinates file (assuming it's a tab-delimited text file)
coordinates = pd.read_csv("./merfish/pnas.1912459116.sd12.txt", sep='\t', index_col=0)

# Print rows and columns separately
print("Number of rows:", coordinates.shape[0])
print("Number of columns:", coordinates.shape[1])

print(coordinates.head(3))


# In[24]:


counts1 = pd.read_csv("./merfish/pnas.1912459116.sd15.sheet1.txt", sep='\t', index_col=0)
counts1 = counts1.transpose()

# Print rows and columns separately
print("Number of rows:", counts1.shape[0])
print("Number of columns:", counts1.shape[1])
print(counts1.head(3))


# In[25]:


counts2 = pd.read_csv("./merfish/pnas.1912459116.sd15.sheet2.txt", sep='\t', index_col=0)
counts2 = counts2.transpose()

# Print rows and columns separately
print("Number of rows:", counts2.shape[0])
print("Number of columns:", counts2.shape[1])
print(counts2.head(3))


# In[26]:


counts3 = pd.read_csv("./merfish/pnas.1912459116.sd15.sheet3.txt", sep='\t', index_col=0)
counts3 = counts3.transpose()

# Print rows and columns separately
print("Number of rows:", counts3.shape[0])
print("Number of columns:", counts3.shape[1])
print(counts3.head(3))


# In[27]:


# Get the current working directory
current_folder = os.getcwd()
print("Current folder path:", current_folder)


# In[28]:


# Merfish example in : https://scanpy.readthedocs.io/en/stable/tutorials/spatial/basic-analysis.html


# In[ ]:




