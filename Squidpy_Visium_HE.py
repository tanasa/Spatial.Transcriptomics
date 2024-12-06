#!/usr/bin/env python
# coding: utf-8

# In[1]:


# https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html


# In[2]:


import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq


# In[3]:


sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

# load the pre-processed dataset
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()


# In[ ]:





# In[4]:


# To visualize cluster annotation in spatial context with squidpy.pl.spatial_scatter().

sq.pl.spatial_scatter(adata, color="cluster")


# In[5]:


print('''

Image features : 

Visium datasets contain high-resolution images of the tissue that was used for the gene extraction. 
squidpy.im.calculate_image_features() computes image features for each Visium spot and creates an 
obs x features matrix in adata that can then be analyzed together with the obs x gene gene expression matrix.
''')

# Squidpy contains several feature extractors and a flexible pipeline of calculating features of different scales and size


# In[6]:


# calculate features for different scales (higher value means more context)
for scale in [1.0, 2.0]:
    feature_name = f"features_summary_scale{scale}"
    sq.im.calculate_image_features(
        adata,
        img.compute(),
        features="summary",
        key_added=feature_name,
        n_jobs=4,
        scale=scale,
    )


# combine features in one dataframe
adata.obsm["features"] = pd.concat(
    [adata.obsm[f] for f in adata.obsm.keys() if "features_summary" in f],
    axis="columns",
)
# make sure that we have no duplicated feature names in the combined table
adata.obsm["features"].columns = ad.utils.make_index_unique(
    adata.obsm["features"].columns
)


# In[7]:


# helper function returning a clustering

def cluster_features(features: pd.DataFrame, like=None) -> pd.Series:
    """
    Calculate leiden clustering of features.

    Specify filter of features using `like`.
    """
    # filter features
    if like is not None:
        features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    return adata.obs["leiden"]


# calculate feature clusters
adata.obs["features_cluster"] = cluster_features(adata.obsm["features"], like="summary")

# compare feature and gene clusters
sq.pl.spatial_scatter(adata, color=["features_cluster", "cluster"])


# In[ ]:





# In[8]:


print("Spatial statistics and graph analysis")


# In[9]:


print('''

NEIGHBOURHOOD ENRICHMENT :

Computing a neighborhood enrichment can help us identify spots clusters that share a common neighborhood structure across the tissue. 
We can compute such score with the following function: squidpy.gr.nhood_enrichment().

It’s an enrichment score on spatial proximity of clusters: if spots belonging to two different clusters are often close to each other, 
then they will have a high score and can be defined as being ENRICHED. 

On the other hand, if they are far apart, and therefore are seldom a neighborhood, the score will be low and they can be defined as DEPLETED. 

This score is based on a permutation-based test, and you can set the number of permutations with the n_perms argument (default is 1000).
''')


# In[10]:


sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")


# In[11]:


# we find high neighborhood enrichment the Hippocampus region: Pyramidal_layer_dentate_gyrus and Pyramidal_layer clusters seems to be often neighbors 
# with the larger Hippocampus cluster.


# In[12]:


sq.gr.co_occurrence(adata, cluster_key="cluster")
sq.pl.co_occurrence(
    adata,
    cluster_key="cluster",
    clusters="Hippocampus",
    figsize=(8, 4),
)


# In[ ]:





# In[13]:


print("Ligand-receptor interaction analysis :")


# In[14]:


print("Available clusters:", adata.obs["cluster"].unique())
print(adata.var.index[:3])  # Check the first 3 gene names

# Step 1: Load your ligand-receptor interaction file
file_LR = "interaction_input_cellphonedb_version_github.filtered2"  
interactions = pd.read_csv(file_LR, sep="\t")
print(interactions.head(2))

# Step 2: Rename columns to 'ligand' and 'receptor' for Squidpy compatibility

# interactions = interactions.rename(columns={
#    'protein_name_a': 'ligand',
#    'protein_name_b': 'receptor'
# })

interactions.rename(columns = {
               "protein_name_a": "source", 
               "protein_name_b": "target"}, inplace=True)

# Step 3: Inspect the DataFrame to ensure it looks correct
print(interactions.head(2))

if "source" not in interactions.columns or "target" not in interactions.columns:
    raise KeyError("The interactions DataFrame must include 'source' and 'target' columns.")


# In[15]:


sq.gr.ligrec(
    adata,
    n_perms = 100,                             # Number of permutations
    threshold = 0.01,                          # Minimum expression threshold
    use_raw = False,                           # Use normalized data
    interactions = interactions,               # Properly formatted DataFrame
    cluster_key = "cluster"    # Cluster key in AnnData
)


# In[16]:


# Visualize ligand-receptor heatmap
sq.pl.ligrec(
    adata,
    key="ligrec",  # Key where results are stored
    cluster_key="cluster",  # Cluster key in AnnData
    top_n_ligands=5,  # Number of top ligands to show
    top_n_receivers=5,  # Number of top receptors to show
    save="Squidpy_Visium+HE_plot.png"  # Saves the plot
)


# In[17]:


print("Spatially variable genes with Moran's I")
# finding genes that show spatial patterns


# In[18]:


genes = adata[:, adata.var.highly_variable].var_names.values[:1000]
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=genes,
    n_perms=100,
    n_jobs=1,
)


# In[19]:


adata.uns["moranI"].head(10)


# In[20]:


sq.pl.spatial_scatter(
    adata,
    shape=None,
    color=["Olfm1", "Plp1", "Itpka", "Snap25", "Nnat", "Ppp3ca", "Chn1", "Mal", "Tmsb4x", "Cldn11"],
    size=0.1,
)


# In[ ]:





# In[21]:


print("Spatially variable genes with Geary's C")
# finding genes that show spatial patterns


# In[23]:


genes = adata[:, adata.var.highly_variable].var_names.values[:1000]
sq.gr.spatial_autocorr(
    adata,
    mode="geary",
    genes=genes,
    n_perms=100,
    n_jobs=1,
)


# In[24]:


adata.uns["gearyC"].head(10)


# In[ ]:





# In[38]:


print('''

Adding Ripley’s statistics : 

It allows us to evaluate whether a discrete annotation (e.g. cell-type) appears to be clustered, dispersed or randomly distributed 
on the area of interest. 

In Squidpy, there are three closely related Ripley’s statistics, that can be easily computed with squidpy.gr.ripley(). 
Ripley’s L statistics is a variance-stabilized version of the Ripley’s K statistics. 
''')


# In[39]:


# mode = "L"
# sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
# sq.pl.ripley(adata, cluster_key="cluster", mode=mode)


# In[40]:


mode = "F"
sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
sq.pl.ripley(adata, cluster_key="cluster", mode=mode)


# In[41]:


# mode = "G"
# sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
# sq.pl.ripley(adata, cluster_key="cluster", mode=mode)


# In[42]:


import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(15, 7))
mode = "F"

# Assuming 'sq.gr.ripley' and 'sq.pl.ripley' are valid methods
sq.gr.ripley(adata, cluster_key="leiden", mode=mode)
sq.pl.ripley(adata, cluster_key="leiden", mode=mode, ax=ax[0])

sq.pl.spatial_scatter(
    adata,
    color="leiden",
    groups=["0", "1", "3"],
    shape=None,
    size=2,
    ax=ax[1],
)


# In[ ]:




