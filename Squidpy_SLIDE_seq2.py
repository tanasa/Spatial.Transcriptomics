#!/usr/bin/env python
# coding: utf-8

# In[1]:


# https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_slideseqv2.html


# In[2]:


# Robert R Stickels, Evan Murray, Pawan Kumar, Jilong Li, Jamie L Marshall, Daniela J Di Bella, Paola Arlotta, Evan Z Macosko, and Fei Chen. 
# Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2. Nat. Biotechnol., 2020. doi:10.1038/s41587-020-0739-1.


# In[21]:


import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq

print(f"squidpy=={sq.__version__}")


# In[4]:


# load the pre-processed dataset
adata = sq.datasets.slideseqv2()
adata


# In[5]:


# Let’s visualize cluster annotation in spatial context.


# In[6]:


sq.pl.spatial_scatter(adata, color="cluster", size=1, shape=None)


# In[ ]:





# In[7]:


print('''
NEIGHBOURHOOD ENRICHMENT ANALYSIS :

We can investigate spatial organization of clusters in a quantitative way, by computing a NEIGHBOURHOOD ENRICHMENT SCORE. 

It’s an enrichment score on spatial proximity of clusters: if spots belonging to two different clusters are often close to each other, 
then they will have a high score and can be defined as being ENRICHED. 

On the other hand, if they are far apart, the score will be low and they can be defined as DEPLETED. 

This score is based on a permutation-based test, and you can set the number of permutations with the n_perms argument (default is 1000).
''')


# In[8]:


sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(
    adata, 
    cluster_key="cluster", 
    method="single", 
    cmap="inferno", 
    vmin=-50, vmax=100
)


# In[9]:


# interestingly, there seems to be an enrichment between the Endothelial_Tip, the Ependymal cells. 
# another putative enrichment is between the Oligodendrocytes and Polydendrocytes cells.


# In[10]:


sq.pl.spatial_scatter(
    adata,
    shape=None,
    color="cluster",
    groups=["Endothelial_Tip", "Ependymal", "Oligodendrocytes", "Polydendrocytes"],
    size=3,
)


# In[ ]:





# In[11]:


print("Ripley’s statistics :")


# In[12]:


print('''
Ripley’s statistics allow analyst to evaluate whether a discrete annotation (e.g. cell-type) appears to be clustered, dispersed or randomly distributed 
on the area of interest. 

In Squidpy, there are three closely related Ripley’s statistics, that can be easily computed with squidpy.gr.ripley(). 
Ripley’s L statistics is a variance-stabilized version of the Ripley’s K statistics. 
''')


# In[13]:


mode = "L"
sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
sq.pl.ripley(adata, cluster_key="cluster", mode=mode)


# In[ ]:


mode = "K"
sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
sq.pl.ripley(adata, cluster_key="cluster", mode=mode)


# In[14]:


# some cell-types have a more clustered pattern, like Astrocytes and CA11_CA2_CA3_Subiculum cells, 
# whereas other have a more dispersed pattern, like Mural cells.


# In[15]:


sq.pl.spatial_scatter(
    adata,
    color="cluster",
    groups=["CA1_CA2_CA3_Subiculum", "Astrocytes"],
    size=3,
    shape=None,
)


# In[16]:


sq.pl.spatial_scatter(
    adata,
    color="cluster",
    groups=["Mural", "Astrocytes"],
    size=3,
    shape=None,
)


# In[17]:


print("Ligand-receptor interaction analysis :")


# In[18]:


# res = sq.gr.ligrec(
#    adata,
#    cluster_key="cluster",
#    interactions_params={"resources": "CellPhoneDB"},
#    n_perms=1000,
#    threshold=0.1,
#    copy=True
# )


# In[19]:


print("Available clusters:", adata.obs["cluster"].unique())
print(adata.var.index[:10])  # Check the first 10 gene names


# In[22]:


# Step 1: Load your ligand-receptor interaction file
file_LR = "interaction_input_cellphonedb_version_github.filtered2"  
interactions = pd.read_csv(file_LR, sep="\t")
print(interactions.head(2))


# In[23]:


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


# In[26]:


sq.gr.ligrec(
    adata,
    n_perms = 1000,                             # Number of permutations
    threshold = 0.001,                          # Minimum expression threshold
    use_raw = False,                           # Use normalized data
    interactions = interactions,               # Properly formatted DataFrame
    cluster_key = "cluster"                    # Cluster key in AnnData
)
#   clusters=["Polydendrocytes", "Oligodendrocytes"],


# In[27]:


# adata.uns['ligrec']

# Check if ligand-receptor data is stored in `adata.uns`

# if "ligrec_means" in adata.uns and "ligrec_pvalues" in adata.uns:
#    print("Ligand-receptor interaction data is available.")
#    print("Interaction means matrix:")
#    print(adata.uns["ligrec_means"].head())

#    print("\nInteraction p-values matrix:")
#    print(adata.uns["ligrec_pvalues"].head())
# else:
#    print("Ligand-receptor interaction data is missing. Ensure the database was uploaded or interactions were computed.")


# In[28]:


# Visualize ligand-receptor heatmap
sq.pl.ligrec(
    adata,
    key="ligrec",  # Key where results are stored
    cluster_key="cluster",  # Cluster key in AnnData
    top_n_ligands=5,  # Number of top ligands to show
    top_n_receivers=5,  # Number of top receptors to show
    save="Squidpy_SLIDEseq2_plot.png"  # Saves the plot
)

#   clusters=["Polydendrocytes", "Oligodendrocytes"]


# In[ ]:





# In[29]:


print("Spatially variable genes with spatial autocorrelation statistics")


# In[30]:


print('''
We can investigate spatial variability of gene expression. squidpy.gr.spatial_autocorr() conveniently wraps two spatial autocorrelation statistics: 
Moran’s I and Geary’s C*. They provide a score on the degree of spatial variability of gene expression.
''')


# In[31]:


sq.gr.spatial_autocorr(adata, mode="moran")
adata.uns["moranI"].head(10)


# In[32]:


sq.pl.spatial_scatter(
    adata,
    shape=None,
    color=["Ttr", "Plp1", "Mbp", "Hpca", "Enpp2", "Pcp4", "Sst", "Ptgds", "Nrgn"],
    size=0.1,
)


# In[33]:


sq.gr.spatial_autocorr(adata, mode="geary")
adata.uns["gearyC"].head(20)


# In[ ]:




