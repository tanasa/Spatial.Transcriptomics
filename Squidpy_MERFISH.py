#!/usr/bin/env python
# coding: utf-8

# In[1]:


# https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_merfish.html
# https://www.science.org/doi/10.1126/science.aau5324


# In[2]:


import squidpy as sq
import matplotlib.pyplot as plt
import scanpy as sc


# In[3]:


print("MERFISH analysis")
print("The dataset consists of consecutive slices from the mouse hypothalamic preoptic region.")


# In[4]:


sc.logging.print_header()
print(f"squidpy=={sq.__version__}")


# In[ ]:


# load the pre-processed dataset
adata = sq.datasets.merfish()
adata


# In[5]:


print("Visualize the 3D of the slice.")


# In[6]:


sc.pl.embedding(adata, basis="spatial3d", projection="3d", color="Cell_class")


# In[7]:


print('Visualize a single slide: the slide identifier is stored in adata.obs["Bregma"]')


# In[8]:


sq.pl.spatial_scatter(
    adata[adata.obs.Bregma == -9], shape=None, color="Cell_class", size=1
)


# In[9]:


print('Neighborhood enrichment analysis in 3D')


# In[10]:


print("It is important to consider whether the analysis should be performed on the 3D spatial coordinates or the 2D coordinates for a single slice.")


# In[11]:


# To compute the neighbor graph on the 3D coordinate space, we need to specify spatial_key = "spatial3d". 


# In[12]:


sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial3d")
sq.gr.nhood_enrichment(adata, cluster_key="Cell_class")

plt.figure(figsize=(4, 4))
sq.pl.nhood_enrichment(
    adata, 
    cluster_key="Cell_class", 
    method="single", 
    cmap="inferno", 
    vmin=-50, vmax=100, size=2
)


# In[13]:


# To visualize some of the co-enriched clusters with scanpy.pl.embedding(). We will set na_colors=(1,1,1,0)


# In[14]:


sc.pl.embedding(
    adata,
    basis="spatial3d",
    groups=["OD Mature 1", "OD Mature 2", "OD Mature 4"],
    na_color=(1, 1, 1, 0),
    projection="3d",
    color="Cell_class",
)


# In[15]:


print("To perform differential expression testing with scanpy.tl.rank_genes_groups() and visualize the results")


# In[16]:


sc.tl.rank_genes_groups(adata, groupby="Cell_class")
sc.pl.rank_genes_groups(adata, groupby="Cell_class")


# In[17]:


print("Their expression in 3D")

sc.pl.embedding(adata, basis="spatial3d", projection="3d", color=["Gad1", "Mlc1"])


# In[18]:


print("If the same analysis should be performed on a single slice, then it is advisable to copy the sample of interest in a new anndata. AnnData and use it as a standard 2D spatial data object.")


# In[19]:


adata_slice = adata[adata.obs.Bregma == -9].copy()
sq.gr.spatial_neighbors(adata_slice, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="Cell_class")
sq.pl.spatial_scatter(
    adata_slice,
    color="Cell_class",
    shape=None,
    groups=[
        "Ependymal",
        "Pericytes",
        "Endothelial 2",
    ],
    size=10,
)


# In[20]:


print("Spatially variable genes with spatial autocorrelation statistics")


# In[ ]:





# In[21]:


print(''' 
Two spatial autocorrelation statistics: Moran’s I and Geary’s C. 
They provide a score on the degree of spatial variability of gene expression. 
The statistic as well as the p-value are computed for each gene, and FDR correction is performed. 
''')


# In[22]:


sq.gr.spatial_autocorr(adata_slice, mode="moran")
adata_slice.uns["moranI"].head()
sq.pl.spatial_scatter(
    adata_slice, 
    shape=None, 
    color=["Cd24a", "Necab1", "Mlc1"], 
    size=3
)


# In[ ]:




