import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy
import igraph
import scvelo as scv
import loompy as lmp
import anndata
import os
#load loom
loom_data = scv.read('/home/shixi7/zhaochuang/project/pymt/PyMT/download_data/loom/PyMT.loom', cache=False)

#rename the barcodes
loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('PyMT:','').replace('x', ''))

#read the metadata
sample_obs = pd.read_csv('./cellID_obs.csv')
cell_umap= pd.read_csv('./cell_embeddings.csv', header=0, names=["Cell ID", "UMAP_1", "UMAP_2"])
cell_clusters = pd.read_csv('./cell_clusters.csv', header=0, names=["Cell ID", "cluster"])
cell_celltype = pd.read_csv('./cell_celltype.csv', header=0, names=["Cell ID", "celltype"])

#filter cell
sample_one = loom_data[np.isin(loom_data.obs.index, sample_obs)]

#create the object
sample_one_index = pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {'CellID':'Cell ID'})
umap_ordered = sample_one_index.merge(cell_umap, on = 'Cell ID')
celltype_ordered = sample_one_index.merge(cell_celltype, on = "Cell ID")
celltype_ordered.head()
clusters_ordered = sample_one_index.merge(cell_clusters, on = "Cell ID")
clusters_ordered.head()
umap_ordered = umap_ordered.iloc[:,1:]
clusters_ordered = clusters_ordered.iloc[:,1:]
celltype_ordered = celltype_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
sample_one.uns['clusters'] = clusters_ordered.values
adata = sample_one
adata.obs['celltype'] = np.ravel(celltype_ordered.values)

#unique Ensembl IDs
adata.var_names_make_unique()
adata.write('Allcelltype_dynamicModel.h5ad', compression = 'gzip')

#filter and normalize
scv.pp.filter_and_normalize(adata,min_shared_counts=30, n_top_genes=2000)

#remove duplicate
scv.pp.remove_duplicate_cells(adata)#去除重名的细胞

#reduction and cluster
scv.pp.pca(adata)
scanpy.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)

#plot
ident_colours = ["#9E0142","#D73027","#F46D43","#5E4FA2","#B3DE69","#41AE76","#08519C"]
scv.pl.velocity_embedding_stream(adata, basis='umap',color = "celltype", palette = ident_colours, size = 20,alpha =1,fontsize=0.5,save='embedding_grid1.svg')

scv.tl.recover_dynamics(adata)
scv.tl.recover_latent_time(adata)
scanpy.pl.violin(adata, keys='latent_time',groupby="celltype",save='scVelo-violin-latent_time1.png')

# Generate plot with UMAP and latent time
scv.pl.velocity_embedding_stream(adata,basis="umap",color="latent_time",title='Myeloid',fontsize=20,legend_fontsize=20,min_mass=2,color_map="plasma",save='scVelo-umap-latent_time1.png')

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5,save='scVelo-uamp-paga1.png')



# Estimate RNA velocity and latent time
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.tl.recover_latent_time(adata)

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5,save='scVelo-uamp-paga2.svg')

#plot
ident_colours = ["#9E0142","#D73027","#F46D43","#5E4FA2","#B3DE69","#41AE76","#08519C"]
scv.pl.velocity_embedding_stream(adata, basis='umap',color = "celltype", palette = ident_colours, size = 20,alpha =1,save='embedding_grid2.svg')

scanpy.pl.violin(adata, keys='latent_time',groupby="celltype",save='scVelo-violin-latent_time2.png')

# Generate plot with UMAP and latent time
scv.pl.velocity_embedding_stream(adata,basis="umap",color="latent_time",title='Myeloid',fontsize=20,legend_fontsize=20,min_mass=2,color_map="plasma",save='scVelo-umap-latent_time2.png')

