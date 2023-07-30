#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
get_ipython().system('{sys.executable} -m pip -q install palantir fa2')


# In[2]:


import warnings
warnings.filterwarnings("ignore")
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import palantir
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
## fix palantir breaking down some plots
import seaborn
seaborn.reset_orig()
get_ipython().run_line_magic('matplotlib', 'inline')

sc.set_figure_params()
scf.set_figure_pubready()


# In[3]:


mg = sc.read_h5ad(r"/mnt/g/Levi_Jon/combined.h5ad")


# In[26]:


mgtotal = sc.read_h5ad(r"/mnt/g/Levi_Jon/combined.h5ad")


# In[4]:


mg.obs


# In[5]:


sc.pp.pca(mg)
pca_projections = pd.DataFrame(mg.obsm["X_pca"],index=mg.obs_names)


# In[6]:


dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# In[7]:


# generate neighbor draph in multiscale diffusion space
mg.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg,n_neighbors=30,use_rep="X_palantir")


# In[8]:


# draw ForceAtlas2 embedding using 3 first PCs as initial positions
mg.obsm["X_pca2d"]=mg.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg, layout = 'fa',init_pos='X_pca2d')


# In[9]:


sc.pl.draw_graph(mg,color=['orig.ident'], layout = 'fa') #WE USE THAT


# In[21]:


sc.pl.draw_graph(mg,color=['seurat_clusters'], layout = 'fa') #WE USE THAT


# In[22]:


old_to_new = {
'0':'Transitory',
'1':'Transitory',
'2':'Transitory',
'3':'Transitory',
'4':'Partially muller glia',
'5':'Microglia',
'6':'RGC',
'7':'Cones',
'8':'Mature muller glia',
'9':'AC',
'10':'Early neurons'}
mg.obs['EK_anno'] = (
mg.obs['seurat_clusters']
.map(old_to_new).astype('category')
)


# In[15]:


mg.obs


# In[16]:


import pandas as pd


# In[20]:


mg.obs['seurat_clusters'] = mg.obs['seurat_clusters'].astype(str)


# In[23]:


sc.pl.draw_graph(mg,color=['EK_anno'], layout = 'fa') #WE USE THAT


# In[24]:


def cluster_small_multiples(
    adata, clust_key, size=60, frameon=False, legend_loc=None, **kwargs
):
    tmp = adata.copy()

    for i, clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype("category")
        tmp.uns[clust + "_colors"] = ["#d3d3d3", adata.uns[clust_key + "_colors"][i]]

    sc.pl.draw_graph(
        tmp,
        groups=tmp.obs[clust].cat.categories[1:].values,
        color=adata.obs[clust_key].cat.categories.tolist(),
        size=size,
        frameon=frameon,
        legend_loc=legend_loc,
        **kwargs,
    )


# In[25]:


cluster_small_multiples(mg, clust_key = 'EK_anno', layout = 'fa')


# In[27]:


sc.pp.highly_variable_genes(mgtotal, n_top_genes=1500, flavor='cell_ranger')


# In[28]:


mgtotal.var


# In[30]:


sc.pp.pca(mgtotal)
pca_projections = pd.DataFrame(mgtotal.obsm["X_pca"],index=mgtotal.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)
mgtotal.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mgtotal,n_neighbors=30,use_rep="X_palantir")
mgtotal.obsm["X_pca2d"]=mgtotal.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mgtotal, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(mgtotal,color=['seurat_clusters'], layout = 'fa') #WE USE THAT


# In[31]:


mgtotal.obs['seurat_clusters'] = mgtotal.obs['seurat_clusters'].astype(str)


# In[32]:


old_to_new = {
'0':'Transitory',
'1':'Transitory',
'2':'Transitory',
'3':'Transitory',
'4':'Partially muller glia',
'5':'Microglia',
'6':'RGC',
'7':'Cones',
'8':'Mature muller glia',
'9':'AC',
'10':'Early neurons'}
mgtotal.obs['EK_anno'] = (
mgtotal.obs['seurat_clusters']
.map(old_to_new).astype('category')
)


# In[33]:


sc.pl.draw_graph(mgtotal,color=['seurat_clusters'], layout = 'fa') #WE USE THAT


# In[34]:


sc.pl.draw_graph(mgtotal,color=['EK_anno'], layout = 'fa') #WE USE THAT


# In[ ]:


mgtotal.obsm['X_draw_graph_fa'] = mg.obsm['X_draw_graph_fa']


# In[35]:


cluster_small_multiples(mgtotal, clust_key = 'EK_anno', layout = 'fa')


# In[36]:


import scvelo as scv


# In[38]:


I = scv.read("/mnt/g/Levi_Jon/early/velocyto/early.loom", cache = True)
UI = scv.read("/mnt/g/Levi_Jon/late/Ascl1_Atoh1_late_timepoint/velocyto/Ascl1_Atoh1_late_timepoint.loom", cache = True)


# In[40]:


early = mg[mg.obs['orig.ident'].isin(['early'])]


# In[41]:


late = mg[mg.obs['orig.ident'].isin(['late'])]


# In[60]:


loom = I.concatenate([UI])


# In[42]:


early_velo = scv.utils.merge(early, I)
late_velo = scv.utils.merge(late, UI)


# In[61]:


mg_velo = scv.utils.merge(mg, loom)


# In[43]:


scv.pl.proportions(early_velo)


# In[44]:


scv.pl.proportions(late_velo)


# In[62]:


scv.pl.proportions(mg_velo)


# In[45]:


scv.pp.filter_and_normalize(early_velo, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(early_velo, n_pcs=30, n_neighbors=30)
scv.pp.filter_and_normalize(late_velo, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(late_velo, n_pcs=30, n_neighbors=30)


# In[63]:


scv.pp.filter_and_normalize(mg_velo, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(mg_velo, n_pcs=30, n_neighbors=30)


# In[46]:


scv.tl.velocity(early_velo)
scv.tl.velocity(late_velo)


# In[64]:


scv.tl.velocity(mg_velo)


# In[49]:


scv.tl.velocity_graph(early_velo)
scv.tl.velocity_graph(late_velo)


# In[65]:


scv.tl.velocity_graph(mg_velo)


# In[48]:


sc.pp.neighbors(late_velo)


# In[53]:


sc.set_figure_params(figsize = [8,8])
sc.pl.draw_graph(early_velo,color=['EK_anno'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[66]:


sc.set_figure_params(figsize = [8,8])
sc.pl.draw_graph(mg_velo,color=['EK_anno'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[54]:


sc.set_figure_params(figsize = [8,8])
sc.pl.draw_graph(late_velo,color=['EK_anno'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[55]:


scv.pl.velocity_embedding_grid(early_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[56]:


scv.pl.velocity_embedding_grid(late_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[58]:


scv.pl.velocity_embedding_stream(early_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, density = 2, add_margin=0.2, arrow_size = 2, legend_loc = 'right margin')


# In[59]:


scv.pl.velocity_embedding_stream(late_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, density = 2, add_margin=0.2, arrow_size = 2, legend_loc = 'right margin')


# In[86]:


scv.pl.velocity_embedding_grid(mg_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[77]:


scv.pl.velocity_embedding_stream(mg_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, density = 2, add_margin=0.2, arrow_size = 2, legend_loc = 'right margin')


# In[78]:


scv.tl.velocity_confidence(mg_velo)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(mg_velo, c=keys, cmap='coolwarm', perc=[5, 95], basis='draw_graph_fa' )


# In[79]:


scv.tl.velocity_confidence(early_velo)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(early_velo, c=keys, cmap='coolwarm', perc=[5, 95], basis='draw_graph_fa' )


# In[80]:


scv.tl.velocity_confidence(late_velo)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(late_velo, c=keys, cmap='coolwarm', perc=[5, 95], basis='draw_graph_fa' )


# In[81]:


scv.tl.velocity_pseudotime(mg_velo)
scv.pl.scatter(mg_velo, color='velocity_pseudotime', cmap='gnuplot', basis='draw_graph_fa')


# In[84]:


scv.tl.velocity_pseudotime(early_velo)
scv.pl.scatter(early_velo, color='velocity_pseudotime', cmap='gnuplot', basis='draw_graph_fa')


# In[83]:


sc.pp.neighbors(early_velo)


# In[85]:


scv.tl.velocity_pseudotime(late_velo)
scv.pl.scatter(late_velo, color='velocity_pseudotime', cmap='gnuplot', basis='draw_graph_fa')


# In[87]:


scv.tl.recover_dynamics(early_velo)
scv.tl.recover_dynamics(late_velo)
scv.tl.recover_dynamics(mg_velo)


# In[93]:


scv.tl.recover_dynamics(late_velo)


# In[94]:


scv.tl.velocity(late_velo, mode='dynamical')
scv.tl.velocity_graph(late_velo)


# In[88]:


scv.tl.velocity(early_velo, mode='dynamical')
scv.tl.velocity_graph(early_velo)
scv.tl.velocity(late_velo, mode='dynamical')
scv.tl.velocity_graph(late_velo)
scv.tl.velocity(mg_velo, mode='dynamical')
scv.tl.velocity_graph(mg_velo)


# In[89]:


scv.pl.velocity_embedding_grid(mg_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[90]:


scv.pl.velocity_embedding_grid(early_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[95]:


scv.pl.velocity_embedding_grid(late_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[96]:


scv.tl.latent_time(early_velo)
scv.pl.scatter(early_velo, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[97]:


scv.tl.latent_time(late_velo)
scv.pl.scatter(late_velo, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[98]:


scv.tl.latent_time(mg_velo)
scv.pl.scatter(mg_velo, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[99]:


mg_velo.uns['neighbors']['distances'] = mg_velo.obsp['distances']
mg_velo.uns['neighbors']['connectivities'] = mg_velo.obsp['connectivities']

scv.tl.paga(mg_velo, groups='EK_anno')
df = scv.get_df(mg_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[100]:


early_velo.uns['neighbors']['distances'] = early_velo.obsp['distances']
early_velo.uns['neighbors']['connectivities'] = early_velo.obsp['connectivities']

scv.tl.paga(early_velo, groups='EK_anno')
df = scv.get_df(early_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[101]:


late_velo.uns['neighbors']['distances'] = late_velo.obsp['distances']
late_velo.uns['neighbors']['connectivities'] = late_velo.obsp['connectivities']

scv.tl.paga(late_velo, groups='EK_anno')
df = scv.get_df(late_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[102]:


scv.pl.paga(mg_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[103]:


scv.pl.paga(early_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[104]:


scv.pl.paga(late_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[105]:


sc.tl.leiden(early_velo)


# In[106]:


early_velo.uns['neighbors']['distances'] = early_velo.obsp['distances']
early_velo.uns['neighbors']['connectivities'] = early_velo.obsp['connectivities']

scv.tl.paga(early_velo, groups='leiden')
df = scv.get_df(early_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[107]:


scv.pl.paga(early_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[108]:


import cellrank as cr


# In[116]:


cr.tl.terminal_states(early_velo, cluster_key="leiden", weight_connectivities=0.2)


# In[117]:


cr.pl.terminal_states(early_velo, basis = 'draw_graph_fa')


# In[118]:


cr.tl.initial_states(early_velo, cluster_key="leiden")
cr.pl.initial_states(early_velo, discrete=True, basis = 'draw_graph_fa')


# In[119]:


cr.tl.lineages(early_velo)
cr.pl.lineages(early_velo, same_plot=False, basis = 'draw_graph_fa')


# In[120]:


scv.tl.recover_latent_time(
    early_velo, root_key="initial_states_probs", end_key="terminal_states_probs"
)


# In[121]:


scv.tl.paga(
    early_velo,
    groups="leiden",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime",
)


# In[122]:


cr.pl.cluster_fates(
    early_velo,
    mode="paga_pie",
    cluster_key="leiden",
    basis="draw_graph_fa",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=1,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA",
)


# In[123]:


# compue DPT, starting from CellRank defined root cell
root_idx = np.where(early_velo.obs["initial_states"] == "1")[0][0]
early_velo.uns["iroot"] = root_idx
sc.tl.dpt(early_velo)

scv.pl.scatter(
    early_velo, basis = 'draw_graph_fa',
    color=["leiden", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["leiden", "root cell", "latent time", "dpt pseudotime"],
)


# In[124]:


cr.tl.lineage_drivers(early_velo)


# In[128]:


cr.pl.lineage_drivers(early_velo, lineage="Lineage 0", n_genes=5, basis = 'draw_graph_fa')


# In[129]:


cr.pl.lineage_drivers(early_velo, lineage="Lineage 1", n_genes=5, basis = 'draw_graph_fa')


# In[130]:


cr.pl.lineage_drivers(early_velo, lineage="Lineage 2", n_genes=5, basis = 'draw_graph_fa')


# In[134]:


model = cr.ul.models.GAM(early_velo)
cr.pl.gene_trends(
    early_velo,
    model=model,
    data_key="X",
    genes=["Mef2c"],
    ncols=3,
    time_key="latent_time",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 4),
    n_test_points=200
)


# In[136]:


cr.pl.heatmap(
    early_velo,
    model,
    genes=early_velo.varm['terminal_lineage_drivers']["Lineage_0_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="Lineage 0",
    n_jobs=1,
    backend="loky",
)


# In[141]:


cr.pl.circular_projection(early_velo, keys="EK_anno", legend_loc="right")


# In[ ]:




