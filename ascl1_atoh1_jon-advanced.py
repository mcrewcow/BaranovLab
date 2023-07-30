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


# In[345]:


mg = sc.read_h5ad(r"/mnt/c/Users/Emil/10X/Levi_Jon/combined.h5ad")


# In[346]:


mgtotal = sc.read_h5ad(r"/mnt/c/Users/Emil/10X/Levi_Jon/combined.h5ad")


# In[347]:


import pandas as pd
mg.obs['seurat_clusters'] = mg.obs['seurat_clusters'].astype(str)
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
'9':'Transitory',
'10':'Transitory'}
mg.obs['EK_anno'] = (
mg.obs['seurat_clusters']
.map(old_to_new).astype('category')
)


# In[348]:


mgtotal.obs['seurat_clusters'] = mgtotal.obs['seurat_clusters'].astype(str)
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
'9':'Transitory',
'10':'Early neurons'}
mgtotal.obs['EK_anno'] = (
mgtotal.obs['seurat_clusters']
.map(old_to_new).astype('category')
)


# In[349]:


mg.obs


# In[350]:


mg = mg[mg.obs['EK_anno'].isin(['Transitory','Partially muller glia','RGC','Cones','Mature muller glia','Early neurons'])]


# In[351]:


mgtotal = mgtotal[mgtotal.obs['EK_anno'].isin(['Transitory','Partially muller glia','RGC','Cones','Mature muller glia','Early neurons'])]


# In[334]:


sc.pp.filter_genes(mg,min_counts=3)


# In[352]:


sc.pp.pca(mg)
pca_projections = pd.DataFrame(mg.obsm["X_pca"],index=mg.obs_names)


# In[353]:


dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# In[354]:


# generate neighbor draph in multiscale diffusion space
mg.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg,n_neighbors=30,use_rep="X_palantir")


# In[355]:


# draw ForceAtlas2 embedding using 3 first PCs as initial positions
mg.obsm["X_pca2d"]=mg.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg, layout = 'fa',init_pos='X_pca2d')


# In[356]:


sc.pl.draw_graph(mg,color=['orig.ident'], layout = 'fa') #WE USE THAT


# In[357]:


sc.pl.draw_graph(mg,color=['seurat_clusters'], layout = 'fa') #WE USE THAT


# In[358]:


sc.pl.draw_graph(mg,color=['EK_anno'], layout = 'fa') #WE USE THAT


# In[359]:


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


# In[360]:


cluster_small_multiples(mg, clust_key = 'EK_anno', layout = 'fa')


# In[678]:


sc.pp.filter_genes(mgtotal,min_counts=3)
sc.pp.highly_variable_genes(mgtotal, n_top_genes=1500, flavor='cell_ranger')


# In[679]:


mgtotal.var


# In[680]:


sc.pp.pca(mgtotal)
pca_projections = pd.DataFrame(mgtotal.obsm["X_pca"],index=mgtotal.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)
mgtotal.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mgtotal,n_neighbors=30,use_rep="X_palantir")
mgtotal.obsm["X_pca2d"]=mgtotal.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mgtotal, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(mgtotal,color=['seurat_clusters'], layout = 'fa') #WE USE THAT


# In[681]:


sc.pl.draw_graph(mgtotal,color=['seurat_clusters'], layout = 'fa') #WE USE THAT


# In[682]:


sc.pl.draw_graph(mgtotal,color=['EK_anno'], layout = 'fa') #WE USE THAT


# In[345]:


sc.pl.draw_graph(mgtotal,color=['EK_anno'], layout = 'fa', save = 'ascl1_total_fa.pdf') #WE USE THAT


# In[613]:


sc.pl.draw_graph(mg,color=['EK_anno'], layout = 'fa', save = 'umap_total.pdf')


# In[614]:


sc.pl.draw_graph(early,color=['EK_anno'], layout = 'fa', save = 'early_umap_total.pdf')


# In[615]:


sc.pl.draw_graph(late,color=['EK_anno'], layout = 'fa', save = 'late_umap_total.pdf')


# In[347]:


mg.obsm['X_draw_graph_fa'] = mgtotal.obsm['X_draw_graph_fa']


# In[348]:


cluster_small_multiples(mgtotal, clust_key = 'EK_anno', layout = 'fa')


# In[362]:


mg.layers['scaled'] = sc.pp.scale(mg, copy=True).X


# In[363]:


sc.pl.matrixplot(mg, ['Trpm1','Nrxn3','Pde6h','Rcvrn','Rlbp1','Id1','Id3','Caln1'], 'EK_anno',
                 colorbar_title='mean z-score',  vmin=-2, vmax=2, cmap='RdBu_r', layer = 'scaled', categories_order = ['RGC','Cones','Mature muller glia','Partially muller glia',
                                                                                                                   'Transitory'], figsize = [4,1], save = 'heatmap_levi1.pdf')


# In[124]:


import scvelo as scv


# In[125]:


I = scv.read("/mnt/c/Users/Emil/10X/Levi_Jon/early.loom", cache = True)
UI = scv.read("/mnt/c/Users/Emil/10X/Levi_Jon/late.loom", cache = True)


# In[126]:


early = mg[mg.obs['orig.ident'].isin(['early'])]


# In[127]:


late = mg[mg.obs['orig.ident'].isin(['late'])]


# In[128]:


early_velo = scv.utils.merge(early, I)


# In[129]:


late_velo = scv.utils.merge(late, UI)


# In[130]:


scv.pl.proportions(early_velo)


# In[131]:


scv.pl.proportions(late_velo)


# In[132]:


scv.pp.neighbors(late_velo)


# In[133]:


scv.pp.neighbors(early_velo)


# In[134]:


scv.pp.filter_and_normalize(late_velo, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(late_velo, n_pcs=30, n_neighbors=30)


# In[ ]:





# In[135]:


scv.pp.filter_and_normalize(early_velo, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(early_velo, n_pcs=30, n_neighbors=30)


# In[136]:


scv.tl.velocity(late_velo)


# In[137]:


scv.tl.velocity(early_velo)


# In[138]:


scv.tl.velocity_graph(early_velo)


# In[139]:


scv.tl.velocity_graph(late_velo)


# In[140]:


sc.set_figure_params(figsize = [8,8])
sc.pl.draw_graph(early_velo,color=['EK_anno'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[141]:


sc.set_figure_params(figsize = [8,8])
sc.pl.draw_graph(late_velo,color=['EK_anno'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[142]:


scv.pl.velocity_embedding_grid(early_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[143]:


scv.pl.velocity_embedding_grid(late_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[157]:


scv.pl.velocity_embedding_stream(early_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, density = 2, add_margin=0.2, arrow_size = 2, legend_loc = 'right margin', linewidth = 2, max_length = 8, save = 'velo_early_new.svg')


# In[158]:


scv.pl.velocity_embedding_stream(late_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, density = 2, add_margin=0.2, arrow_size = 2, legend_loc = 'right margin', linewidth = 2, max_length = 8, save = 'velo_late_new.svg')


# In[146]:


scv.tl.velocity_confidence(early_velo)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(early_velo, c=keys, cmap='coolwarm', perc=[5, 95], basis='draw_graph_fa' )


# In[147]:


scv.tl.velocity_confidence(late_velo)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(late_velo, c=keys, cmap='coolwarm', perc=[5, 95], basis='draw_graph_fa' )


# In[148]:


scv.tl.velocity_pseudotime(early_velo)
scv.pl.scatter(early_velo, color='velocity_pseudotime', cmap='gnuplot', basis='draw_graph_fa')


# In[149]:


scv.tl.velocity_pseudotime(late_velo)
scv.pl.scatter(late_velo, color='velocity_pseudotime', cmap='gnuplot', basis='draw_graph_fa')


# In[150]:


scv.tl.recover_dynamics(early_velo)


# In[151]:


scv.tl.recover_dynamics(late_velo)


# In[152]:


scv.tl.velocity(early_velo, mode='dynamical')
scv.tl.velocity_graph(early_velo)
scv.tl.velocity(late_velo, mode='dynamical')
scv.tl.velocity_graph(late_velo)


# In[153]:


scv.pl.velocity_embedding_grid(early_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[154]:


scv.pl.velocity_embedding_grid(late_velo, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[765]:


scv.tl.latent_time(early_velo)
scv.pl.scatter(early_velo, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[766]:


scv.tl.latent_time(late_velo)
scv.pl.scatter(late_velo, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[767]:


early_velo.uns['neighbors']['distances'] = early_velo.obsp['distances']
early_velo.uns['neighbors']['connectivities'] = early_velo.obsp['connectivities']

scv.tl.paga(early_velo, groups='EK_anno')
df = scv.get_df(early_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[768]:


late_velo.uns['neighbors']['distances'] = late_velo.obsp['distances']
late_velo.uns['neighbors']['connectivities'] = late_velo.obsp['connectivities']

scv.tl.paga(late_velo, groups='EK_anno')
df = scv.get_df(late_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[769]:


scv.pl.paga(early_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[770]:


scv.pl.paga(late_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[771]:


sc.tl.leiden(early_velo)


# In[772]:


sc.tl.leiden(late_velo)


# In[773]:


early_velo.uns['neighbors']['distances'] = early_velo.obsp['distances']
early_velo.uns['neighbors']['connectivities'] = early_velo.obsp['connectivities']

scv.tl.paga(early_velo, groups='leiden')
df = scv.get_df(early_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[774]:


late_velo.uns['neighbors']['distances'] = late_velo.obsp['distances']
late_velo.uns['neighbors']['connectivities'] = late_velo.obsp['connectivities']

scv.tl.paga(late_velo, groups='leiden')
df = scv.get_df(late_velo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[775]:


scv.pl.paga(early_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[776]:


scv.pl.paga(late_velo, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[777]:


import cellrank as cr


# In[778]:


cr.tl.terminal_states(early_velo, cluster_key="leiden", weight_connectivities=0.2, force_recompute = True, n_states = 1)


# In[786]:


cr.tl.terminal_states(late_velo, cluster_key="leiden", weight_connectivities=0.2, n_states = 3)


# In[780]:


cr.pl.terminal_states(early_velo, basis = 'draw_graph_fa')


# In[787]:


cr.pl.terminal_states(late_velo, basis = 'draw_graph_fa')


# In[790]:


cr.tl.initial_states(early_velo, cluster_key="leiden", n_states = 4)
cr.pl.initial_states(early_velo, discrete=True, basis = 'draw_graph_fa')


# In[792]:


cr.tl.initial_states(late_velo, cluster_key="leiden", n_states = 2)
cr.pl.initial_states(late_velo, discrete=True, basis = 'draw_graph_fa')


# In[793]:


cr.tl.lineages(early_velo, backward= False)
cr.pl.lineages(early_velo, same_plot=False, basis = 'draw_graph_fa')


# In[794]:


cr.tl.lineages(late_velo, backward= False)
cr.pl.lineages(late_velo, same_plot=False, basis = 'draw_graph_fa')


# In[795]:


scv.tl.recover_latent_time(  
    early_velo, root_key="initial_states_probs", end_key="terminal_states_probs"  
)


# In[796]:


scv.tl.recover_latent_time(  
    late_velo, root_key="initial_states_probs", end_key="terminal_states_probs"  
)


# In[797]:


scv.tl.paga(
    early_velo,
    groups="leiden",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="latent_time",
)


# In[798]:


scv.tl.paga(
    late_velo,
    groups="leiden",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="latent_time",
)


# In[799]:


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


# In[800]:


cr.pl.cluster_fates(
    late_velo,
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


# In[801]:


# compue DPT, starting from CellRank defined root cell
root_idx = np.where(early_velo.obs["initial_states"] == "Transitory")[0][0]
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


# In[471]:





# In[248]:


sc.pl.draw_graph(early_velo,color="leiden", layout = 'fa', legend_loc = 'on data')


# In[606]:


cr.tl.lineage_drivers(early_velo, lineages='5', use_raw=False, backward=False)


# In[257]:


cr.tl.lineage_drivers(early_velo, cluster_key = 'leiden',clusters = ['7','2','8'])


# In[256]:


cr.tl.lineage_drivers(early_velo, cluster_key = 'leiden',clusters = ['5','10','8'])


# In[607]:


cr.pl.lineage_drivers(early_velo, lineage="5", n_genes=5, basis = 'draw_graph_fa')


# In[603]:


early_velo.var_names


# In[802]:


model = cr.ul.models.GAM(early_velo)
cr.pl.gene_trends(
    early_velo,
    model=model,
    data_key="X",
    genes=["Trpm3"],
    ncols=3,
    time_key="latent_time",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 4),
    n_test_points=200
)


# In[804]:


cr.pl.heatmap(
    late_velo,
    model,
    genes=late_velo.varm['terminal_lineage_drivers']["5_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="5",
    n_jobs=1,
    backend="loky",
)


# In[807]:


late_velo


# In[809]:


cr.pl.circular_projection(late_velo, keys="EK_anno", legend_loc="right")


# In[157]:


cr.pl.circular_projection(early_velo, keys="EK_anno", legend_loc="right", save = 'figures/triangle.pdf') #do reversed to find genes


# In[178]:


import sys
get_ipython().system('{sys.executable} -m pip -q install palantir fa2')


# In[179]:


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


# In[199]:


sc.pp.pca(early_velo)
pca_projections = pd.DataFrame(early_velo.obsm["X_pca"],index=early_velo.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)
early_velo.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(early_velo,n_neighbors=30,use_rep="X_palantir")
early_velo.obsm["X_pca2d"]=early_velo.obsm["X_pca"][:,:2]
sc.tl.draw_graph(early_velo,init_pos='X_pca2d')


# In[211]:


sc.set_figure_params()
sc.pl.draw_graph(early_velo,color="EK_anno", layout = 'fa')


# In[180]:


scf.tl.tree(mg,method="ppt",Nodes=400,use_rep="palantir",
            device="cpu",seed=1,ppt_lambda=1,ppt_sigma=0.07,ppt_nsteps=100)


# In[182]:


sc.set_figure_params(figsize = [8,8])
scf.pl.graph(mg)


# In[183]:


sc.set_figure_params()
scf.tl.root(mg,395)


# In[184]:


scf.tl.pseudotime(mg,n_jobs=20,n_map=100,seed=42)


# In[285]:


scf.pl.trajectory(mg)


# In[255]:


sc.pl.draw_graph(early,color=["seg","milestones"])
#scf.tl.rename_milestones(early_velo,["bifA","Ery","DC","Root","BifB","Mono"])
# we change the color of the root milestone for better visualisations
early.uns["milestones_colors"][3]="#17bece"


# In[276]:


sc.pl.draw_graph(mg, layout = 'fa')


# In[275]:


mg


# In[ ]:


scf.tl.test_association(early,n_jobs=20)


# In[ ]:


scf.pl.test_association(early)


# In[ ]:


scf.tl.fit(early,n_jobs=20)


# In[312]:


top_genes = early_velo.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.pl.heatmap(early_velo, var_names=top_genes, sortby='latent_time', col_color='EK_anno', n_convolve=300,yticklabels=True, figsize = (8,32))


# In[313]:


early_velo.var


# In[160]:


sc.tl.leiden(early)


# In[161]:


sc.set_figure_params(figsize = [16,16])
sc.pl.draw_graph(early,color="leiden", layout = 'fa', legend_loc = 'on data')


# In[234]:


early_rgc = early[early.obs['leiden'].isin(['20','2','12','19','11','8','5','13','22','21','24','7','10','15','26','27'])]


# In[261]:


early_cones= early[early.obs['leiden'].isin(['18','16','8','5','13','22','21','24','7','10','15','26','27'])]


# In[262]:


early_mg= early[early.obs['leiden'].isin(['17','6','9','25','1','14','0','2','23','3','8','5','13','22','21','24','7','10','15','26','27'])]


# In[263]:


early_mg = early_mg[early_mg.obs['EK_anno'].isin(['Transitory','Partially muller glia','Mature muller glia'])]


# In[264]:


sc.set_figure_params()
sc.pl.draw_graph(early_mg,color="leiden", layout = 'fa', legend_loc = 'on data')


# In[265]:


sc.set_figure_params()
sc.pl.draw_graph(early_cones,color="leiden", layout = 'fa', legend_loc = 'on data')


# In[235]:


sc.set_figure_params()
sc.pl.draw_graph(early_rgc,color="leiden", layout = 'fa', legend_loc = 'on data')


# In[236]:


sc.set_figure_params()
sc.pl.draw_graph(early_rgc,color="EK_anno", layout = 'fa', legend_loc = 'on data')


# In[266]:


sc.set_figure_params()
sc.pl.draw_graph(early_mg,color="EK_anno", layout = 'fa', legend_loc = 'on data')


# In[267]:


sc.set_figure_params()
sc.pl.draw_graph(early_cones,color="EK_anno", layout = 'fa', legend_loc = 'on data')


# In[237]:


early_rgc = early_rgc[early_rgc.obs['EK_anno'].isin(['RGC','Transitory','Partially muller glia','Mature muller glia'])]


# In[238]:


sc.set_figure_params()
sc.pl.draw_graph(early_rgc,color="EK_anno", layout = 'fa', legend_loc = 'on data')


# In[198]:


import scFates as scf


# In[199]:


early_mg.var


# In[200]:


sc.pl.pca(early_mg,color="Clu",cmap="RdBu_r", use_raw = False)


# In[241]:


early_rgc.var_names


# In[242]:


sc.pl.pca(early_rgc,color="Clu",cmap="RdBu_r", use_raw = False)


# In[268]:


sc.pp.pca(early_mg)


# In[269]:


sc.pp.pca(early_cones)


# In[243]:


sc.pp.pca(early_rgc)


# In[270]:


scf.tl.curve(early_mg,Nodes=30,use_rep="X_draw_graph_fa",ndims_rep=2,)


# In[271]:


scf.tl.curve(early_cones,Nodes=30,use_rep="X_draw_graph_fa",ndims_rep=2,)


# In[244]:


scf.tl.curve(early_rgc,Nodes=30,use_rep="X_draw_graph_fa",ndims_rep=2,)


# In[302]:


scf.pl.graph(early_mg,basis="draw_graph_fa", save = 'mg1.pdf')


# In[310]:


scf.pl.graph(early_cones,basis="draw_graph_fa", save = 'mg2.pdf')


# In[311]:


scf.pl.graph(early_rgc,basis="draw_graph_fa", save = 'mg3.pdf')


# In[274]:


scf.tl.root(early_mg, root = 2)


# In[275]:


scf.tl.root(early_cones, root = 2)


# In[246]:


scf.tl.root(early_rgc, root = 2)


# In[276]:


scf.tl.pseudotime(early_mg,seed=42)


# In[277]:


scf.tl.pseudotime(early_cones,seed=42)


# In[247]:


scf.tl.pseudotime(early_rgc,seed=42)


# In[301]:


scf.pl.trajectory(early_mg,basis="draw_graph_fa",arrows=True,arrow_offset=3, color_cells = 'EK_anno', alpha = 1, scale_path = 0.5, save = 'mg.pdf')


# In[279]:


scf.pl.trajectory(early_cones,basis="draw_graph_fa",arrows=True,arrow_offset=3, color_cells = 'EK_anno', alpha = 1, scale_path = 0.5)


# In[248]:


scf.pl.trajectory(early_rgc,basis="draw_graph_fa",arrows=True,arrow_offset=3, color_cells = 'EK_anno', alpha = 1, scale_path = 0.5)


# In[207]:


early_mg.obs


# In[280]:


scf.tl.linearity_deviation(early_mg,start_milestone="2",end_milestone="0",n_jobs=20,plot=True,basis="draw_graph_fa")


# In[281]:


scf.tl.linearity_deviation(early_cones,start_milestone="2",end_milestone="0",n_jobs=20,plot=True,basis="draw_graph_fa")


# In[249]:


scf.tl.linearity_deviation(early_rgc,start_milestone="2",end_milestone="0",n_jobs=20,plot=True,basis="draw_graph_fa")


# In[282]:


scf.pl.linearity_deviation(early_mg,start_milestone="2",end_milestone="0")


# In[283]:


scf.pl.linearity_deviation(early_cones,start_milestone="2",end_milestone="0")


# In[250]:


scf.pl.linearity_deviation(early_rgc,start_milestone="2",end_milestone="0")


# In[210]:


import sys
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout


# In[284]:


sc.pl.pca(early_mg,color=["Igfbpl1"], use_raw = False)


# In[285]:


sc.pl.pca(early_mg,color=["EK_anno"], use_raw = False)


# In[286]:


sc.pl.pca(early_cones,color=["EK_anno"], use_raw = False)


# In[252]:


sc.pl.pca(early_rgc,color=["EK_anno"], use_raw = False)


# In[287]:


scf.tl.test_association(early_mg,n_jobs=3)


# In[288]:


scf.tl.test_association(early_cones,n_jobs=3)


# In[253]:


scf.tl.test_association(early_rgc,n_jobs=3)


# In[289]:


scf.tl.fit(early_mg,n_jobs=3)


# In[290]:


scf.tl.fit(early_cones,n_jobs=3)


# In[254]:


scf.tl.fit(early_rgc,n_jobs=3)


# In[255]:


early_rgc.var


# In[215]:


early_mg.var


# In[216]:


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


# In[291]:


scf.pl.single_trend(early_mg,"Mlc1",basis="draw_graph_fa",color_exp="k", use_raw = False)


# In[256]:


scf.pl.single_trend(early_rgc,"Hes1",basis="draw_graph_fa",color_exp="k", use_raw = False)


# In[292]:


scf.pl.single_trend(early_mg,"Igfbpl1",basis="draw_graph_fa",color_exp="k", use_raw = False)


# In[315]:


scf.pl.single_trend(early_cones,"Pde6g",basis="draw_graph_fa",color_exp="k", use_raw = False, save = 'pde.pdf')


# In[316]:


scf.pl.single_trend(early_rgc,"C1ql1",basis="draw_graph_fa",color_exp="k", use_raw = False, save = 'c1q.pdf')


# In[294]:


scf.tl.cluster(early_mg, n_neighbors = 20, metric="correlation")


# In[295]:


scf.tl.cluster(early_cones, n_neighbors = 20, metric="correlation")


# In[257]:


scf.tl.cluster(early_rgc, n_neighbors = 20, metric="correlation")


# In[296]:


early_mg.var.cluters.unique()


# In[297]:


early_cones.var.cluters.unique()


# In[258]:


early_rgc.var.cluters.unique()


# In[308]:


'mgmg' + c + '.pdf'


# In[309]:


for c in early_mg.var["cluters"].unique():
    filename = 'mgmg' + c + '.pdf'
    scf.pl.trends(early_mg,features=early_mg.var_names[early_mg.var.cluters==c],basis="draw_graph_fa", save = filename)


# In[312]:


for c in early_cones.var["cluters"].unique():
    filename = 'mgc' + c + '.pdf'
    scf.pl.trends(early_cones,features=early_cones.var_names[early_cones.var.cluters==c],basis="draw_graph_fa", save = filename)


# In[313]:


for c in early_rgc.var["cluters"].unique():
    filename = 'mgr' + c + '.pdf'
    scf.pl.trends(early_rgc,features=early_rgc.var_names[early_rgc.var.cluters==c],basis="draw_graph_fa", save = filename)


# In[ ]:




