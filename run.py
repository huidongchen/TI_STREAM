#!/usr/local/bin/python

import dynclipy
task = dynclipy.main()
# task = dynclipy.main(
#   ["--dataset", "/code/example.h5", "--output", "/mnt/output"],
#   "/code/definition.yml"
# )

import stream as st
import os
import sys
import json
import pandas as pd

import time
checkpoints = {}


#   ____________________________________________________________________________
#   Load data                                                               ####
task["counts"].to_csv("/counts.csv")

p = task["parameters"]

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Create trajectory                                                       ####
# normalise data

adata=st.read(file_name="/counts.csv")
st.add_cell_labels(adata)
st.add_cell_colors(adata)

if(p["norm"]):
	st.normalize_per_cell(adata)
if(p["log2"]):
	st.log_transform(adata)

st.filter_genes(adata,min_num_cells = max(5,int(round(adata.shape[0]*0.001))),min_pct_cells = None,expr_cutoff = 1)
st.select_variable_genes(adata,loess_frac=p["loess_frac"])

st.dimension_reduction(adata,nb_pct = p["nb_pct"],n_components = p["n_components"],n_jobs = p["n_jobs"],method = p["method"])
st.plot_dimension_reduction(adata,n_components = p["n_components"],save_fig=True)
st.plot_visualization_2D(adata,save_fig=True,nb_pct=p["nb_pct"])

st.seed_elastic_principal_graph(adata,clustering=p["clustering"],n_clusters=p["clustering"],nb_pct=p["nb_pct"])
st.elastic_principal_graph(adata,epg_alpha=p["epg_alpha"],save_fig=True)

if(not p["disable_EPG_optimize"]):
	st.optimize_branching(adata)

if(not p["disable_EPG_ext"]):
	st.extend_elastic_principal_graph(adata)

st.plot_visualization_2D(adata,color_by='branch',save_fig=True)

st.subwaymap_plot(adata,root=p["root"],fig_size=(8,6), save_fig=True) 
st.stream_plot(adata,root=p["root"],fig_size=(8,8),save_fig=True)

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process output & save                                                   ####

# grouping
grouping = pd.DataFrame({"cell_id": adata.obs.index, "group_id": adata.obs.node})

# milestone network
for i, edge in enumerate(adata.uns['epg'].edges):
    milestone_network.loc[i] = ""
    milestone_network.loc[i]['from'] = edge[0]
    milestone_network.loc[i]['to'] = edge[1]
    milestone_network.loc[i]['length'] = np.sqrt(((dict_nodes_pos[edge[0]] - dict_nodes_pos[edge[1]])**2).sum())
    milestone_network.loc[i]['directed'] = False

# dimred
dimred = pd.DataFrame([x for x in adata.obsm['X_vis_umap'].T]).T
dimred.columns = ["comp_" + str(i+1) for i in range(dimred.shape[1])]
dimred["cell_id"] = adata.obs.index

# progressions
progressions = adata.obs[p["root"]+'_pseudotime'].reset_index()
progressions.columns = ["cell_id", "pseudotime"]

# save
dataset = dynclipy.wrap_data(cell_ids = adata.obs.index)
dataset.add_trajectory(
  grouping = grouping,
  milestone_network = milestone_network,
  progressions = progressions,
)
dataset.add_dimred(dimred = dimred)
dataset.add_timings(timings = checkpoints)
dataset.write_output(task["output"])
