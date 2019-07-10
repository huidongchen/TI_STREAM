#!/usr/bin/env python

import h5py
import sys
import scipy.sparse
import sklearn.decomposition
import pandas as pd
import numpy as np
import os

import json
import stream as st
import yaml
import networkx as nx


import time
checkpoints = {}

##### Read data #####

# parse location of dataset and output folder
dataset_location = sys.argv[1]
# dataset_location = './bifurcating.h5'
output_folder = sys.argv[2]
# output_folder = './output/'


# read in sparse matrix
dataset_h5 = h5py.File(dataset_location)
expression_h5 = dataset_h5["data"]["expression"]
counts_h5 = dataset_h5["data"]["counts"]
expression = scipy.sparse.csc_matrix((
  expression_h5["x"][()],
  expression_h5["i"][()],
  expression_h5["p"][()]),
  expression_h5["dims"][()]
)
counts = scipy.sparse.csc_matrix((
  counts_h5["x"][()],
  counts_h5["i"][()],
  counts_h5["p"][()]),
  counts_h5["dims"][()]
)
cell_ids = expression_h5["rownames"][()].astype(str)
gene_ids = expression_h5["colnames"][()].astype(str)


# Infer trajectories

# read in parameters
definition = open('./definition.yml', 'r')
task = yaml.safe_load(definition)
p = dict()
for x in task["parameters"]:
    p[x['id']] = x['default']

pd.DataFrame(counts.toarray(),index=cell_ids,columns=gene_ids).T.to_csv('./counts.tsv',sep='\t')

checkpoints["method_afterpreproc"] = time.time()


adata=st.read(file_name="./counts.tsv")
st.add_cell_labels(adata)
st.add_cell_colors(adata)

if(p["norm"]):
    st.normalize_per_cell(adata)
if(p["log2"]):
    st.log_transform(adata)


st.filter_genes(adata,min_num_cells = max(5,int(round(adata.shape[0]*0.001))),min_pct_cells = None,expr_cutoff = 1)
if(adata.shape[1]<1000):
    adata.uns['var_genes'] = gene_ids
    adata.obsm['var_genes'] = adata.X
else:
    st.select_variable_genes(adata,loess_frac=p["loess_frac"])

st.dimension_reduction(adata,nb_pct = p["nb_pct"],n_components = p["n_components"],n_jobs = p["n_jobs"],method = p["method"])
st.plot_dimension_reduction(adata,n_components = p["n_components"],save_fig=p["save_fig"])
st.plot_visualization_2D(adata,save_fig=p["save_fig"],nb_pct=p["nb_pct"])


st.seed_elastic_principal_graph(adata,clustering=p["clustering"],n_clusters=p["n_clusters"],nb_pct=p["nb_pct"])
st.elastic_principal_graph(adata,epg_alpha=p["epg_alpha"],save_fig=p["save_fig"])

if(len(adata.obs['branch_id'].unique())>1):
    if(not p["disable_EPG_optimize"]):
        st.optimize_branching(adata)

if(not p["disable_EPG_ext"]):
    st.extend_elastic_principal_graph(adata)

st.plot_branches_with_cells(adata,n_components=p["n_components"])
st.plot_visualization_2D(adata,color_by='branch',save_fig=p["save_fig"])


st.subwaymap_plot(adata,root=p["root"],fig_size=(8,6), save_fig=p["save_fig"]) 
st.stream_plot(adata,root=p["root"],fig_size=(8,8),save_fig=p["save_fig"])



checkpoints["method_aftermethod"] = time.time()


# grouping
grouping = pd.DataFrame({"cell_id": adata.obs.index, "group_id": adata.obs.node})

# milestone network
milestone_network = pd.DataFrame(columns=['from','to','length','directed'])
for i, edge in enumerate(adata.uns['epg'].edges):
    dict_nodes_pos = nx.get_node_attributes(adata.uns['epg'],'pos')
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


# Save output -------------------------------------------------------------
milestone_network.to_csv(output_folder + "milestone_network.csv", index = False)
progressions.to_csv(output_folder + "progressions.csv", index = False)