###################################################
# Workflow for Enhanced Experimental Signal (EES) #
#                                                 #
# Much of this code came from or is based on      #
# code from Savell et al. 2020 (Sci. Adv.) and    #
# this Krishnaswamy Lab vignette:                 #
# https://nbviewer.org/github/KrishnaswamyLab/    #
# MELD/blob/main/notebooks/Wagner2018_Chordin_    #
# Cas9_Mutagenesis.ipynb                          #
#                                                 #
# See README for instructions on how to run this  #
# code in x86 Python and how to generate the      #
# input files.                                    #
###################################################

# import packages
import pandas as pd
import numpy as np
import graphtools as gt
import magic
import scprep
import meld
import matplotlib.pyplot as plt
import phate
# Required only for differential expression
import diffxpy.api as de

# Set defaults for matplotlib font sizes
plt.rc('font', size = 12)

# Make sure plots & clusters are reproducible
np.random.seed(638)

# Load CSV file containing counts data
# Cells are on columns here
data = scprep.io.load_csv('allRats_counts_min10.csv', sparse = False, cell_axis = 'column')

# Load CSV file containing metadata
metadata = scprep.io.load_csv('allRats_min10_metadata.csv')
 
# Create indicator array for samples
# This is the "raw experimental signal" (RES)
# Sets lacz = -1, reln = 1
metadata['RES'] = np.array([-1 if label.startswith('lacz') else 1 for label in metadata['target']])

# Map cluster names to cluster IDs
# Check your metadata first and change to your own links between seurat_clusters and cluster name
def cluster (row):
   if row['seurat_clusters'] == 0 :
      return 'Drd1'
   if row['seurat_clusters'] == 1 :
      return 'Drd2-1'
   if row['seurat_clusters'] == 2 :
      return 'Olig'
   if row['seurat_clusters'] == 3 :
      return 'Astro'  
   if row['seurat_clusters'] == 4 :
      return 'Grm8'
   if row['seurat_clusters'] == 5 :
      return 'Micro'
   if row['seurat_clusters'] == 6 :
      return 'GABA'
   if row['seurat_clusters'] == 7 :
      return 'Poly' 
   if row['seurat_clusters'] == 8 :
      return 'Drd2-2'
   if row['seurat_clusters'] == 9 :
      return 'Pvalb'
   if row['seurat_clusters'] == 10 :
      return 'Mural'
   if row['seurat_clusters'] == 11 :
      return 'Drd3'  
   if row['seurat_clusters'] == 12 :
      return 'Sst'
   return 'Other'

metadata['cluster'] = metadata.apply (lambda row: cluster(row), axis = 1)

# Sort clusters by abundance
cluster_abundance = metadata.groupby('cluster').aggregate({'RES':np.mean}).sort_values('RES')
cluster_abundance['newClusterID'] = np.arange(cluster_abundance.shape[0])

# Relabel existing clusters in metadata
new_clusters = cluster_abundance.loc[metadata['cluster']]['newClusterID']
new_clusters.index = metadata.index
metadata['clusterID'] = new_clusters

## ------------- TRANSFORM DATA IN MAGIC & GENERATE EES ---------------------------------------------

# Create primary graph
G = gt.Graph(data, knn = 9, decay = 10, n_pca = 100, use_pygsp = True, n_jobs = -2, verbose = True)

# Impute gene expression using MAGIC
# (Not entirely sure what this is doing...)
magic_op = magic.MAGIC(knn = G.knn, decay = G.decay)
data_magic = magic_op.fit_transform(data, graph = G)

# Calculate densities using main graph and MAGIC-transformed data
meld_op = meld.MELD()
reln_densities = meld_op.fit_transform(G, metadata['RES'])

# Calculate relative likelihoods (EES values) from densities
reln_likelihoods = meld.utils.normalize_densities(reln_densities)

# Mean-center likelihoods
reln_likelihoods = reln_likelihoods - np.mean(reln_likelihoods)

# Add Reln likelihoods column (1) to metadata
metadata['EES'] = reln_likelihoods[1]

# Save new metadata with EES values to CSV file
metadata.to_csv("metadata_with_EES.csv", header = True)

## ------------- GENERATE DOTPLOT FOR CLUSTER EES VALUES --------------------------------------------

# RELN treatment = orange
# LACZ treatment = blue
samples_cdict = {'reln' : '#fb6a4a', 'lacz' : '#08519c'}

samples_cvec = np.array([samples_cdict[s] for s in metadata['target']])

def plot_EES_in_clusters(ax, clusters, EES, c = None):
    
    clusters_idx = np.arange(len(set(clusters)))
    n_clusts = len(set(clusters_idx))
    
    # Calculate means
    means = np.zeros(n_clusts)
    for i, cl in enumerate(np.unique(clusters)):
        means[i] = np.mean(EES[clusters == cl])
    
    # Plotting cells
    x = clusters + np.random.normal(0, 0.1, len(clusters))
    y = EES
    r = np.random.choice(len(y), len(y), replace = False)
    ax.scatter(x[r], y[r], c = c[r], s = 2)
    
    # Plotting means
    ax.scatter(np.arange(n_clusts), means, c = '#cbc9ff', edgecolors = 'black', lw = 1.5, marker = 'o', zorder = 3, s = 100)
    
    # Plotting vetical lines
    for i in np.unique(clusters):
        ax.axvline(i, c = 'k', lw = 0.1, zorder = 0)

fig, ax = plt.subplots(1, figsize = (14,10))
plot_EES_in_clusters(ax, clusters = metadata['clusterID'], EES = metadata['EES'], c = samples_cvec)

cluster_labels = metadata.groupby('cluster').aggregate({'RES':np.mean}).sort_values('RES').index.values

ax.set_xticks(np.arange(len(cluster_labels)))
ax.set_xticklabels(cluster_labels, rotation = 45, ha = 'right', fontsize = 14)

scprep.plot.utils.shift_ticklabels(ax.xaxis, dx = 0.1)

ax.set_ylabel('Enhanced Experimental Signal', fontsize = 18)
ax.tick_params(labelsize = 14)
fig.tight_layout()

# To show in plot window
plt.show()

## ------------- PEARSON CORRELATION BETWEEN GENE EXPRESSION & EES ----------------------------------

# Build EES gene correlation dataframe and tables. Calculate and include Log2FC in tables.
# Note: Will give an expected divide by zero error

# Disable display of divide by zero warning
import warnings;
warnings.filterwarnings('ignore');

# List of all clusters
cluster_list = ['Astro','Drd1','Drd2-1','Drd2-2','Drd3','GABA','Grm8',
                'Micro','Mural','Olig','Poly','Pvalb','Sst']

# Do Pearson correlation for each cluster and save data tables
for j in cluster_list:
    cluster_name = j
    
    # Un-normalize data (nonlogp1)
    is_in_cluster = metadata['cluster'] == cluster_name
    data_cluster, metadata_cluster = scprep.select.select_rows(data, metadata, idx = is_in_cluster)
    data_cluster = (np.e ** data_cluster) -1
    tx_name = 'lacz'
    is_in_tx = metadata_cluster['target'] == tx_name
    data_tx, metadata_tx = scprep.select.select_rows(data_cluster, metadata_cluster, idx = is_in_tx)
    tx_name = 'reln'
    is_in_tx = metadata_cluster['target'] == tx_name
    data_tx2, metadata_tx2 = scprep.select.select_rows(data_cluster, metadata_cluster, idx = is_in_tx)
    
    # Calculate fold change (non log plus one)
    mean_lacz = data_tx.mean(axis = 0)
    mean_reln = data_tx2.mean(axis = 0)
    FC = mean_reln / mean_lacz
    FC_log = np.log2(FC)
    FC_log = FC_log.replace([np.inf, -np.inf], np.nan)
    FC_log = FC_log.dropna()
    fc_genes = FC_log.reset_index()
    fc_genes.columns = ['gene','logFC']
    
    # Calculate Pearson's r
    is_in_cluster = metadata['cluster'] == cluster_name
    data_cluster, metadata_cluster = scprep.select.select_rows(data, metadata, idx = is_in_cluster)
    
    gene_list = fc_genes['gene']
    gene_list = np.array(gene_list)
    
    fc_list = fc_genes['logFC']
    fc_list = np.array(fc_list)
    
    dremi_list = []
    for i in gene_list:
        value = scprep.stats.pairwise_correlation(metadata_cluster['EES'], data_cluster[i])
        dremi_list.append(value)
        
    # Rank by Pearson and build dataframe
    dremi_cluster = pd.DataFrame({'Gene':gene_list, 'r_val':dremi_list, 'log2FC':fc_list})
    dremi_cluster = dremi_cluster.sort_values(by = ['r_val'], ascending = False)
    dremi_cluster['r_Rank'] = [1 * i for i in range(len(dremi_cluster))]
    dremi_cluster['r_val'] = dremi_cluster['r_val'].str[0] #remove brackets
    dremi_cluster['r_val'] = dremi_cluster['r_val'].str[0]
    
    # Save cluster data
    dremi_cluster.to_csv('Reln_'+cluster_name+'_Pearson_FC_by_Gene_no_cutoff.csv', index = False, header = True)

## ------------- SUBCLUSTERING WITH VERTEX FREQUENCY CLUSTERING (VFC) -------------------------------

# Apply PHATE to data
phate_op = phate.PHATE(knn = 10, decay = 10, n_jobs = -1)
data_phate = phate_op.fit_transform(data)

# Get cluster indexes and number of cells per cluster
clusters, counts = np.unique(metadata['clusterID'], return_counts=True)

# Keep cluster labels with at least 1% of the data
clusters = clusters[counts > data.shape[0] * 0.01]

data_cluster_phate = {}

for clust in clusters:
    curr_data = data.loc[metadata['clusterID'] == clust]
    data_cluster_phate[clust] = phate.PHATE(verbose = 0).fit_transform(curr_data)

# Visualizing EES values for each cluster using PHATE
fig, axes = plt.subplots(3,4, figsize = (4*3, 5*3))

# Play with vmin & vmax for best color scaling
for i , ax in enumerate(axes.flatten()):
    if not i < len(clusters):
        ax.axis('off')
        continue
    curr_cluster = clusters[i]
    curr_phate = data_cluster_phate[curr_cluster]
    
    scprep.plot.scatter2d(curr_phate, 
                          c = metadata['EES'].loc[metadata['clusterID'] == curr_cluster], 
                          cmap = meld.get_meld_cmap(), vmin = -0.4, vmax = 0.2,
                         ax = ax, ticks = False, 
                          title = 'Cluster {} ({})'.format(curr_cluster, curr_phate.shape[0]), 
                          legend = False, fontsize = 10, filename = "subplots_n0.4_0.2.png", dpi = 100)

# Build a VFC operator for each cluster
np.random.seed(0)
vfc_op_per_cluster = {}

for clust in np.unique(metadata['clusterID']):
    curr_G = gt.Graph(data.loc[metadata['clusterID'] == clust], use_pygsp=True)
    curr_G.compute_fourier_basis()
    curr_sample_labels = metadata['RES'].loc[metadata['clusterID'] == clust]
    curr_likelihood = metadata['EES'].loc[metadata['clusterID'] == clust]
    curr_vfc = meld.VertexFrequencyCluster(n_clusters = 3)
    curr_vfc.fit_transform(curr_G, curr_sample_labels, curr_likelihood)
    vfc_op_per_cluster[clust] = curr_vfc

# Build subcluster dictionary
subclustering_results = {}
for clust in np.unique(metadata['clusterID']):
    curr_vfc = vfc_op_per_cluster[clust]
    clusters_by_n = {}
    for n in [2,3,4,5]:
        clusters_by_n[n] = curr_vfc.predict(n)
    subclustering_results[clust] = clusters_by_n

# Which boxes to highlight in big plot
# Cluster:box (1 = EES; 2-5 = k of 2-5 respectively)
picked_clusters = {2:2, 3:2, 4:4, 5:3, 10:3, 11:4}

# Big plot
# For each cluster: EES, k = 2, 3, 4, 5
fig, axes= plt.subplots(10,5, figsize=(10, 10*3))

for i, r_ax in enumerate(axes):
    clust = clusters[i]
    curr_phate = data_cluster_phate[clust]
    
    for i, n in enumerate([0,2,3,4,5]):
        if i == 0:
            cvec = metadata['EES'].loc[metadata['clusterID'] == clust]
            ylabel= str(clust) + ' ({})'.format(len(cvec))
            
            cmap= meld.get_meld_cmap()
            vmin = -0.1
            vmax = 0.2           
        else:
            cvec = subclustering_results[clust][n]
            cmap = ylabel = vmin = vmax = None
            
        if np.sum(metadata['clusterID'] == clust) < (metadata.shape[0] * 0.01):
            plt.setp(r_ax[i].spines.values(), color = 'lightgrey', linewidth = 4)
        scprep.plot.scatter2d(curr_phate, c = cvec, cmap = cmap, vmin = vmin, vmax = vmax,
                          ax = r_ax[i], ticks = False, ylabel = ylabel, fontsize = 8, legend = False,
                          filename = "subplots_kmeans.png", dpi = 100)
        # Red outline for clusters that we decided to highlight
        if clust in picked_clusters:
            if picked_clusters[clust] == n:
                plt.setp(r_ax[i].spines.values(), color = 'r', linewidth = 4)

# Add VFC number to metadata
metadata['VFC'] = metadata['clusterID'].copy()
for curr_cluster in clusters:
    if curr_cluster not in picked_clusters:
        continue
    new_clusters = subclustering_results[curr_cluster][picked_clusters[curr_cluster]]
    metadata.loc[metadata['clusterID'] == curr_cluster, 'VFC'] = new_clusters + 1 + np.max(metadata['VFC'])

# Uncomment this only if you want subcluster jitter plot sorted by mean EES instead of grouped by cluster
# Omit this step to keep jitter plot sorted by sub/cluster
# metadata['VFC'] = scprep.utils.sort_clusters_by_values(metadata['VFC'], metadata['EES'])

# Set colors for jitter plot
sample_cmap = {-1: "#3182bd", 1: "#de2d26"}

# Jitter plot separated by subcluster
ax = scprep.plot.jitter(metadata['VFC'], metadata['EES'], c = metadata['RES'],
                        cmap = sample_cmap, figsize = (18,7), legend = False,
                        filename = "subcluster_jitter.pdf")

## ------------- CLOSER LOOK AT DRD1 CLUSTER --------------------------------------------------------

# Subset data from this cluster (Drd1 = cluster 4)
curr_cluster = 4
curr_cluster_data = data.loc[np.isin(metadata['clusterID'], curr_cluster)].copy()
curr_cluster_data_magic = data_magic.loc[np.isin(metadata['clusterID'], curr_cluster)].copy()
curr_cluster_metadata = metadata.loc[np.isin(metadata['clusterID'], curr_cluster)].copy()
curr_data_phate = data_cluster_phate[curr_cluster]

# Plot D1 using Phate
scprep.plot.scatter2d(curr_data_phate, c = curr_cluster_metadata['EES'], cmap = meld.get_meld_cmap(), 
                      figsize = (6,5), s = 5, label_prefix = "PHATE", ticks = False,
                      filename = "D1_phate.png", dpi = 100)

# Vertex frequency clustering (VFC) plot for D1
scprep.plot.scatter2d(curr_data_phate, c = curr_cluster_metadata['VFC'], 
                      figsize = (5,5), s = 5, label_prefix = "PHATE", ticks = False,
                      filename = "D1_VFC_plot.png", dpi = 100)

# Jitter plot for D1 subclusters
scprep.plot.jitter(curr_cluster_metadata['VFC'], curr_cluster_metadata['EES'], 
                   c = curr_cluster_metadata['RES'], cmap = sample_cmap,
                  legend_anchor = (1,1), filename = "D1_jitter.png", dpi = 100)

# Import "marker genes" (DEGs for D1 from DESeq2 analysis)
marker_genes = pd.read_csv("marker_genes.csv")

## Differential expression for each D1 subcluster vs. all other D1 subclusters
subcluster_marker_genes = {}

# Iterate over VFC subclusters
for curr_subcluster in np.unique(curr_cluster_metadata['VFC']):
    # Create mask for current subcluster
    is_curr_subcluster = curr_cluster_metadata['VFC'] == curr_subcluster
    
    # Rank test comparing each subcluster to all other cells in cluster
    curr_de_results = de.test.two_sample(curr_cluster_data.values, 
                                         gene_names = data.columns, 
                                         grouping = is_curr_subcluster,
                                         test = 'rank').summary()
    
    # Keep only top 1% of DEGs (does this do that??)
    top_genes = curr_de_results.sort_values('qval').loc[curr_de_results['qval'] < 0.05, 'gene'].values

    # Get list of DEGs that are also on marker gene list from DESeq2
    curr_marker_genes = []
    for gene in top_genes:
        if gene in marker_genes['gene'].values:
            curr_marker_genes.append(gene)
    subcluster_marker_genes[curr_subcluster] = curr_marker_genes

# Build dot plot for DEGs
all_markers = np.unique(np.hstack(subcluster_marker_genes.values()))
curr_marker_genes = marker_genes.loc[np.isin(marker_genes['gene'], all_markers), ['gene']]
curr_marker_genes_list = curr_marker_genes['gene'].values.tolist()

scprep.plot.marker_plot(data = curr_cluster_data.values,
                        gene_names = curr_cluster_data.columns,
                        clusters = curr_cluster_metadata['VFC'], 
                        markers = curr_marker_genes_list,
                        cmap = "PuRd",
                        figsize = (12,8))
