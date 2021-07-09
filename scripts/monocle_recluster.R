#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/',dataset,'_doublet_clusters_removed.rds'))
cds$monocle_cluster = clusters(cds)
cds$curated_cluster = as.integer(NA)

# Vector of clusters to recluster
lumped.clusters = readRDS(paste0('checkpoints/lumped_clusters_',dataset,'.rds'))

lumped.cluster.results = list()
for (i in 1:length(lumped.clusters)) {
	x = lumped.clusters[i]
	this = cds[,clusters(cds) == x] %>% { cluster_cells(.,cluster_method='leiden',k=cluster.k.recluster) }
	colData(cds)[rownames(colData(this)),]$curated_cluster = max(c(0,cds$curated_cluster),na.rm=TRUE) + as.integer(clusters(this))
}

# Add back unaffected clusters and assign them new cluster numbers
colData(cds)$curated_cluster[is.na(colData(cds)$curated_cluster)] = max(c(0,cds$curated_cluster),na.rm=TRUE) + as.integer(colData(cds)$monocle_cluster[is.na(colData(cds)$curated_cluster)][,drop=TRUE])

# Make into a factor based on decreasing cell counts
colData(cds)$curated_cluster = factor(colData(cds)$curated_cluster,levels=names(sort(table(colData(cds)$curated_cluster),decreasing=TRUE)))

# Cluster names are arbitrary
levels(colData(cds)$curated_cluster) = 1:nlevels(colData(cds)$curated_cluster)

out = colData(cds)$curated_cluster
names(out) = rownames(colData(cds))

cds@clusters[['UMAP']]$clusters = out

saveRDS(cds,file=paste0('checkpoints/',dataset,'_reclustered.rds'))
