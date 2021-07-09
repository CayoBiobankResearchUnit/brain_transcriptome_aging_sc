#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/',dataset,'_doublets_removed.rds'))

# saveRDS(as.integer(which(with(within(as.data.frame(colData(cds)),{ cluster = factor(clusters(cds),levels=sort(unique(clusters(cds)))) }),tapply(scrublet_score,cluster,mean)) > 0.165)),file=paste0('checkpoints/doublet_clusters_',dataset,'.rds'))

doublet.clusters = readRDS(paste0('checkpoints/doublet_clusters_',dataset,'.rds'))

cds = cds[,!clusters(cds) %in% doublet.clusters] %>% { reduce_dim_and_cluster(., k=cluster.k.final) }

write(paste0(nrow(colData(cds)),' cells remain after removing ',length(doublet.clusters),' doublet clusters.'),file=paste0('reports/',dataset,'_doublet_cluster_qc.txt'))
write(paste0('Clustering with k = ',cluster.k.final),file=paste0('reports/',dataset,'_doublet_cluster_qc.txt'),append=TRUE)

saveRDS(cds,file=paste0('checkpoints/',dataset,'_doublet_clusters_removed.rds'))
