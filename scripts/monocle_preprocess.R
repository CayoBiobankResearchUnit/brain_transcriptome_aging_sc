#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/',dataset,'_doublets_included.rds'))

# Read in manually set thresholds
doublet.thresholds = readRDS(paste0('checkpoints/doublet_thresholds_',dataset,'.rds'))

write(c('Thresholding samples by the following doublet scores:',apply(format(doublet.thresholds),1,function(x) paste(x,collapse='\t')),'\n'),file=paste0('reports/',dataset,'_doublet_score_qc.txt'))

cds = mark_doublets(cds,doublet.thresholds) %>% remove_doublets %>% { reduce_dim_and_cluster(., k=cluster.k.initial) }

write(paste0(nrow(colData(cds)),' cells pass doublet thresholds.'),file=paste0('reports/',dataset,'_doublet_score_qc.txt'),append=TRUE)
write(paste0('Clustering with k = ',cluster.k.initial),file=paste0('reports/',dataset,'_doublet_score_qc.txt'),append=TRUE)

saveRDS(cds,file=paste0('checkpoints/',dataset,'_doublets_removed.rds'))
