#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/',dataset,'_reclustered.rds'))

cell.assignments = readRDS(paste0('checkpoints/cell_assignments_',dataset,'.rds'))

cds$assigned_cell_type = 'Unknown'
for (i in 1:nrow(cell.assignments)) {
	j = eval(parse(text=paste0('c(',cell.assignments$clusters[i],')')))
	cds$assigned_cell_type[clusters(cds) %in% j] = cell.assignments$cell_type[i]
}

# Merge metadata file
cds = merge_metadata(cds,metadata.file=paste0('data/metadata_',dataset,'.txt'))

cds$assigned_cell_type = factor(cds$assigned_cell_type,levels=cell.levels)

assigned.cell.clusters = unlist(lapply(cell.levels,function(i) {
	this = cds[,cds$assigned_cell_type == i]
	out = if (nlevels(this$curated_cluster[,drop=TRUE]) > 1) paste(this$assigned_cell_type,as.integer(this$curated_cluster[,drop=TRUE])) else as.character(this$assigned_cell_type)
	names(out) = rownames(colData(this))
	out
}))

cds$assigned_cell_cluster = factor(assigned.cell.clusters[rownames(colData(cds))],levels=unlist(lapply(cell.levels,function(i) {
	this = cds[,cds$assigned_cell_type == i]
	if (nlevels(this$curated_cluster[,drop=TRUE]) > 1) paste(i,1:nlevels(this$curated_cluster[,drop=TRUE])) else as.character(i)
})))

saveRDS(cds,file=paste0('checkpoints/',dataset,'_classified_manual.rds'))
saveRDS(rownames(rowData(cds)),file=paste0('checkpoints/',dataset,'_mmul_genes.rds'))
