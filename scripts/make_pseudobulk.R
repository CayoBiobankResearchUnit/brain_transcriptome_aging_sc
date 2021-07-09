#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/',dataset,'_classified_manual.rds'))

meta = read.delim(paste0('data/metadata_',dataset,'.txt'),stringsAsFactors=FALSE)

ids = sort(unique(cds$sample))

meta = subset(meta,sample %in% ids)

meta = subset(meta,age >= min.age)
ids = meta$sample

rownames(meta) = NULL

pseudobulk = lapply(cell.levels,function(x) {
	this = cds[,cds$assigned_cell_type == x]
	out = do.call(cbind,lapply(ids,function(i) {
		rowSums(matrix(assays(this)$counts[,which(colData(this)$sample == i)],nrow=nrow(assays(this)$counts),dimnames=list(rownames(assays(this)$counts),colnames(assays(this)$counts[,which(colData(this)$sample == i)]))))
	}))
	colnames(out) = ids
	out
})

keep.genes = lapply(pseudobulk,function(x) {
	tpm = x / colSums(x) * 1e6
	names(which(rowMeans(tpm) >= tpm.cutoff))
})

names(pseudobulk) = names(keep.genes) = cell.levels

saveRDS(pseudobulk,file=paste0('checkpoints/',dataset,'_pseudobulk.rds'))
saveRDS(keep.genes,file=paste0('checkpoints/',dataset,'_keep_genes.rds'))
saveRDS(meta,file=paste0('checkpoints/',dataset,'_metadata.rds'))
