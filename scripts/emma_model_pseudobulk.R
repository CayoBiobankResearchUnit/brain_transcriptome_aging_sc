#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)
library(EMMREML)
library(parallel)
library(doParallel)
library(limma)
library(GMPR)

dataset = arguments[1]

pseudobulk = readRDS(paste0('checkpoints/',dataset,'_pseudobulk.rds'))
keep.genes = readRDS(paste0('checkpoints/',dataset,'_keep_genes.rds'))

meta = readRDS(paste0('checkpoints/',dataset,'_metadata.rds'))

kinship = readRDS(paste0('checkpoints/',dataset,'_relatedness.rds'))
k = kinship[['k']][sort(meta$sample),sort(meta$sample)]
z = kinship[['z']][sort(meta$sample),sort(meta$sample)]

meta = subset(meta,sample %in% rownames(z))

pseudobulk.keep.voom = lapply(cell.levels,function(x) {
	v = voom(pseudobulk[[x]])
	v$E[keep.genes[[x]],]
})
names(pseudobulk.keep.voom) = cell.levels
saveRDS(pseudobulk.keep.voom,paste0('checkpoints/',dataset,'_pseudobulk_normalized_voom.rds'))

pseudobulk.keep.gmpr = lapply(cell.levels,function(x) {
	m = pseudobulk[[x]]
	size.factor = suppressWarnings(GMPR(t(m),min_ct=1,intersect_no=1))
	g = t(t(m) / size.factor)
	g[keep.genes[[x]],]
#	v = voom(pseudobulk[[x]])
#	v$E[keep.genes[[x]],]
})
names(pseudobulk.keep.gmpr) = cell.levels
saveRDS(pseudobulk.keep.gmpr,paste0('checkpoints/',dataset,'_pseudobulk_normalized_gmpr.rds'))

pseudobulk.keep = pseudobulk.keep.gmpr

model.output = list()
for (i in 1:length(cell.levels)) {
	message(i)

	m = meta
	z.this = z
	e.this = pseudobulk.keep[[i]]

	design = model.matrix(as.formula(paste('~',paste(model.covariates,collapse=' + '))),data=m)

	clus = makeCluster(n.cores)
	registerDoParallel(cores=n.cores)  
	clusterExport(clus,varlist=c('e.this','k','z.this','design'),envir=environment())

	model.output[[cell.levels[i]]] = t(parApply(clus,e.this,1,function(y) {
		require(EMMREML)

		emma = emmreml(y = y,X = design,Z = z.this,K = k,varbetahat = T,varuhat = T,PEVuhat = T,test = T)

		p = emma$pvalbeta
		varb = emma$varbetahat
		b = emma$betahat

		c(b,varb,p[,'none'])
	}))

	colnames(model.output[[cell.levels[i]]])[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
	colnames(model.output[[cell.levels[i]]])[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
	colnames(model.output[[cell.levels[i]]])[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')

	stopCluster(clus)
}

# Now super-pseudobulk (combine all cells for an animal)

superbulk = Reduce(`+`,pseudobulk)

v = voom(superbulk)
e.keep = v$E[rowMeans(superbulk/colSums(superbulk) * 1e6) > tpm.cutoff,]

message('Superbulk')

m = meta
z.this = z
e.this = e.keep

design = model.matrix(as.formula(paste('~',paste(model.covariates,collapse=' + '))),data=m)

clus = makeCluster(n.cores)
registerDoParallel(cores=n.cores)  
clusterExport(clus,varlist=c('e.this','k','z.this','design'),envir=environment())

superbulk.output = t(parApply(clus,e.this,1,function(y) {
	require(EMMREML)

	emma = emmreml(y = y,X = design,Z = z.this,K = k,varbetahat = T,varuhat = T,PEVuhat = T,test = T)

	p = emma$pvalbeta
	varb = emma$varbetahat
	b = emma$betahat

	c(b,varb,p[,'none'])
}))

colnames(superbulk.output)[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
colnames(superbulk.output)[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
colnames(superbulk.output)[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')

stopCluster(clus)

saveRDS(model.output,file=paste0('checkpoints/',dataset,'_model_emma.rds'))
saveRDS(superbulk.output,file=paste0('checkpoints/',dataset,'_model_emma_superbulk.rds'))
saveRDS(superbulk,file=paste0('checkpoints/',dataset,'_superbulk.rds'))
saveRDS(pseudobulk.keep,paste0('checkpoints/',dataset,'_pseudobulk_normalized.rds'))
