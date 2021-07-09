#!/usr/bin/env Rscript

options(future.globals.maxSize = 100000 * 1024^2)

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(Matrix)
library(Seurat)

dataset = arguments[1]

macmul.sct = readRDS(paste0('checkpoints/',dataset,'_seurat_sct.rds'))

homsap.sct = readRDS(paste0('checkpoints/homsap_seurat_sct.rds'))

cds.list = list(homsap.sct,macmul.sct)
names(cds.list) = c('human',dataset)

# Select integration features
integration.genes = SelectIntegrationFeatures(object.list = cds.list, nfeatures = length(intersect(rownames(homsap.sct),rownames(macmul.sct))))

saveRDS(integration.genes,file=paste0('checkpoints/homsap_integration_genes.rds'))

cds.features = integration.genes

message('Starting integration')

cds.list = PrepSCTIntegration(object.list = cds.list, anchor.features = cds.features, verbose = TRUE)

message('Prep SCT integration finished')

reference_dataset = which(names(cds.list) == 'human')

cds.anchors = FindIntegrationAnchors(object.list = cds.list,
					normalization.method = 'SCT',
					anchor.features = cds.features,
					reference = reference_dataset,
					verbose = TRUE)

saveRDS(cds.anchors,file=paste0('checkpoints/',dataset,'_macmul_homsap_anchors.rds'))

message('Find integration anchors finished')

cds.integrated = IntegrateData(anchorset = cds.anchors, normalization.method = 'SCT', verbose = TRUE)

message('Integrate data finished')

saveRDS(cds.integrated,file=paste0('checkpoints/',dataset,'_macmul_homsap_integrated.rds'))
