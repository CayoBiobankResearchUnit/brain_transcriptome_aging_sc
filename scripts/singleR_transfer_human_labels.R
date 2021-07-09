#!/usr/bin/env Rscript

options(future.globals.maxSize = 100000 * 1024^2)

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

dataset = arguments[1]

library(monocle3)
library(Seurat)
library(SingleR)
library(Matrix)

cds = readRDS(paste0('checkpoints/',dataset,'_macmul_homsap_integrated.rds'))

cds.matrix = cds@assays$integrated@data

mac.meta = colData(readRDS(paste0('checkpoints/',dataset,'_classified_manual.rds')))
hom.meta = readRDS('checkpoints/homsap_metadata.rds')

# Ignore cells that have no reference labels
hom.meta = subset(hom.meta,class_label != '')

mac = cds.matrix[,rownames(mac.meta)]
hom = cds.matrix[,rownames(hom.meta)]

macmul.integrated = Matrix(mac,sparse=TRUE)
homsap.integrated = Matrix(hom,sparse=TRUE)

macmul.singler.integrated = SingleR(
	test=macmul.integrated,
	ref=homsap.integrated,
	labels=hom.meta$cell_type_accession_label,
	method='cluster',
	clusters=mac.meta$assigned_cell_cluster)
saveRDS(macmul.singler.integrated,file=paste0('checkpoints/',dataset,'_homsap_singleR_integrated.rds'))

macmul.integrated.results = data.frame(mmul.cluster = rownames(macmul.singler.integrated), hsap.cluster=macmul.singler.integrated$labels,stringsAsFactors=FALSE)

macmul.join = unique(subset(mac.meta,select=c('assigned_cell_cluster','assigned_cell_type')))
homsap.class = unique(subset(hom.meta,select=c('cell_type_accession_label','class_label','subclass_label','cell_type_alias_label')))


# Go through macmul.singler.integrated and quantify the strenth of the assignment
macmul.integrated.hits = do.call(rbind,lapply(1:nrow(macmul.singler.integrated),function(i) {
 	x = macmul.singler.integrated[i,]
 	best.score = x$scores[which.max(x$scores)]
 	these.scores = data.frame(
 		cell_type_accession_label = colnames(x$scores)[order(x$scores,decreasing=TRUE)],
 		score = x$scores[order(x$scores,decreasing=TRUE)]
 	)
 	these.scores = merge(these.scores,homsap.class,by='cell_type_accession_label')
 	these.scores = these.scores[order(these.scores$score,decreasing=TRUE),]
 	class.diff = with(these.scores,max(score) - max(score[class_label != class_label[1]]))
 	subclass.diff = with(these.scores,max(score) - max(score[subclass_label != subclass_label[1]]))
 	out = data.frame(
 		subset(macmul.join,assigned_cell_cluster == rownames(x)),
 		these.scores[1,],
 		class.diff,
 		subclass.diff,
 		class2 = with(these.scores,class_label[class_label != class_label[1]][1]),
 		subclass2 = with(these.scores,subclass_label[subclass_label != subclass_label[1]][1])
 	)
 	out
 }))
 
saveRDS(macmul.integrated.hits,file=paste0('checkpoints/',dataset,'_homsap_singleR_integrated_results_matches.rds'))
